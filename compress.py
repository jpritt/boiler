#! /usr/bin/env python3
import alignments
import re
import read
import pairedread
import junction
import bisect
import binaryIO
import math
import readNode
import time
import os

class Compressor:
    aligned = None

    sectionLen = 100000
    exonChunkSize = 100
    junctionChunkSize = 50

    # 0 - zlib
    # 1 - lzma
    # 2 - bz2
    compressMethod = 0

    def __init__(self):   
        if self.compressMethod == 0:
            self.zlib = __import__('zlib')
        elif self.compressMethod == 1:
            self.lzma = __import__('lzma')
        elif self.compressMethod == 2:
            self.bz2 = __import__('bz2')

    def compress(self, samFilename, compressedFilename, min_filename=None, binary=False, debug=False):
        ''' Compresses the alignments to 2 files, one for unspliced and one for spliced

            file_prefix: Prefix for all output file names
        '''

        self.debug = debug

        # Read header
        header = ''
        with open(samFilename, 'r') as f:
            for line in f:
                if line[0] == '@':
                    header += line
                else:
                    break
        self.chromosomes = self.parseSAMHeader(header)
        self.aligned = alignments.Alignments(self.chromosomes, self.debug)

        self.compressByCluster(samFilename, compressedFilename, min_filename)

    def add_to_partition(self, r, max_len, partitions, debug=False):
        '''

        :param r: Read to add
        :param max_len: Maximum fragment length seen so far
        :param partitions: List of previously computed partitions
        :return: Key where r is stored, updated partitions list, updated max_len
        '''

        exonIds = r.exonIds

        if not r.xs == None:
            # XS is defined for this read
            key = '\t'.join([str(e) for e in exonIds]) + '\t' + r.xs + '\t' + str(r.NH)

            if not key in partitions:
                covLength = 0
                #boundaries = []
                for i in range(len(exonIds)):
                    e = exonIds[i]
                    subexon_length = self.aligned.exons[e+1] - self.aligned.exons[e]

                    covLength += subexon_length
                    #if i == 0:
                    #    boundaries.append(subexon_length)
                    #elif i < len(exonIds)-1:
                    #    boundaries.append(boundaries[-1] + subexon_length)

                partitions[key] = junction.Junction(exonIds, covLength)
                partitions[key].xs = r.xs
                partitions[key].NH = r.NH

                partitions[key].countBefore = 0

                if debug:
                    print(key)
                    print(partitions[key].coverage)
                    print('')
            j = partitions[key]
        else:
            # XS not defined for this read, so see which one (+/-) is in the dict already or add + by default
            key1 = '\t'.join([str(e) for e in exonIds]) + '\t+\t' + str(r.NH)
            key2 = '\t'.join([str(e) for e in exonIds]) + '\t-\t' + str(r.NH)
            if key1 in partitions:
                read.xs = '+'
                j = partitions[key1]

                key = key1
            elif key2 in partitions:
                read.xs = '-'
                j = partitions[key2]

                key = key2
            else:
                read.xs = '+'
                covLength = 0
                #boundaries = []

                for i in range(len(exonIds)):
                    e = exonIds[i]
                    subexon_length = self.aligned.exons[e+1] - self.aligned.exons[e]

                    covLength += subexon_length
                    #if i == 0:
                    #    boundaries.append(subexon_length)
                    #elif i < len(exonIds)-1:
                    #    boundaries.append(boundaries[-1] + subexon_length)

                partitions[key1] = junction.Junction(exonIds, covLength)
                partitions[key1].xs = '+'
                partitions[key1].NH = r.NH

                partitions[key1].countBefore = 0

                j = partitions[key1]

                if debug:
                    print(key1)
                    print(partitions[key1].coverage)
                    print('')

                key = key1


        # update junction coverage vector in dictionary
        if (r.lenLeft == 0 and r.lenRight == 0):
            j.coverage = self.updateRLE(j.coverage, r.startOffset, j.length-r.endOffset-r.startOffset, 1)
        else:
            j.coverage = self.updateRLE(j.coverage, r.startOffset, r.lenLeft, 1)
            j.coverage = self.updateRLE(j.coverage, j.length-r.endOffset-r.lenRight, r.lenRight, 1)

        # If this is a paired read, lenLeft and lenRight are both > 0
        if r.lenLeft > 0 or r.lenRight > 0:
            if not r.startOffset+r.readLen+r.endOffset == j.length:
                print('Error in read/junction length!')
                exit()

            # update pairedLens
            if r.readLen in j.pairedLens:
                j.pairedLens[r.readLen] += 1
            else:
                j.pairedLens[r.readLen] = 1

            if r.lenLeft in j.lensLeft:
                j.lensLeft[r.lenLeft] += 1
            else:
                j.lensLeft[r.lenLeft] = 1

            if r.lenRight in j.lensRight:
                j.lensRight[r.lenRight] += 1
            else:
                j.lensRight[r.lenRight] = 1

            # Update maximum fragment length
            if r.lenLeft > max_len:
                max_len = r.lenLeft
            if r.lenRight > max_len:
                max_len = r.lenRight
        else:
            # update unpairedLens
            if r.readLen in j.unpairedLens:
                j.unpairedLens[r.readLen] += 1
            else:
                j.unpairedLens[r.readLen] = 1

            if r.readLen > max_len:
                max_len = r.readLen

        return key, partitions, max_len


    def computeJunctions(self, debug=False):

        # Compute coverage levels across every exon junction
        partitions = dict()

        # Maximum fragment length
        max_len = 0

        for r in self.aligned.spliced:
            _, partitions, max_len = self.add_to_partition(r, max_len, partitions)

        # Junction index: list of junctions and length of each chunk
        self.sortedJuncs = sorted(partitions.keys(), key=lambda x: [int(n) for n in x.split('\t')[:-2]])
        self.junctionChunkLens = [0] * int(math.ceil(len(self.sortedJuncs) / self.junctionChunkSize))

        return partitions, max_len


    def compressSpliced(self, junctions, maxReadLen, filehandle, binary=False, debug=False):
        ''' Compress the spliced alignments to a single file

            filehandle: File to write to
        '''

        # Determine the number of bytes for read lengths
        readLenBytes = binaryIO.findNumBytes(maxReadLen)
        s = binaryIO.valToBinary(1, readLenBytes)

        s += binaryIO.writeJunctionsList(self.sortedJuncs, 2)

        for key in self.sortedJuncs:
            junc = junctions[key]

            if binary:
                s += binaryIO.writeJunction(readLenBytes, junc)
            else:
                print('Non-binary spliced encoding is not supported')
                exit()

        # Write to file
        start = filehandle.tell()
        filehandle.write(self.compressString(s))

        # return length of chunk in file
        return filehandle.tell()-start


    def compressUnspliced(self, filehandle, binary=False, debug=False):
        ''' Compress the unspliced alignments as a run-length-encoded coverage vector

            filename: Name of file to compress to
        '''

        st = time.time()

        maxFragLen = 0
        maxReadLen = 0

        # divide reads by NH value
        reads_by_NH = dict()
        coverages = dict()
        fragment_coverages = dict()

        range_start = self.aligned.exons[0]
        range_end = self.aligned.exons[-1]

        for i in range(len(self.aligned.unspliced)):
            r = self.aligned.unspliced[len(self.aligned.unspliced) - i - 1]

            if binary:
                # find maximum fragment and read lengths
                if r.lenLeft == 0 and r.lenRight == 0:
                    if r.readLen > maxReadLen:
                        maxReadLen = r.readLen
                else:
                    if r.lenLeft > maxReadLen:
                        maxReadLen = r.lenLeft
                    if r.lenRight > maxReadLen:
                        maxReadLen = r.lenRight
                    if r.readLen > maxFragLen:
                        maxFragLen = r.readLen

            if not r.NH in reads_by_NH:
                reads_by_NH[r.NH] = []
                for n in range(len(self.aligned.exons)-1):
                    reads_by_NH[r.NH] += [[]]

                if r.NH == 1:
                    coverages[r.NH] = [0] * (range_end - range_start)
                    #fragment_coverages[r.NH] = [0] * self.aligned.exons[-1]
                else:
                    coverages[r.NH] = [ [0, range_end-range_start] ]
                    #fragment_coverages[r.NH] = [ [0, self.aligned.exons[-1]] ]

            j = bisect.bisect_right(self.aligned.exons, r.exons[0][0])-1
            reads_by_NH[r.NH][j] += [len(self.aligned.unspliced) - i - 1]

            if self.debug:
                if len(r.exons) > 1:
                    print('Error! %d exons' % len(r.exons))
                if not r.readLen == (r.exons[0][1]-r.exons[0][0]):
                    print('Lengths %d and %d do not match!' % (r.readLen, r.exons[0][1]-r.exons[0][0]))

            if r.NH == 1:
                #for n in range(r.exons[0][0], r.exons[0][1]):
                #    fragment_coverages[r.NH][n] += 1

                # update coverage vector
                if r.lenLeft == 0 and r.lenRight == 0:
                    for n in range(r.exons[0][0]-range_start, r.exons[0][1]-range_start):
                        coverages[r.NH][n] += 1
                        #fragment_coverages[r.NH][n] += 1
                else:
                    for n in range(r.exons[0][0]-range_start, r.exons[0][0]-range_start+r.lenLeft):
                        coverages[r.NH][n] += 1
                    for n in range(r.exons[0][1]-range_start-r.lenRight, r.exons[0][1]-range_start):
                        coverages[r.NH][n] += 1
            else:
                #fragment_coverages[r.NH] = self.updateRLE(fragment_coverages[r.NH], r.exons[0][0], r.exons[0][1]-r.exons[0][0], 1)

                # update coverage vector
                if r.lenLeft == 0 and r.lenRight == 0:
                    coverages[r.NH] = self.updateRLE(coverages[r.NH], r.exons[0][0]-range_start, r.exons[0][1]-r.exons[0][0], 1)
                else:
                    coverages[r.NH] = self.updateRLE(coverages[r.NH], r.exons[0][0]-range_start, r.lenLeft, 1)
                    coverages[r.NH] = self.updateRLE(coverages[r.NH], r.exons[0][1]-range_start-r.lenRight, r.lenRight, 1)


        if 1 in coverages:
            coverages[1] = self.RLE(coverages[1])
            #fragment_coverages[1] = self.RLE(fragment_coverages[1])

        if binary:
            # find length in bytes to fit fragment and read lengths
            fragLenBytes = binaryIO.findNumBytes(maxFragLen)
            readLenBytes = binaryIO.findNumBytes(maxReadLen)

            length = self.writeUnsplicedBinary(filehandle, reads_by_NH, coverages, fragment_coverages, fragLenBytes, readLenBytes, debug)
            return length
        else:
            print('Non-binary unspliced compression not supported')
            exit()
            self.writeUnspliced(filehandle, reads_by_NH, coverages)


    def writeUnsplicedBinary(self, filehandle, reads_by_NH, coverages, fragment_coverages, fragLenBytes, readLenBytes, debug=False):
        ''' Compress the unspliced alignments as a run-length-encoded coverage vector

            filename: Name of file to compress to
        '''

        s = binaryIO.valToBinary(2, len(coverages))
        s += binaryIO.valToBinary(1, fragLenBytes)
        s += binaryIO.valToBinary(1, readLenBytes)

        if debug:
            print('Exons:')
            print(self.aligned.exons)
            print('')

        for NH in sorted(reads_by_NH.keys()):
            if debug:
                print('NH: ' + str(NH))
            s += binaryIO.valToBinary(2, NH)

            #cov = fragment_coverages[NH]
            cov = coverages[NH]
            reads = reads_by_NH[NH]

            if debug:
                print('Reads:')
                print(reads)

            s += binaryIO.writeCov(cov)


            # Write read length distributions
            for i in range(len(reads)):
                if i == len(reads)-1 and debug:
                    all_r = []
                    for id in reads[i]:
                        r = self.aligned.unspliced[id]
                        if r.lenLeft == 0 and r.lenRight == 0:
                            all_r.append(r.exons[0])
                        else:
                            all_r.append([r.exons[0][0], r.exons[0][0]+r.lenLeft])
                            all_r.append([r.exons[0][1]-r.lenRight, r.exons[0][1]])
                    print('Reads:')
                    print(all_r) 

                    print('Fragments:')
                    print([self.aligned.unspliced[id].exons for id in reads[i]])
                    cov_start = self.aligned.exons[0]
                    print('Coverage:')
                    c = self.RLEtoVector(coverages[NH])
                    print(c[self.aligned.exons[i]-cov_start:self.aligned.exons[i+1]-cov_start])

                start = self.aligned.exons[i]

                # Distribution of all read lengths
                unpairedLens = dict()

                pairedLens = dict()

                # Distribution of all gap lengths in paired-end reads
                lensLeft = dict()
                lensRight = dict()

                offset = self.aligned.exons[i]
                #unpaired = []
                #paired = []
                for readId in reads[i]:
                    read = self.aligned.unspliced[readId]

                    # update left and right read lengths
                    if read.lenLeft > 0 or read.lenRight > 0:
                        #paired.append([[read.exons[0][0]-start, read.exons[0][0]+read.lenLeft-start], [read.exons[-1][1]-read.lenRight-start, read.exons[-1][1]-start]])

                        # update paired read lengths distribution
                        length = read.readLen
                        if length in pairedLens:
                            pairedLens[length] += 1
                        else:
                            pairedLens[length] = 1

                        if read.lenLeft in lensLeft:
                            lensLeft[read.lenLeft] += 1
                        else:
                            lensLeft[read.lenLeft] = 1

                        if read.lenRight in lensRight:
                            lensRight[read.lenRight] += 1
                        else:
                            lensRight[read.lenRight] = 1
                    else:
                        #unpaired.append([read.exons[0][0]-start, read.exons[-1][1]-start])

                        # update unpaired read lengths distribution
                        length = read.readLen
                        if length in unpairedLens:
                            unpairedLens[length] += 1
                        else:
                            unpairedLens[length] = 1
                
                if i == len(reads)-1 and debug:
                    print(unpairedLens)
                    print(pairedLens)
                    print(lensLeft)
                    print(lensRight)
                    print('')

                s += binaryIO.writeLens(readLenBytes, unpairedLens)
                s += binaryIO.writeLens(fragLenBytes, pairedLens)

                if len(pairedLens) > 0:
                    s += binaryIO.writeLens(readLenBytes, lensLeft)
                    if len(lensLeft) > 0:
                        s += binaryIO.writeLens(readLenBytes, lensRight)

        startPos = filehandle.tell()
        filehandle.write(self.compressString(s))
        return filehandle.tell() - startPos


    def compressByCluster(self, input_name, compressed_name, intermediate_name=None):
        '''
        Read a sorted SAM file and compress in segments determined by clusters of reads

        :param filename:
        :return:
        '''

        readLens = dict()

        # If coverage is 0 for at least this many bases end of a potential gene
        overlapRadius = 50

        if self.debug:
            start_time = time.time()

        unspliced_index = []
        spliced_index = []
        clusters = []

        first = True

        with open(input_name, 'r') as filehandle:
            for line in filehandle:
                row = line.strip().split('\t')

                # Check if header line
                if len(row) < 6:
                    continue

                if not row[2] in self.chromosomes.keys():
                    print('Error! Chromosome ' + str(row[2]) + ' not found!')
                    exit()

                start = self.aligned.chromOffsets[row[2]] + int(row[3])

                if self.aligned.gene_bounds and start > (self.aligned.gene_bounds[-1] + overlapRadius):
                    # Compress most recent cluster
                    self.aligned.finalizeUnmatched()
                    self.aligned.finalizeExons()

                    clusters.append(self.aligned.exons)


                    if intermediate_name:
                        if first:
                            with open(intermediate_name, 'w') as f1:
                                self.aligned.writeSAM(f1, True)
                            first = False
                        else:
                            with open(intermediate_name, 'a') as f1:
                                self.aligned.writeSAM(f1, False)

                    #self.aligned.finalizeReads()

                    #junctions, maxReadLen = self.computeJunctions(debug=False)
                    junctions, maxReadLen = self.aligned.computeJunctions()
                    self.sortedJuncs = sorted(junctions.keys(), key=lambda x: [int(n) for n in x.split('\t')[:-2]])

                    #for j in junctions:
                    #    print(j)
                    #    print(junctions[j].pairedLens)
                    #    print(junctions[j].unpairedLens)
                    #    print('')
                    #exit()

                    with open('temp.bin', 'ab') as f:
                        #l = self.compressUnspliced(f, binary=True, debug=False)
                        #unspliced_index.append(l)

                        l = self.compressSpliced(junctions, maxReadLen, f, binary=True, debug=False)
                        spliced_index.append(l)

                    # Start new cluster
                    self.aligned.resetCluster()
                    self.aligned.exons.add(start)


                exons = self.parseCigar(row[5], int(row[3]))

                # find XS and NH values
                xs = None
                NH = 1
                for r in row[11 : len(row)]:
                    if r[0:5] == 'XS:A:' or r[0:5] == 'XS:a:':
                        xs = r[5]
                    elif r[0:3] == 'NH:':
                        NH = int(r[5:])


                r = read.Read(row[2], exons, xs, NH)
                if row[6] == '=':
                    r.pairOffset = int(row[7])
                    self.aligned.processRead(r, row[0], paired=True)
                else:
                    self.aligned.processRead(r, row[0], paired=False)

            # Compress final cluster

            self.aligned.finalizeUnmatched()
            self.aligned.finalizeExons()

            clusters.append(self.aligned.exons)

            if intermediate_name:
                if first:
                    with open(intermediate_name, 'w') as f1:
                        self.aligned.writeSAM(f1, True)
                    first = False
                else:
                    with open(intermediate_name, 'a') as f1:
                        self.aligned.writeSAM(f1, False)

            #self.aligned.finalizeReads()

            #junctions, maxReadLen = self.computeJunctions(debug=False)
            junctions, maxReadLen = self.aligned.computeJunctions()
            self.sortedJuncs = sorted(junctions.keys(), key=lambda x: [int(n) for n in x.split('\t')[:-2]])

            with open('temp.bin', 'ab') as f:
                #l = self.compressUnspliced(f, binary=True, debug=False)
                #unspliced_index.append(l)

                l = self.compressSpliced(junctions, maxReadLen, f, binary=True, debug=False)
                spliced_index.append(l)


        # Write index information and append spliced and unspliced files
        with open(compressed_name, 'wb') as f:
            s = binaryIO.writeChroms(self.aligned.chromosomes)
            s += binaryIO.writeClusters(clusters)
            #s += binaryIO.writeList(unspliced_index)
            s += binaryIO.writeList(spliced_index)
            f.write(s)

            with open('temp.bin', 'rb') as f2:
                f.write(f2.read())

        os.remove('temp.bin')


        if self.debug:
            end_time = time.time()
            print('Compression time: %0.3fs' % (end_time - start_time))


    def parseAlignments(self, filename):
        ''' Parse a file in SAM format
            Add the exonic region indices for each read to the correct ChromosomeAlignments object
            
            filehandle: SAM filehandle containing aligned reads
        '''

        if self.debug:
            start = time.time()

        with open(filename, 'r') as filehandle:
            for line in filehandle:
                row = line.strip().split('\t')

                # Check if header line
                if len(row) < 5:
                    continue

                if not row[2] in self.chromosomes.keys():
                    print('Chromosome ' + str(row[2]) + ' not found!')
                    continue

                exons = self.parseCigar(row[5], int(row[3]))

                # find XS value:
                xs = None
                NH = 1
                for r in row[11 : len(row)]:
                    if r[0:5] == 'XS:A:' or r[0:5] == 'XS:a:':
                        xs = r[5]
                    elif r[0:3] == 'NH:':
                        NH = int(r[5:])

                
                r = read.Read(row[2], exons, xs, NH)
                if row[6] == '=':
                    r.pairOffset = int(row[7])
                    self.aligned.processRead(r, row[0], paired=True)
                else:
                    self.aligned.processRead(r, row[0], paired=False)

        if self.debug:
            t = time.time() - start
            print('%0.3fs' % t)
            start = time.time()
            print('Finalizing unmatched')

        self.aligned.finalizeUnmatched()

        if self.debug:
            t = time.time() - start
            print('%0.3fs' % t)
            start = time.time()
            print('Finalizing exons')

        self.aligned.finalizeExons()

        if self.debug:
            t = time.time() - start
            print('%0.3fs' % t)

    def parseCigar(self, cigar, offset):
        ''' Parse the cigar string starting at the given index of the genome
            Returns a list of offsets for each exonic region of the read [(start1, end1), (start2, end2), ...]
        '''

        exons = []
        newExon = True

        # Parse cigar string
        match = re.search("\D", cigar)
        while match:
            index = match.start()
            length = int(''.join(cigar[:index]))

            if cigar[index] == 'N':
                # Separates contiguous exons, so set boolean to start a new one
                newExon = True
            elif cigar[index] == 'M':
                # If in the middle of a contiguous exon, append the length to it, otherwise start a new exon
                if newExon:
                    exons.append([offset, offset+length])
                    newExon = False
                else:
                    exons[-1][1] += length
            elif cigar[index] == 'D':
                # If in the middle of a contiguous exon, append the deleted length to it
                if not newExon:
                    exons[-1][1] += length

            offset += length
            cigar = cigar[index+1:]
            match = re.search("\D", cigar)

        return exons

    def parseSAMHeader(self, header):
        chromosomes = dict()
        for line in header.split('\n'):
            if line[0:3] == '@SQ':
                row = line.strip().split('\t')
                chromosomes[row[1][3:]] = int(row[2][3:])
        return chromosomes

    def RLE(self, vector):
        rle = []

        val = vector[0]
        length = 0
        for v in vector:
            if v == val:
                length += 1
            else:
                rle.append([val, length])
                val = v
                length = 1

        rle.append([val, length])

        return rle

    def updateRLE(self, RLE, start, length, value):
        i = 0

        while start > 0 and start >= RLE[i][1]:
            start -= RLE[i][1]
            i += 1

        if start > 0:
            RLE = RLE[:i] + [ [RLE[i][0], start], [RLE[i][0], RLE[i][1]-start] ] + RLE[i+1:]
            i += 1

        while length > 0 and length >= RLE[i][1]:
            RLE[i][0] += value
            if RLE[i][0] < 0:
                RLE[i][0] = 0

            length -= RLE[i][1]
            #if i > 0 and RLE[i][0] == RLE[i-1][0]:
            #    RLE = RLE[:i-1] + [ [RLE[i][0], RLE[i-1][1]+RLE[i][1]] ] + RLE[i+1:]
            #else:
            i += 1

        if length > 0:
            RLE = RLE[:i] + [ [max(RLE[i][0]+value,0), length], [RLE[i][0], RLE[i][1]-length] ] + RLE[i+1:]
        #elif i > 0 and RLE[i][0]  == RLE[i-1][0]:
        #    RLE = RLE[:i-1] + [ [RLE[i][0], RLE[i-1][1]+RLE[i][1]] ] + RLE[i+1:]

        return RLE

    def writeIndexBinary(self, f):
        s = binaryIO.valToBinary(4, self.exonChunkSize)
        s += binaryIO.valToBinary(4, self.junctionChunkSize)

        # Write junction names
        s += binaryIO.writeJunctionsList(self.sortedJuncs, self.exonBytes)

        # Write length of junction chunks
        s += binaryIO.writeList(self.junctionChunkLens)

        # Number of NH values
        s += binaryIO.valToBinary(4, len(self.unsplicedIndex.items()))
        for k in self.unsplicedIndex.keys():
            # NH value
            s += binaryIO.valToBinary(4, k)

            # Write unspliced coverage index
            s += binaryIO.writeList(self.unsplicedIndex[k])
        
            # Write unspliced exon reads index
            s += binaryIO.writeList(self.unsplicedExonsIndex[k])

        # Compress and write to file
        s = self.compressString(s)
        binaryIO.writeVal(f, 4, len(s))
        f.write(s)

    def compressString(self, s):
        ''' Use a predefined python library to compress the given string.
            Return the compressed string '''

        if self.compressMethod == 0:
            return self.zlib.compress(s)
        elif self.compressMethod == 1:
            return self.lzma.compress(s)
        elif self.compressMethod == 2:
            return self.bz2.compress(s)

    def RLEtoVector(self, rle):
        ''' Convert a run length encoded vector to the original coverage vector '''

        vector = []
        for row in rle:
            vector += [row[0]] * row[1]
        return vector














    '''
    TODO
    These functions are out of date and need to be revised/rewritten
    '''

    def writeUnspliced(self, filehandle, readExons, coverages):
        breakpoints = range(0, self.aligned.exons[-1], self.sectionLen)

        # For each NH value and for each exon, store the byte offset in the file for the given exon in the coverage vector
        self.unsplicedIndex = dict()

        for NH in sorted(readExons.keys()):
            cov = coverages[NH]
            reads = readExons[NH]

            filehandle.write('#' + str(NH) + '\n')

            offsets = [0] * len(breakpoints)

            pos = 0
            currBreakpoint = 0
            segmentRLE = []
            for c in cov:
                if breakpoints[currBreakpoint] == pos:
                    offsets[currBreakpoint] = filehandle.tell()
                    currBreakpoint += 1

                segmentEnd = pos+c[1]
                while currBreakpoint < len(breakpoints) and breakpoints[currBreakpoint] < segmentEnd:
                    length = breakpoints[currBreakpoint] - pos
                    if length == 1:
                        filehandle.write(str(c[0]) + '\n')
                    else:
                        filehandle.write(str(c[0]) + '\t' + str(length) + '\n')

                    c[1] -= length
                    pos += length
                    offsets[currBreakpoint] = filehandle.tell()
                    currBreakpoint += 1

                if c[1] == 1:
                    filehandle.write(str(c[0]) + '\n')
                else:
                    filehandle.write(str(c[0]) + '\t' + str(c[1]) + '\n')
                pos += c[1]

            self.unsplicedIndex[NH] = offsets

            # Write lengths for each exon
            for i in range(len(reads)):
                if len(reads[i]) > 0:
                    length = self.aligned.exons[i+1] - self.aligned.exons[i]

                    exonStart = self.aligned.exons[i]

                    # Distribution of all read lengths
                    unpairedLens = dict()
                    pairedLens = dict()

                    # Distribution of all gap lengths in paired-end reads
                    lensLeft = dict()
                    lensRight = dict()

                    for readId in reads[i]:
                        read = self.aligned.unspliced[readId]
                        alignment = read.exons

                        # update left and right read lengths
                        if read.lenLeft > 0 or read.lenRight > 0:
                            # update read lengths distribution
                            length = alignment[0][1] - alignment[0][0]
                            if length in pairedLens:
                                pairedLens[length] += 1
                            else:
                                pairedLens[length] = 1

                            if read.lenLeft in lensLeft:
                                lensLeft[read.lenLeft] += 1
                            else:
                                lensLeft[read.lenLeft] = 1

                            if read.lenRight in lensRight:
                                lensRight[read.lenRight] += 1
                            else:
                                lensRight[read.lenRight] = 1
                        else:
                            # update read lengths distribution
                            length = alignment[0][1] - alignment[0][0]
                            if length in unpairedLens:
                                unpairedLens[length] += 1
                            else:
                                unpairedLens[length] = 1

                    # Write bounds
                    filehandle.write('>' + str(i) + '\n')

                    # Write read lengths
                    filehandle.write('\t'.join([str(k)+','+str(v) for k,v in readLens.items()]) + '\n')

                    # Write left and right read lengths
                    filehandle.write('\t'.join([str(k)+','+str(v) for k,v in lensLeft.items()]) + '\n')
                    filehandle.write('\t'.join([str(k)+','+str(v) for k,v in lensRight.items()]) + '\n')

    def writeRLE(self, vector, filehandle):
        ''' Write vector that is already run-length encoded '''
        for v in vector:
            if v[1] == 1:
                filehandle.write(str(v[0]) + '\n')
            else:
                filehandle.write(str(v[0]) + '\t' + str(v[1]) + '\n')

    def writeIndex(self, f):
        # Write junctions index
        f.write(bytes('#' + str(self.exonChunkSize) + '\n', 'UTF-8'))
        f.write(bytes('#' + ' '.join([str(k)+','+str(v) for (k,v) in self.junctionIndex]) + '\n', 'UTF-8'))

        for k,v in self.unsplicedIndex.items():
            # Write unspliced coverage index
            f.write(bytes('#' + str(k) + '\t' + '\t'.join([str(i) for i in v]) + '\n', 'UTF-8'))

            # Write unspliced exon reads index
            f.write(bytes('#' + str(k) + '\t' + '\t'.join([str(i) for i in self.unsplicedExonsIndex[k]]) + '\n', 'UTF-8'))






    '''
    These functions are currently unused
    '''

    def generatePairs(self, reads):
        '''
        :param reads: A list of mixed unpaired and paired reads
        :return: A list of index pairs
        '''

        if len(reads) == 0:
            return []

        # Construct linked list of reads, sorted by start index
        root = None
        for r in reads:
            start1 = r.exons[0][0]
            node1 = readNode.ReadNode(start1)
            if root == None:
                root = node1
            else:
                root = root.addRead(node1)

            if r.lenLeft > 0 or r.lenRight > 0:
                start2 = start1 + r.readLen - r.lenLeft
                node2 = readNode.ReadNode(start2)
                root = root.addRead(node2)

                node1.pair = node2
                node2.pair = node1

        # Index reads in list
        root.index()

        # Find pairs
        pairs = []
        node = root
        while node:
            if node.pair and node.pair.id > node.id:
                pairs.append((node.id, node.pair.id))
            node = node.next
        return pairs
