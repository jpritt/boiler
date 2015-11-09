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


    def compressByCluster(self, input_name, compressed_name, intermediate_name=None):
        '''
        Read a sorted SAM file and compress in segments determined by clusters of reads

        :param filename:
        :return:
        '''

        # If coverage is 0 for at least this many bases end of a potential gene
        overlapRadius = 50

        if self.debug:
            start_time = time.time()

        spliced_index = []
        clusters = []

        first = True

        with open(input_name, 'r') as filehandle:
            for line in filehandle:
                row = line.strip().split('\t')

                # Check if header line
                if len(row) < 6 or row[5] == '*':
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

                    junctions, maxReadLen = self.aligned.computeJunctions()
                    self.sortedJuncs = sorted(junctions.keys(), key=lambda x: [int(n) for n in x.split('\t')[:-2]])

                    with open('temp.bin', 'ab') as f:
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

            junctions, maxReadLen = self.aligned.computeJunctions()
            self.sortedJuncs = sorted(junctions.keys(), key=lambda x: [int(n) for n in x.split('\t')[:-2]])

            with open('temp.bin', 'ab') as f:
                l = self.compressSpliced(junctions, maxReadLen, f, binary=True, debug=False)
                spliced_index.append(l)

        # Write index information and append spliced and unspliced files
        with open(compressed_name, 'wb') as f:
            s = binaryIO.writeChroms(self.aligned.chromosomes)
            s += binaryIO.writeClusters(clusters)
            s += binaryIO.writeList(spliced_index)
            f.write(s)

            with open('temp.bin', 'rb') as f2:
                f.write(f2.read())

        os.remove('temp.bin')


        if self.debug:
            end_time = time.time()
            print('Compression time: %0.3fs' % (end_time - start_time))

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

    def compressString(self, s):
        ''' Use a predefined python library to compress the given string.
            Return the compressed string '''

        if self.compressMethod == 0:
            return self.zlib.compress(s)
        elif self.compressMethod == 1:
            return self.lzma.compress(s)
        elif self.compressMethod == 2:
            return self.bz2.compress(s)
