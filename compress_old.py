#! /usr/bin/env python
import alignments
import re
import read
import pairedread
import junction
import bisect

class Compressor:
    aligned = None

    def compress(self, samFilename, compressedFilename, binary=False):
        ''' Compresses the alignments to 2 files, one for unspliced and one for spliced

            file_prefix: Prefix for all output file names
        '''

        # Read header
        header = ''
        with open(samFilename, 'r') as f:
            for line in f:
                if line[0] == '@':
                    header += line
                else:
                    break
        self.chromosomes = self.parseSAMHeader(header)
        self.aligned = alignments.Alignments(self.chromosomes)

        with open(samFilename, 'r') as f:
            self.parseAlignments(f)


        if binary:
            import shelve
            f = shelve.open(compressedFilename, 'n', 2)

            f['chroms'] = self.aligned.chromosomes
            f['exons'] = self.aligned.exons

            self.compressSpliced(f, binary)
            self.compressUnspliced(f, binary)

            f.close()

        else:
            with open(compressedFilename, 'w') as f:
                if len(self.aligned.spliced) > 0:
                    self.compressSpliced(f)

                f.write('/\n')

                if len(self.aligned.unspliced) > 0:
                    self.compressUnspliced(f)

            # Write exons and index information
            compressed = None
            with open(compressedFilename, 'r') as f:
                compressed = f.read()
            with open(compressedFilename, 'w') as f:
                chroms = []
                for k,v in self.aligned.chromosomes.items():
                    chroms.append(str(k)+','+str(v))
                f.write('\t'.join(chroms) + '\n')

                f.write('\t'.join([str(e) for e in self.aligned.exons]) + '\n')
                f.write('#' + '\t' + '\t'.join([str(i) for i in self.splicedIndex]) + '\n')

                for k,v in self.unsplicedIndex.items():
                    f.write('#' + str(k))
                    for i in xrange(len(v[0])):
                        if v[1][i] == 0:
                            f.write('\t' + str(v[0][i]))
                        else:
                            f.write('\t' + str(v[0][i]) + ',' + str(v[1][i]))
                    f.write('\n')

                f.write(compressed)

    def compressSpliced(self, filehandle, binary=False):
        ''' Compress the spliced alignments to a single file

            filehandle: File to write to
        '''

        # Compute coverage levels across every exon junction
        junctions = dict()

        # For each exon, offset in file of first junction containing that exon
        self.splicedIndex = [-1] * len(self.aligned.exons)

        for r in self.aligned.spliced:
            exonIds = r.exonIds

            if not r.xs == None:
                # XS is defined for this read
                key = '\t'.join([str(e) for e in exonIds]) + '\t' + r.xs + '\t' + str(r.NH)
                if not key in junctions:
                    covLength = 0
                    for e in exonIds:
                        covLength += self.aligned.exons[e+1] - self.aligned.exons[e]
                    junctions[key] = junction.Junction(exonIds, covLength)
                j = junctions[key]
            else:
                # XS not defined for this read, so see which one (+/-) is in the dict already or add + by default
                key1 = '\t'.join([str(e) for e in exonIds]) + '\t+\t' + str(r.NH)
                key2 = '\t'.join([str(e) for e in exonIds]) + '\t-\t' + str(r.NH)
                if key1 in junctions:
                    read.xs = '+'
                    j = junctions[key1]

                    key = key1
                elif key2 in junctions:
                    read.xs = '-'
                    j = junctions[key2]

                    key = key2
                else:
                    read.xs = '+'
                    covLength = 0
                    for e in exonIds:
                        covLength += self.aligned.exons[e+1] - self.aligned.exons[e]
                    junctions[key1] = junction.Junction(exonIds, covLength)
                    j = junctions[key1]

                    key = key1

            
            # update junction coverage vector in dictionary
            if (r.lenLeft == 0 and r.lenRight == 0):# or (r.startOffset+r.lenLeft > len(j.coverage)-r.endOffset-r.lenRight):
                for i in xrange(r.startOffset, len(j.coverage)-r.endOffset):
                    j.coverage[i] += 1
            else:
                for i in xrange(r.startOffset, r.startOffset+r.lenLeft):
                    j.coverage[i] += 1
                for i in xrange(len(j.coverage)-r.endOffset-r.lenRight, len(j.coverage)-r.endOffset):
                    j.coverage[i] += 1

            # update readLens
            if r.readLen in j.readLens:
                j.readLens[r.readLen] += 1
            else:
                j.readLens[r.readLen] = 1

            # update left and right read lengths
            if r.lenLeft > 0 or r.lenRight > 0:
                if r.lenLeft in j.lensLeft:
                    j.lensLeft[r.lenLeft] += 1
                else:
                    j.lensLeft[r.lenLeft] = 1

                if r.lenRight in j.lensRight:
                    j.lensRight[r.lenRight] += 1
                else:
                    j.lensRight[r.lenRight] = 1

        # Write junction information
        if not binary:
            #for key in sorted(junctions.keys(), key=lambda x: min([int(n) for n in x.split('\t')[:-2]])):
            for key in sorted(junctions.keys(), key=lambda x: [int(n) for n in x.split('\t')[:-2]]):
                # get offset in file
                offset = filehandle.tell()
                exons = [int(e) for e in key.split('\t')[:-2]]
                for e in exons:
                    if self.splicedIndex[e] == -1:
                        self.splicedIndex[e] = offset

                junc = junctions[key]

                # write junction information
                filehandle.write('>' + key + '\n')

                # write read lengths
                filehandle.write('\t'.join( [str(k)+','+str(v) for k,v in junc.readLens.items()] ) + '\n')

                # Write left and right read lengths
                filehandle.write('\t'.join([str(k)+','+str(v) for k,v in junc.lensLeft.items()]) + '\n')
                filehandle.write('\t'.join([str(k)+','+str(v) for k,v in junc.lensRight.items()]) + '\n')

                # write coverage
                self.writeRLE(junc.coverage, filehandle)
        else:
            filehandle['juncs'] = sorted(junctions.keys(), key=lambda x: min([int(n) for n in x.split('\t')[:-2]]))
            for key, junc in junctions.items():
                filehandle[key + 'reads'] = junc.readLens
                filehandle[key + 'left'] = junc.lensLeft
                filehandle[key + 'right'] = junc.lensRight
                filehandle[key + 'cov'] = self.RLE(junc.coverage)

    def compressUnspliced(self, filehandle, binary=False):
        ''' Compress the unspliced alignments as a run-length-encoded coverage vector

            filename: Name of file to compress to
        '''

        # For each NH value and for each exon, store the byte offset in the file for the given exon in the coverage vector
        self.unsplicedIndex = dict()
        
        # sort reads into exons
        readExons = dict()
        coverages = dict()
        for i in xrange(len(self.aligned.unspliced)):
            #r = self.aligned.unspliced[i]
            r = self.aligned.unspliced[len(self.aligned.unspliced) - i - 1]

            if not r.NH in readExons:
                readExons[r.NH] = []
                for n in xrange(len(self.aligned.exons)-1):
                    readExons[r.NH] += [[]]

                if r.NH == 1:
                    coverages[r.NH] = [0] * self.aligned.exons[-1]
                else:
                    coverages[r.NH] = [ [0, self.aligned.exons[-1]] ]

            j = bisect.bisect_right(self.aligned.exons, r.exons[0][0])-1
            readExons[r.NH][j] += [len(self.aligned.unspliced) - i - 1]

            # update coverage vector
            if r.NH == 1:
                
                if r.lenLeft == 0 and r.lenRight == 0:
                    for base in xrange(r.exons[0][0], r.exons[0][1]):
                        coverages[r.NH][base] += 1
                else:
                    for base in xrange(r.exons[0][0], r.exons[0][0]+r.lenLeft):
                        coverages[r.NH][base] += 1
                    for base in xrange(r.exons[0][1]-r.lenRight, r.exons[0][1]):
                        coverages[r.NH][base] += 1
            else:
                
                if (r.lenLeft == 0 and r.lenRight == 0):# or (r.exons[0][0]+r.lenLeft > r.exons[0][1]-r.lenRight):
                    coverages[r.NH] = self.updateRLE(coverages[r.NH], r.exons[0][0], r.exons[0][1]-r.exons[0][0], 1)
                else:
                    coverages[r.NH] = self.updateRLE(coverages[r.NH], r.exons[0][0], r.lenLeft, 1)
                    coverages[r.NH] = self.updateRLE(coverages[r.NH], r.exons[0][1]-r.lenRight, r.lenRight, 1)

        if binary:
            filehandle['NHs'] = [str(nh) for nh in readExons.keys()]

        for NH in readExons.keys():
            cov = coverages[NH]
            reads = readExons[NH]

            if not binary:
                filehandle.write('#' + str(NH) + '\n')

                # Write coverage vector
                if NH == 1:
                    line, offsetInLine = self.writeRLE(cov, filehandle, self.aligned.exons)
                    self.unsplicedIndex[NH] = (line, offsetInLine)
                else:
                    line = [-1] * len(self.aligned.exons)
                    offsetInLine = [0] * len(self.aligned.exons)
                    exonId = 0
                    baseId = 0
                    numExons = len(self.aligned.exons)

                    for c in cov:
                        while exonId < numExons and self.aligned.exons[exonId] < baseId + c[1]:
                            line[exonId] = filehandle.tell()
                            offsetInLine[exonId] = self.aligned.exons[exonId] - baseId
                            exonId += 1

                        if c[1] == 1:
                            filehandle.write(str(c[0]) + '\n')
                        else:
                            filehandle.write(str(c[0]) + '\t' + str(c[1]) + '\n')

                        baseId += c[1]

                    self.unsplicedIndex[NH] = (line, offsetInLine)
            else:
                key = str(NH)
                # Write coverage vector
                if NH == 1:
                    filehandle[key + 'cov'] = self.RLE(cov)
                else:
                    filehandle[key + 'cov'] = cov

            # Write lengths for each exon
            for i in xrange(len(reads)):
                if len(reads[i]) > 0:
                    length = self.aligned.exons[i+1] - self.aligned.exons[i]

                    exonStart = self.aligned.exons[i]

                    # Distribution of all read lengths
                    readLens = dict()

                    # Distribution of all gap lengths in paired-end reads
                    lensLeft = dict()
                    lensRight = dict()

                    for readId in reads[i]:
                        read = self.aligned.unspliced[readId]

                        alignment = read.exons

                        start = alignment[0][0] - exonStart
                        end = alignment[0][1] - exonStart

                        # update read lengths distribution
                        length = end - start
                        if length in readLens:
                            readLens[length] += 1
                        else:
                            readLens[length] = 1

                        # update left and right read lengths
                        if read.lenLeft > 0 or read.lenRight > 0:
                            #if NH > 1:
                            #    print str(NH) + ', ' + str(self.exons[i]) + '-' + str(self.exons[i+1])
                            if read.lenLeft in lensLeft:
                                lensLeft[read.lenLeft] += 1
                            else:
                                lensLeft[read.lenLeft] = 1

                            if read.lenRight in lensRight:
                                lensRight[read.lenRight] += 1
                            else:
                                lensRight[read.lenRight] = 1

                    # Write bounds
                    if not binary:
                        filehandle.write('>' + str(i) + '\n')

                        # Write read lengths
                        filehandle.write('\t'.join([str(k)+','+str(v) for k,v in readLens.items()]) + '\n')

                        # Write left and right read lengths
                        filehandle.write('\t'.join([str(k)+','+str(v) for k,v in lensLeft.items()]) + '\n')
                        filehandle.write('\t'.join([str(k)+','+str(v) for k,v in lensRight.items()]) + '\n')
                    else:
                        newKey = key + str(i)
                        filehandle[newKey + 'reads'] = readLens
                        filehandle[newKey + 'left'] = lensLeft
                        filehandle[newKey + 'right'] = lensRight

    def parseAlignments(self, filehandle):
        ''' Parse a file in SAM format
            Add the exonic region indices for each read to the correct ChromosomeAlignments object
            
            filehandle: SAM filehandle containing aligned reads
        '''

        # paired reads that have not yet found their partner
        unpaired = dict()

        for line in filehandle:
            row = line.strip().split('\t')
            if len(row) < 5:
                continue

            chromosome = str(row[2])
            
            if not row[2] in self.chromosomes.keys():
                print 'Chromosome ' + str(row[2]) + ' not found!'
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

            if not row[6] == '*':
                if row[6] == '=':
                    pair_chrom = chromosome
                else:
                    pair_chrom = row[6]
                pair_index = int(row[7])
                self.aligned.processRead(read.Read(chromosome, exons, xs, NH), pair_chrom, pair_index)
            else:
                self.aligned.processRead(read.Read(chromosome, exons, xs, NH))

        self.aligned.finalizeExons()
        self.aligned.finalizeReads()

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
                if length == 1:
                    rle.append([val])
                else:
                    rle.append([val, length])
                val = v
                length = 1

        if length == 1:
            rle.append([val])
        else:
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

    def writeRLE(self, vector, filehandle, exons=None):
        '''
            exons: If included, create an index storing the position in the file for each exon
        '''

        index = False
        numExons = 0
        if not exons == None:
            index = True
            numExons = len(exons)

            line = [-1]*numExons
            offsetInLine = [0]*numExons
            exonId = 0
            baseId = 0

        val = vector[0]
        length = 0

        for v in vector:
            if v == val:
                length += 1
            else:
                if index:
                    while exonId < numExons and exons[exonId] < baseId:
                        line[exonId] = filehandle.tell()
                        offsetInLine[exonId] = exons[exonId] - baseId + length
                        if offsetInLine[exonId] < 0:
                            print '%d, %d, %d' % (exons[exonId], baseId, length)
                            exit()
                        exonId += 1

                if length == 1:
                    filehandle.write(str(val) + '\n')
                else:
                    filehandle.write(str(val) + '\t' + str(length) + '\n')
                val = v
                length = 1

            if index:
                baseId += 1

        if index:
            while exonId < numExons and exons[exonId] < baseId:
                line[exonId] = filehandle.tell()
                offsetInLine[exonId] = exons[exonId] - baseId + length
                if offsetInLine[exonId] < 0:
                    print '%d, %d, %d' % (exons[exonId], baseId, length)
                    exit()
                exonId += 1
        if length == 1:
            filehandle.write(str(val) + '\n')
        else:
            filehandle.write(str(val) + '\t' + str(length) + '\n')

        if index:
            return line, offsetInLine
        