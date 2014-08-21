#! /usr/bin/env python

import bisect

class Junction:
    readLens = None
    lensLeft = None
    lensRight = None
    coverage = []
    NH = 1

class Chromosome:
    # RLE coverage vector
    coverage = []
    NH = 1


class CompressedAlignment:
    '''
        Contains the data from a compressed set of aligned reads
    '''
    exons = []
    junctions = []

    def __init__(self, header, spliced, unspliced):
        
        self.chromosomes = dict()
        for line in header.split('\n'):
            if line[0:3] == '@SQ':
                row = line.strip().split('\t')
                self.chromosomes[row[1][3:]] = int(row[2][3:])

        with open(spliced, 'r') as f:
            self.exons = None
            junc = None

            self.spliced = []
            for line in f:
                if self.exons == None:
                    # First line contains exons
                    self.exons = [int(e) for e in line.rstrip().split('\t')]
                else:
                    # Remaining lines grouped by junction

                    if line[0] == '>':
                        if not junc == None:
                            self.junctions.append(junc)

                        # Start of a new junction
                        key = line[1:].rstrip().split('\t')
                        junctionExons = [int(e) for e in key[:-2]]
                        length = 0
                        for e in junctionExons:
                            length += self.exons[e+1] - self.exons[e]

                        junc = Junction()
                        junc.readLens = None
                        junc.lensLeft = None
                        junc.lensRight = None

                        #junc.xs = key[-2]
                        junc.NH = int(key[-1])

                    else:
                        if junc.readLens == None:
                            # First line after '>' contains read length distribution
                            junc.readLens = dict()
                            for readLen in line.rstrip().split('\t'):
                                readLen = readLen.split(',')
                                junc.readLens[int(readLen[0])] = int(readLen[1])
                            junc.coverage = []
                        elif junc.lensLeft == None:
                            # Second line after '>' contains left read length distribution
                            junc.lensLeft = dict()
                            if len(line.rstrip()) > 0:
                                for length in line.rstrip().split('\t'):
                                    length = length.split(',')
                                    junc.lensLeft[int(length[0])] = int(length[1])
                        elif junc.lensRight == None:
                            # Third line after '>' contains right read length distribution
                            junc.lensRight = dict()
                            if len(line.rstrip()) > 0:
                                for length in line.rstrip().split('\t'):
                                    length = length.split(',')
                                    junc.lensRight[int(length[0])] = int(length[1])
                        else:
                            # Rest of lines contain run-length-encoded coverage vector
                            row = line.rstrip().split('\t')
                            if len(row) == 1:
                                length = 1
                            else:
                                length = int(row[1])
                            val = int(row[0])

                            junc.coverage += [val] * length


            # process the final junction
            self.junctions.append(junc)

        with open(unspliced, 'r') as f:
            NH = 1
            readLens = dict()
            exonStart = 0
            exonEnd = 0
            for line in f:
                if line[0] == '#':
                    # process last chromosome

                    chrom = Chromosome()
                    chrom.NH = int(line.rstrip()[1:])

                    lineNum = 0
                    readLens = dict()
                    lensLeft = dict()
                    lensRight = dict()
                    coverage = []
                    RLE_segment = True

                elif RLE_segment and not line[0] == '>':
                    # Rest of lines contain run-length-encoded coverage vector
                    row = line.rstrip().split("\t")

                    if len(row) == 1:
                        length = 1
                    else:
                        length = int(row[1])

                    coverage += [int(row[0])] * length

                elif line[0] == '>':
                    RLE_segment = False
                    if len(readLens) > 0:
                        # Process previous exon
                        unpaired, paired = self.findReads(readLens, lensLeft, lensRight, coverage[exonStart:exonEnd])

                        for r in unpaired:
                            self.unspliced.append(read.Read(self.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))

                        for p in paired:
                            self.paired.append(pairedread.PairedRead(self.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                                     self.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH ))
                    #exonId += 1
                    lineNum = 0
                    bounds = line[1:].rstrip().split('-')
                    exonStart = int(bounds[0])
                    exonEnd = int(bounds[1])

                elif lineNum == 1:
                    # First line after '>' contains read lengths
                    readLens = dict()
                    for row in line.rstrip().split('\t'):
                        if len(row) > 0:
                            readLen = row.split(',')
                            readLens[int(readLen[0])] = int(readLen[1])
                elif lineNum == 2:
                    # Second line after '>' contains left read length distribution
                    lensLeft = dict()
                    if len(line.rstrip()) > 0:
                        for length in line.rstrip().split('\t'):
                            length = length.split(',')
                            lensLeft[int(length[0])] = int(length[1])
                elif lineNum == 3:
                    # Third line after '>' contains right read length distribution
                    lensRight = dict()
                    if len(line.rstrip()) > 0:
                        for length in line.rstrip().split('\t'):
                            length = length.split(',')
                            lensRight[int(length[0])] = int(length[1])
                lineNum += 1

            # Process the final exon
            if len(readLens) > 0:
                unpaired, paired = self.findReads(readLens, lensLeft, lensRight, coverage[exonStart:exonEnd])

                for r in unpaired:
                    self.unspliced.append(read.Read(self.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))

                for p in paired:
                    self.paired.append(pairedread.PairedRead(self.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                             self.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH))

    def getChromosome(self, index):
        ''' Return chromosome name containing the given index from the whole-genome vector
        '''

        for c in self.chromosomes:
            if self.chromosomes[c] > index:
                return c
            else:
                index -= self.chromosomes[c]

    def getCoverage(chrom, start=None, end=None):
        ''' Return coverage vector over the given interval
        '''
        if not chrom in self.chromosomes:
            print 'Error! Chromosome name not recognized!'
            print 'Chromosomes: ' + ', '.join(self.chromosomes.keys())
            return []

        coverage = []
        if start == None or end == None:
            start = 0
            end = len(self.chromosomes[chrom])
            coverage = self.chromosomes[chrom]
        else:
            coverage = self.chromosomes[chrom][start:end]

        # add coverage from junctions
        startExon = 0
        endExon = len(self.exons)
        while start + self.chromosomes[chrom] > self.exons[startExon+1]:
            startExon += 1
        while end + self.chromosomes[chrom] <= self.exons[endExon-1]:
            endExon -= 1

        for junc in self.junctions:
            overlapping = False
            for e in junc.exons:
                if e >= startExon and e < endExon:
                    overlapping = True
                    break

            if overlapping:
                indexInCov = 0
                for e in junc.exons:
                    if e >= startExon and e < endExon:
                        for i in xrange(self.exons[e], self.exons[e+1]):
                            if i >= start+self.chromosomes[chrom] and i < end+self.chromosomes[chrom]:
                                coverage[i-self.chromosomes[chrom]-start] += junc.coverage[indexInCov + i - self.exons[e]]
                    indexInCov += self.exons[e+1] - self.exons[e]

        return coverage

    def getCoverageFromRLE(RLE, start, end):
        ''' Return the coverage vector expanded from the RLE vector between start and end.
        '''
        coverage = []

        index = 0
        offset = RLE[index][1]
        while offset < start:
            index += 1
            offset += RLE[index][1]

        coverage += [RLE[index][0]] * (min(offset, end) - start)

        while offset < end:
            index += 1
            nextOffset = offset + RLE[index][1]

            if nextOffset <= end:
                coverage += [RLE[index][0]] * RLE[index][1]
            else:
                coverage += [RLE[index][0]] * (end - offset)

            offset = nextOffset

        print coverage


    def getAvgCoverage(chrom, start=None, end=None):
        ''' Return average coverage over the given interval
        '''
        if not chrom in self.chromosomes:
            print 'Error! Chromosome name not recognized!'
            print 'Chromosomes: ' + ', '.join(self.chromosomes.keys())
            exit()

        cov = self.getCoverage(chrom, start, end)
        covSum = 0
        for c in cov:
            covSum += c
        return float(covSum) / len(cov)

    def getAlignments(chrom, start=None, end=None):
        ''' Return a list of all alignments that lie entirely within the given interval
        '''
        if not chrom in self.chromosomes:
            print 'Error! Chromosome name not recognized!'
            print 'Chromosomes: ' + ', '.join(self.chromosomes.keys())
            exit()

        startExon = 0
        endExon = len(self.exons)
        while start > self.exons[startExon+1]:
            startExon += 1
        while end <= self.exons[endExon-1]:
            endExon -= 1

    def getOverlappingAlignments(chrom, start=None, end=None):
        ''' Return a list of all alignments that overlap a given interval
        '''
        if not chrom in self.chromosomes:
            print 'Error! Chromosome name not recognized!'
            print 'Chromosomes: ' + ', '.join(self.chromosomes.keys())
            exit()

        return []

    def getGeneBoundaries(chrom):
        ''' Return a list of all gene boundaries in a chromosome
        '''
        if not chrom in self.chromosomes:
            print 'Error! Chromosome name not recognized!'
            print 'Chromosomes: ' + ', '.join(self.chromosomes.keys())
            exit()

        return []

