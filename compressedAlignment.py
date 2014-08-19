#! /usr/bin/env python

import bisect

class Junction:
    readLens = None
    lensLeft = None
    lensRight = None
    coverage = []

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
                        junctionExons = [int(e) for e in key[:-1]]
                        xs = key[-1]
                        length = 0
                        for e in junctionExons:
                            length += self.exons[e+1] - self.exons[e]

                        junc = junction.Junction(junctionExons, length)
                        junc.readLens = None
                        junc.lensLeft = None
                        junc.lensRight = None

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
            lineNum = 0
            segmentNum = 0
            readLens = None
            lensLeft = None
            lensRight = None
            coverage = []

            # Read read lengths and coverage vector
            for line in f:
                if line[0] == '>':
                    if not readLens == None and len(readLens) > 0:

                        unpaired, paired = self.findReads(readLens, lensLeft, lensRight, coverage)
                        start = self.exons[segmentNum-1]


                        for r in unpaired:
                            self.unspliced.append(read.Read(self.getChromosome(r[0]+start), [[r[0]+start, r[1]+start]]))

                        for p in paired:
                            self.paired.append(pairedread.PairedRead(self.getChromosome(p[0][0]+start), [[p[0][0]+start, p[0][1]+start]],   \
                                                                     self.getChromosome(p[1][0]+start), [[p[1][0]+start, p[1][1]+start]] ))

                    lineNum = 0
                    segmentNum += 1

                    coverage = []

                    # First line contains read lengths
                    readLens = dict()
                    for row in line[1:].rstrip().split('\t'):
                        if len(row) > 0:
                            readLen = row.split(',')
                            readLens[int(readLen[0])] = int(readLen[1])
                elif lineNum == 1:
                    # Second line after '>' contains left read length distribution
                    lensLeft = dict()
                    if len(line.rstrip()) > 0:
                        for length in line.rstrip().split('\t'):
                            length = length.split(',')
                            lensLeft[int(length[0])] = int(length[1])
                elif lineNum == 2:
                    # Third line after '>' contains right read length distribution
                    lensRight = dict()
                    if len(line.rstrip()) > 0:
                        for length in line.rstrip().split('\t'):
                            length = length.split(',')
                            lensRight[int(length[0])] = int(length[1])
                else:
                    # Rest of lines contain run-length-encoded coverage vector
                    row = line.rstrip().split("\t")

                    if len(row) == 1:
                        length = 1
                    else:
                        length = int(row[1])

                    coverage += [int(row[0])] * length

                lineNum += 1

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
            exit()

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

