#! /usr/bin/env python
import alignments
import read
import pairedread
import junction
import time

class Expander:
    aligned = None

    def expand(self, compressedFilename, uncompressedFilename, binary=False):
        ''' Expand both spliced and unspliced alignments
        '''

        with open(compressedFilename, 'r') as f:
            self.aligned = alignments.Alignments(self.readHeader(f))

            self.aligned.exons = self.readExons(f)
            self.skipIndex(f)
            self.expandSpliced(f)

            self.expandUnspliced(f)

            with open(uncompressedFilename, 'w') as f2:
                self.aligned.writeSAM(f2)

    def expandSpliced(self, f):
        ''' Expand a file containing compressed spliced alignments
        '''

        #self.aligned.exons = None
        junc = None

        self.aligned.spliced = []
        for line in f:
            # Lines grouped by junction
            if line[0] == '>' or line[0] == '/':
                if not junc == None:
                    # process junction
                    unpaired, paired = self.aligned.findReads(junc.readLens, junc.lensLeft, junc.lensRight, junc.coverage)

                    juncBounds = []
                    for j in junctionExons:
                        juncBounds.append([self.aligned.exons[j], self.aligned.exons[j+1]])

                    # Offset of start of junction from beginning of chromosome
                    juncOffset = juncBounds[0][0]

                    # marks indices in junction coverage vector where exons are
                    mapping = [0]
                    # offsets of start of each exon in junction
                    offsets = [0]
                    for j in xrange(1,len(juncBounds)):
                        mapping.append(mapping[-1] + juncBounds[j-1][1] - juncBounds[j-1][0])
                        offsets.append(juncBounds[j][0] - juncBounds[0][0])

                    for r in unpaired:
                        start = r[0]
                        i = 0
                        while i < (len(mapping)-1) and mapping[i+1] <= start:
                            i += 1
                        start += offsets[i] - mapping[i]

                        end = r[1]
                        j = 0
                        while j < (len(mapping)-1) and mapping[j+1] < end:
                            j += 1
                        end += offsets[j] - mapping[j]
                        
                        readExons = []
                        if i == j:
                            readExons.append( [start+juncOffset, end+juncOffset] )
                        else:
                            readExons.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                            for x in xrange(i+1,j):
                                readExons.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                            readExons.append( [offsets[j]+juncOffset, end+juncOffset] ) 

                        self.aligned.spliced.append(read.Read(self.aligned.getChromosome(readExons[0][0]), readExons, junc.xs, junc.NH))

                    for p in paired:
                        start = p[0][0]
                        i = 0
                        while i < (len(mapping)-1) and mapping[i+1] <= start:
                            i += 1
                        start += offsets[i] - mapping[i]

                        end = p[0][1]
                        j = 0
                        while j < (len(mapping)-1) and mapping[j+1] < end:
                            j += 1
                        end += offsets[j] - mapping[j]
                        
                        readExonsA = []
                        if i == j:
                            readExonsA.append( [start+juncOffset, end+juncOffset] )
                        else:
                            readExonsA.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                            for x in xrange(i+1,j):
                                readExonsA.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                            readExonsA.append( [offsets[j]+juncOffset, end+juncOffset] )

                        start = p[1][0]
                        i = 0
                        while i < (len(mapping)-1) and mapping[i+1] <= start:
                            i += 1
                        start += offsets[i] - mapping[i]

                        end = p[1][1]
                        j = 0
                        while j < (len(mapping)-1) and mapping[j+1] < end:
                            j += 1
                        end += offsets[j] - mapping[j]
                        
                        readExonsB = []
                        if i == j:
                            readExonsB.append( [start+juncOffset, end+juncOffset] )
                        else:
                            readExonsB.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                            for x in xrange(i+1,j):
                                readExonsB.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                            readExonsB.append( [offsets[j]+juncOffset, end+juncOffset] )

                        self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(readExonsA[0][0]), readExonsA,  \
                                                                 self.aligned.getChromosome(readExonsB[0][0]), readExonsB, junc.xs, junc.NH))

                if line[0] == '/':
                    return

                # Start of a new junction
                key = line[1:].rstrip().split('\t')
                junctionExons = [int(e) for e in key[:-2]]
                length = 0

                for e in junctionExons:
                    length += self.aligned.exons[e+1] - self.aligned.exons[e]

                junc = junction.Junction(junctionExons, length)
                junc.readLens = None
                junc.lensLeft = None
                junc.lensRight = None

                junc.xs = key[-2]
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
                    #junc.coverage += [[val, length]]

    def expandUnspliced(self, f):
        ''' Expand a file containing compressed unspliced alignments
        '''

        NH = 1
        readLens = dict()
        exonStart = 0
        exonEnd = 0
        for line in f:
            if line[0] == '#':
                if len(readLens) > 0:
                    unpaired, paired = self.aligned.findReads(readLens, lensLeft, lensRight, coverage[exonStart:exonEnd])

                    for r in unpaired:
                        self.aligned.unspliced.append(read.Read(self.aligned.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))
                    for p in paired:
                        self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                                 self.aligned.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH))

                NH = int(line.rstrip()[1:])

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
                #coverage += [[int(row[0]), length]]

            elif line[0] == '>':
                RLE_segment = False
                if len(readLens) > 0:
                    # Process previous exon
                    unpaired, paired = self.aligned.findReads(readLens, lensLeft, lensRight, coverage[exonStart:exonEnd])

                    for r in unpaired:
                        self.aligned.unspliced.append(read.Read(self.aligned.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))

                    for p in paired:
                        self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                                 self.aligned.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH ))
                #exonId += 1
                lineNum = 0
                exonId = int(line[1:].rstrip())
                bounds = line[1:].rstrip().split('-')
                exonStart = self.aligned.exons[exonId]
                exonEnd = self.aligned.exons[exonId+1]

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
            unpaired, paired = self.aligned.findReads(readLens, lensLeft, lensRight, coverage[exonStart:exonEnd])

            for r in unpaired:
                self.aligned.unspliced.append(read.Read(self.aligned.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))

            for p in paired:
                self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                         self.aligned.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH))


    def getCoverage(self, compressedFile, chrom, start=None, end=None):
        with open(compressedFile, 'r') as f:
            chromosomes = self.readHeader(f)
            if not chrom in chromosomes:
                print 'Error! Chromosome name not recognized!'
                print 'Chromosomes: ' + ', '.join(chromosomes.keys())
                exit()

            self.aligned = alignments.Alignments(chromosomes)
            self.aligned.exons = self.readExons(f)

            splicedIndex, unsplicedIndex = self.readIndex(f)
            startPos = f.tell()

            if start == None or end == None:
                start = 0
                end = chromosomes[chrom]

            for k in sorted(chromosomes.keys()):
                if not k == chrom:
                    start += chromosomes[k]
                    end += chromosomes[k]
                else:
                    break

            coverage = [0.0] * (end-start)
            coverage = self.getSplicedCoverage(f, coverage, start, end, splicedIndex, startPos)
            coverage = self.getUnsplicedCoverage(f, coverage, start, end, unsplicedIndex, startPos)

        return coverage

    def getSplicedCoverage(self, filehandle, coverage, start, end, splicedIndex=None, startPos=0):
        if not splicedIndex == None:
            numExons = len(self.aligned.exons)
            i = 0
            while i < numExons and self.aligned.exons[i+1] <= start:
                i += 1
            j = i
            while j < numExons and self.aligned.exons[j] < end:
                j += 1

            minPos = splicedIndex[i]
            for n in xrange(i,j):
                if splicedIndex[n] < minPos:
                    minPos = splicedIndex[n]
            filehandle.seek(startPos + minPos)

        junc = None
        relevant = False

        relevantExonsStart = 0
        while self.aligned.exons[relevantExonsStart+1] < start:
            relevantExonsStart += 1
        relevantExonsEnd = relevantExonsStart
        while self.aligned.exons[relevantExonsEnd] < end:
            relevantExonsEnd += 1

        for line in filehandle:
            # Lines grouped by junction
            if line[0] == '/':
                return coverage

            elif line[0] == '>':
                # Start of a new junction
                key = line[1:].rstrip().split('\t')
                junctionExons = [int(e) for e in key[:-2]]

                relevant = False
                passedRange = True
                for e in junctionExons:
                    if e < relevantExonsEnd:
                        passedRange = False
                        if e >= relevantExonsStart:
                            relevant = True
                if passedRange:
                    break

                NH = 1.0 / float(key[-1])

                count = 0
                currExon = 0
                currExonLen = self.aligned.exons[junctionExons[0]+1] - self.aligned.exons[junctionExons[0]]
                offsetInExon = 0

            elif relevant:
                if count < 3:
                    count += 1
                else:
                    # Rest of lines contain run-length-encoded coverage vector
                    row = line.rstrip().split('\t')
                    if len(row) == 1:
                        length = 1
                    else:
                        length = int(row[1])

                    prevOffset = offsetInExon+self.aligned.exons[junctionExons[currExon]]
                    offsetInExon += length
                    while offsetInExon > currExonLen:
                        if junctionExons[currExon] >= relevantExonsStart and junctionExons[currExon] < relevantExonsEnd:
                            val = float(row[0]) * NH
                            if val > 0:
                                for i in xrange(max(prevOffset, self.aligned.exons[junctionExons[currExon]]), self.aligned.exons[junctionExons[currExon]+1]):
                                    if i >= start and i < end:
                                        coverage[i-start] += val

                        currExon += 1
                        offsetInExon -= currExonLen
                        currExonLen = self.aligned.exons[junctionExons[currExon]+1] - self.aligned.exons[junctionExons[currExon]]

                    if junctionExons[currExon] >= relevantExonsStart and junctionExons[currExon] < relevantExonsEnd:
                        val = float(row[0]) * NH
                        if val > 0:
                            for i in xrange(max(prevOffset, self.aligned.exons[junctionExons[currExon]]), self.aligned.exons[junctionExons[currExon]]+offsetInExon):
                                if i >= start and i < end:
                                    coverage[i-start] += val
        return coverage
    
    def getUnsplicedCoverage(self, filehandle, coverage, start, end, unsplicedIndex=None, startPos=0):
        segmentLength = len(coverage)
        
        sectionLen = 100000
        breakpointId = int(start / sectionLen)

        #for exonId in xrange(len(self.aligned.exons)):
        #    if self.aligned.exons[exonId] > start:
        #        break
        #exonId -= 1

        for k,v in unsplicedIndex.items():
            filehandle.seek(startPos + v[breakpointId])

            NH = 1.0 / float(k)

            offset = breakpointId * sectionLen
            line = filehandle.readline().rstrip()
            while offset < end and not line[0] == '>':
                row = line.split('\t')
                if len(row) == 1:
                    length = 1
                else:
                    length = int(row[1])

                newOffset = offset+length
                if newOffset > start:
                    val = float(row[0]) * NH
                    if val > 0:
                        for i in xrange(max(offset-start, 0), min(newOffset-start, segmentLength)):
                            coverage[i] += val
                offset = newOffset
                line = filehandle.readline().rstrip()
        
        return coverage

    def readExons(self, filehandle):
        line = filehandle.readline().rstrip()
        return [int(e) for e in line.split('\t')]

    def readHeader(self, filehandle):
        chromInfo = filehandle.readline().rstrip().split('\t')
        chromosomes = dict()
        for c in chromInfo:
            ch = c.split(',')
            chromosomes[ch[0]] = int(ch[1])
        return chromosomes

    def readIndex(self, filehandle):
        '''
            All index lines begin with #
        '''

        line = filehandle.readline().rstrip()
        splicedIndex = [int(i) for i in line.split('\t')[1:]]

        unsplicedIndex = dict()
        while True:
            offset = filehandle.tell()
            line = filehandle.readline().rstrip()
            if not line[0] == '#':
                # reset pointer
                filehandle.seek(offset)
                return splicedIndex, unsplicedIndex
            else:
                row = line[1:].split('\t')

                unsplicedIndex[int(row[0])] = [int(i) for i in row[1:]]

        return splicedIndex, unsplicedIndex

    def skipIndex(self, filehandle):
        '''
            All index lines begin with #
        '''
        while True:
            offset = filehandle.tell()
            if not filehandle.readline()[0] == '#':
                # reset pointer
                filehandle.seek(offset)

                return