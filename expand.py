#! /usr/bin/env python
import alignments
import read
import pairedread
import junction

class Expander:
    aligned = None

    def expand(self, compressedFilename, uncompressedFilename):
        ''' Expand both spliced and unspliced alignments
        '''

        with open(compressedFilename, 'r') as f:
            self.aligned = alignments.Alignments(self.readHeader(f))

            self.expandSpliced(f)

            if self.aligned.exons == None:
                self.aligned.exons = [0, self.aligned.length]

            self.expandUnspliced(f)

            with open(uncompressedFilename, 'w') as f:
                self.aligned.writeSAM(f)

    def expandSpliced(self, f):
        ''' Expand a file containing compressed spliced alignments
        '''
        self.aligned.exons = None
        junc = None

        self.aligned.spliced = []
        for line in f:
            if self.aligned.exons == None:
                # First line contains exons
                self.aligned.exons = [int(e) for e in line.rstrip().split('\t')]
            else:
                # Remaining lines grouped by junction

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

    def readHeader(self, filehandle):
        chromInfo = filehandle.readline().rstrip().split('\t')
        chromosomes = dict()
        for c in chromInfo:
            ch = c.split(',')
            chromosomes[ch[0]] = int(ch[1])
        return chromosomes