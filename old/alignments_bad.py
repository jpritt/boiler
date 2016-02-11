#! /usr/bin/env python

import bisect
import bucket
import read
import pairedread
import os.path
import sys

import time

class Alignments:
    ''' A set of reads aligned to a genome '''

    def __init__(self, chromosomes):
        ''' Initialize a genome for alignments

            chromosomes: A dictionary with keys corresponding to chromosome names and values corresponding to lengths
        '''
        self.chromosomes = chromosomes

        # Initialize exon breaks between all chromosomes
        self.exons = [0]

        # Offset of each chromosome from the start of the genome
        self.chromOffsets = dict()

        nextOffset = 0
        for c in chromosomes.keys():
            self.chromOffsets[c] = nextOffset
            nextOffset += chromosomes[c]
            self.exons += [nextOffset]

        self.exons = set(self.exons)

        self.unspliced = []
        self.spliced = []

        # paired reads for which the mate still needs to be found
        #self.unmatched = []
        self.unmatched = dict()

        self.paired = []

    def addUnspliced(self, read):
        self.unspliced += [read]

    def addSpliced(self, read):
        self.spliced += [read]

        # update list of exons
        alignment = read.exons
        for i in xrange(len(alignment)-1):
            self.exons.add(alignment[i][1])
            self.exons.add(alignment[i+1][0])

    def compress(self, file_prefix):
        ''' Compresses the alignments to 2 files, one for unspliced and one for spliced

            file_prefix: Prefix for all output file names
        '''

        #start = time.time()
        if len(self.spliced) > 0:
            self.compressSpliced(file_prefix + '.spliced.txt')
        #time1 = time.time() - start
        #print '\nSpliced (%d):\t%f s' % (len(self.spliced), time1)

        #start = time.time()
        if len(self.unspliced) > 0:
            self.compressUnspliced(file_prefix + '.unspliced.txt')
        #time1 = time.time() - start
        #print 'Unspliced (%d):\t%f s' % (len(self.unspliced), time1)

    def compressSpliced(self, filename):
        ''' Compress the spliced alignments to a single file

            filename: Name of file to compress to
        '''

        # Compute coverage levels across every exon junction
        junctions = dict()

        for r in self.spliced:
            exonIds = r.exonIds

            if not r.xs == None:
                # XS is defined for this read
                key = '\t'.join([str(e) for e in exonIds]) + '\t' + r.xs + '\t' + str(r.NH)
                if not key in junctions:
                    covLength = 0
                    for e in exonIds:
                        covLength += self.exons[e+1] - self.exons[e]
                    junctions[key] = bucket.Junction(exonIds, covLength)
                j = junctions[key]
                j.pairedLens = dict()
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
                        covLength += self.exons[e+1] - self.exons[e]
                    junctions[key1] = bucket.Junction(exonIds, covLength)
                    j = junctions[key1]
                    j.pairedLens = dict()

                    key = key1

            # update junction coverage vector in dictionary
            if r.lenLeft == 0 and r.lenRight == 0:
                for i in xrange(r.startOffset, len(j.coverage)-r.endOffset):
                    j.coverage[i] += 1
            else:
                for i in xrange(r.startOffset, r.startOffset+r.lenLeft):
                    j.coverage[i] += 1
                for i in xrange(len(j.coverage)-r.endOffset-r.lenRight, len(j.coverage)-r.endOffset):
                    j.coverage[i] += 1

            # update readLens
            totalLen = r.readLen
            if r.lenLeft == 0 and r.lenRight == 0:
                if totalLen in j.readLens:
                    j.readLens[totalLen] += 1
                else:
                    j.readLens[totalLen] = 1
            else:
                if r.lenLeft in j.readLens:
                    j.readLens[r.lenLeft] += 1
                else:
                    j.readLens[r.lenLeft] = 1

                if r.lenRight in j.readLens:
                    j.readLens[r.lenRight] += 1
                else:
                    j.readLens[r.lenRight] = 1

                if totalLen in j.pairedLens:
                    j.pairedLens[totalLen] += 1
                else:
                    j.pairedLens[totalLen] = 1

        with open(filename, 'w') as f:

            # Write exons
            f.write('\t'.join([str(e) for e in self.exons]) + '\n')


            # Write junction information
            for key, junc in junctions.items():
                # write junction information
                f.write('>' + key + '\n')

                # write read lengths
                f.write('\t'.join( [str(k)+','+str(v) for k,v in junc.readLens.items()] ) + '\n')

                # Write left and right read lengths
                #f.write('\t'.join([str(k)+','+str(v) for k,v in junc.lensLeft.items()]) + '\n')
                #f.write('\t'.join([str(k)+','+str(v) for k,v in junc.lensRight.items()]) + '\n')
                f.write('\t'.join([str(k)+','+str(v) for k,v in junc.pairedLens.items()]) + '\n')

                # write coverage
                self.RLE(junc.coverage, f)

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
            length -= RLE[i][1]
            i += 1

        if length > 0:
            RLE = RLE[:i] + [ [RLE[i][0]+value, length], [RLE[i][0], RLE[i][1]-length] ] + RLE[i+1:]

        return RLE

    def compressUnspliced(self, filename):
        ''' Compress the unspliced alignments as a run-length-encoded coverage vector

            filename: Name of file to compress to
        '''

        
        # sort reads into exons
        readExons = dict()
        coverages = dict()
        for i in xrange(len(self.unspliced)):
            #r = self.unspliced[i]
            r = self.unspliced[len(self.unspliced) - i - 1]

            if not r.NH in readExons:
                readExons[r.NH] = []
                for n in xrange(len(self.exons)-1):
                    readExons[r.NH] += [[]]

                if r.NH == 1:
                    coverages[r.NH] = [0] * self.exons[-1]
                else:
                    coverages[r.NH] = [ [0, self.exons[-1]] ]

            j = bisect.bisect_right(self.exons, r.exons[0][0])-1
            readExons[r.NH][j] += [len(self.unspliced) - i - 1]

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
                if r.lenLeft == 0 and r.lenRight == 0:
                    coverages[r.NH] = self.updateRLE(coverages[r.NH], r.exons[0][0], r.exons[0][1]-r.exons[0][0], 1)
                else:
                    coverages[r.NH] = self.updateRLE(coverages[r.NH], r.exons[0][0], r.lenLeft, 1)
                    coverages[r.NH] = self.updateRLE(coverages[r.NH], r.exons[0][1]-r.lenRight, r.lenRight, 1)

        with open(filename, 'w') as f:
            for NH in readExons.keys():
                f.write('#' + str(NH) + '\n')

                cov = coverages[NH]
                reads = readExons[NH]

                # Write coverage vector
                if NH == 1:
                    self.RLE(cov, f)
                else:
                    for c in cov:
                        if c[1] == 1:
                            f.write(str(c[0]) + '\n')
                        else:
                            f.write(str(c[0]) + '\t' + str(c[1]) + '\n')

                # Write lengths for each exon
                for i in xrange(len(reads)):
                    if len(reads[i]) > 0:
                        length = self.exons[i+1] - self.exons[i]

                        exonStart = self.exons[i]

                        # Distribution of all read lengths
                        readLens = dict()

                        # Distribution of all gap lengths in paired-end reads
                        #lensLeft = dict()
                        #lensRight = dict()
                        pairedLens = dict()

                        for readId in reads[i]:
                            read = self.unspliced[readId]

                            alignment = read.exons

                            start = alignment[0][0] - exonStart
                            end = alignment[0][1] - exonStart

                            # update read lengths distribution
                            totalLen = end - start
                            if read.lenLeft == 0 and read.lenRight == 0:
                                if totalLen in readLens:
                                    readLens[totalLen] += 1
                                else:
                                    readLens[totalLen] = 1
                            else:
                                if read.lenLeft in readLens:
                                    readLens[read.lenLeft] += 1
                                else:
                                    readLens[read.lenLeft] = 1

                                if read.lenRight in readLens:
                                    readLens[read.lenRight] += 1
                                else:
                                    readLens[read.lenRight] = 1

                                if totalLen in pairedLens:
                                    pairedLens[totalLen] += 1
                                else:
                                    pairedLens[totalLen] = 1

                            '''
                            # update left and right read lengths
                            if read.lenLeft > 0 or read.lenRight > 0:
                                if read.lenLeft in lensLeft:
                                    lensLeft[read.lenLeft] += 1
                                else:
                                    lensLeft[read.lenLeft] = 1

                                if read.lenRight in lensRight:
                                    lensRight[read.lenRight] += 1
                                else:
                                    lensRight[read.lenRight] = 1
                            '''

                        # Write bounds
                        #f.write('>' + str(self.exons[i]) + '-' + str(self.exons[i+1]) + '\n')
                        f.write('>' + str(i) + '\n')

                        # Write read lengths
                        f.write('\t'.join([str(k)+','+str(v) for k,v in readLens.items()]) + '\n')

                        # Write left and right read lengths
                        #f.write('\t'.join([str(k)+','+str(v) for k,v in lensLeft.items()]) + '\n')
                        #f.write('\t'.join([str(k)+','+str(v) for k,v in lensRight.items()]) + '\n')
                        f.write('\t'.join([str(k)+','+str(v) for k,v in pairedLens.items()]) + '\n')

    def expand(self, file_prefix):
        ''' Expand both spliced and unspliced alignments
        '''
        #print 'Expanding spliced ' + file_prefix
        splicedName = file_prefix + '.spliced.txt'
        if os.path.isfile(splicedName):
            self.expandSpliced(splicedName)
        else:
            self.exons = [0, self.length]

        #print 'Expanding unspliced'
        unsplicedName = file_prefix + '.unspliced.txt'
        if os.path.isfile(unsplicedName):
            self.expandUnspliced(unsplicedName)

    def expandSpliced(self, filename):
        ''' Expand a file containing compressed spliced alignments
        '''
        print 'Expanding spliced'
        with open(filename, "r") as f:
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
                            # process junction
                            unpaired, paired = self.findReads(junc.readLens, junc.pairedLens, junc.coverage)

                            juncBounds = []
                            for j in junctionExons:
                                juncBounds.append([self.exons[j], self.exons[j+1]])

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

                                self.spliced.append(read.Read(self.getChromosome(readExons[0][0]), readExons, junc.xs, junc.NH))

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

                                self.paired.append(pairedread.PairedRead(self.getChromosome(readExonsA[0][0]), readExonsA,  \
                                                                         self.getChromosome(readExonsB[0][0]), readExonsB, junc.xs, junc.NH))

                        # Start of a new junction
                        key = line[1:].rstrip().split('\t')
                        junctionExons = [int(e) for e in key[:-2]]
                        length = 0
                        for e in junctionExons:
                            length += self.exons[e+1] - self.exons[e]

                        junc = bucket.Junction(junctionExons, length)
                        junc.readLens = None
                        junc.pairedLens = None

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
                        elif junc.pairedLens == None:
                            # Second line after '>' contains left read length distribution
                            junc.pairedLens = dict()
                            if len(line.rstrip()) > 0:
                                for length in line.rstrip().split('\t'):
                                    length = length.split(',')
                                    junc.pairedLens[int(length[0])] = int(length[1])
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
            unpaired, paired = self.findReads(junc.readLens, junc.pairedLens, junc.coverage)

            juncBounds = []
            for j in junctionExons:
                juncBounds.append([self.exons[j], self.exons[j+1]])

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
                    readExons.append( [start + juncOffset, end + juncOffset] )
                else:
                    readExons.append( [start + juncOffset, offsets[i]+mapping[i+1]-mapping[i] + juncOffset] )

                    for x in xrange(i+1,j):
                        readExons.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                    readExons.append( [offsets[j]+juncOffset, end+juncOffset] )

                self.spliced.append(read.Read(self.getChromosome(readExons[0][0]), readExons, junc.xs, junc.NH))

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

                self.paired.append(pairedread.PairedRead(self.getChromosome(readExonsA[0][0]), readExonsA,  \
                                                         self.getChromosome(readExonsB[0][0]), readExonsB, junc.xs, junc.NH ))

    def expandUnspliced(self, filename):
        ''' Expand a file containing compressed unspliced alignments
        '''
        print 'Expanding unspliced'

        with open(filename, "r") as f:
            NH = 1
            readLens = dict()
            exonStart = 0
            exonEnd = 0
            for line in f:
                if line[0] == '#':
                    print line[:-1]
                    if len(readLens) > 0:
                        unpaired, paired = self.findReads(readLens, pairedLens, coverage[exonStart:exonEnd])

                        for r in unpaired:
                            self.unspliced.append(read.Read(self.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))
                        for p in paired:
                            self.paired.append(pairedread.PairedRead(self.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                                     self.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH))

                    NH = int(line.rstrip()[1:])

                    lineNum = 0
                    readLens = dict()
                    pairedLens = dict()
                    coverage = []
                    RLE_segment = True
                    segmentId = 0

                elif RLE_segment and not line[0] == '>':
                    # Rest of lines contain run-length-encoded coverage vector
                    row = line.rstrip().split("\t")

                    if len(row) == 1:
                        length = 1
                    else:
                        length = int(row[1])

                    coverage += [int(row[0])] * length

                elif line[0] == '>':
                    segmentId += 1
                    print '  Segment ' + str(segmentId)

                    RLE_segment = False
                    if len(readLens) > 0:
                        # Process previous exon
                        unpaired, paired = self.findReads(readLens, pairedLens, coverage[exonStart:exonEnd], True)

                        for r in unpaired:
                            self.unspliced.append(read.Read(self.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))

                        for p in paired:
                            self.paired.append(pairedread.PairedRead(self.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                                     self.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH ))
                    #exonId += 1
                    lineNum = 0

                    exonId = int(line[1:].rstrip())
                    exonStart = self.exons[exonId]
                    exonEnd = self.exons[exonId+1]

                elif lineNum == 1:
                    # First line after '>' contains read lengths
                    readLens = dict()
                    for row in line.rstrip().split('\t'):
                        if len(row) > 0:
                            readLen = row.split(',')
                            readLens[int(readLen[0])] = int(readLen[1])
                elif lineNum == 2:
                    # Second line after '>' contains left read length distribution
                    pairedLens = dict()
                    if len(line.rstrip()) > 0:
                        for length in line.rstrip().split('\t'):
                            length = length.split(',')
                            pairedLens[int(length[0])] = int(length[1])
                lineNum += 1

            # Process the final exon
            if len(readLens) > 0:
                unpaired, paired = self.findReads(readLens, pairedLens, coverage[exonStart:exonEnd])

                for r in unpaired:
                    self.unspliced.append(read.Read(self.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))

                for p in paired:
                    self.paired.append(pairedread.PairedRead(self.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                             self.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH))
        print 'Finished'

    def finalizeExons(self):
        ''' Convert the set of exon boundaries to a list
        '''
        self.exons = list(sorted(self.exons))

    def finalizeReads(self):
        ''' Now that all exon boundaries are known, fix unspliced regions that cross exon boundaries
            and finalized paired-end reads
        '''

        if len(self.unmatched) > 0:
            # unmatched paired-end reads - shouldn't be any of these, but if there are, just make them unpaired reads
            print str(len(self.unmatched)) + ' unmatched'
            for k,v in self.unmatched.items():
                for r in k:
                    # unpaired read
                    if len(r.exons) == 1:
                        self.addUnspliced(r)
                    else:
                        self.addSpliced(r)

        # Splice any regions of unspliced reads that cross exon boundaries
        x = 0
        while x < len(self.unspliced):
            r = self.unspliced[x]

            exons = r.exons[0]
            spliceSites = self.findSpliceSites(exons[0], exons[1])

            if len(spliceSites) > 0:
                newExons = []
                newExons.append([exons[0], spliceSites[0]])
                for i in xrange(1, len(spliceSites)):
                    newExons.append([spliceSites[i-1], spliceSites[i]])
                newExons.append([spliceSites[-1], exons[1]])

                r.exons = newExons

                del self.unspliced[x]

                self.spliced += [r]
            else:
                x += 1


        # Compute list of junctions crossed by each spliced read
        for r in self.spliced:
            exons = r.exons

            # compute the length of the read
            r.readLen = 0

            # Find exons included in this read
            r.exonIds = []

            for segment in exons:
                r.readLen += segment[1] - segment[0]

                exonId = bisect.bisect_right(self.exons, segment[0])-1
                while self.exons[exonId] < segment[1]:
                    r.exonIds += [exonId]
                    exonId += 1

            # offset of start into first exon
            r.startOffset = exons[0][0] - self.exons[r.exonIds[0]]

            # offset of end from last exon
            r.endOffset =  self.exons[r.exonIds[-1] + 1] - exons[-1][1]

        # Convert paired reads to single long gapped reads for compressing
        for pair in self.paired:
            # the length of the coding portions of the read
            readLenLeft = 0
            readLenRight = 0

            # Find exons included in the first read from this pair
            pair.exonIdsA = []
            exonsA = pair.exonsA

            newExons = []
            for segment in exonsA:
                readLenLeft += segment[1] - segment[0]

                exonId = bisect.bisect_right(self.exons, segment[0])-1
                exonIds = []
                while self.exons[exonId] < segment[1]:
                    exonIds += [exonId]
                    exonId += 1

                pair.exonIdsA += exonIds

                spliceSites = [self.exons[e] for e in exonIds[1:]]
                if len(spliceSites) > 0:
                    newExons.append([segment[0], spliceSites[0]])
                    for i in xrange(1, len(spliceSites)):
                        newExons.append([spliceSites[i-1], spliceSites[i]])
                    newExons.append([spliceSites[-1], segment[1]])
                else:
                    newExons.append(segment)
            exonsA = newExons

            # Find exons included in the second read from this pair
            pair.exonIdsB = []
            exonsB = pair.exonsB

            newExons = []
            for segment in exonsB:
                readLenRight += segment[1] - segment[0]

                exonId = bisect.bisect_right(self.exons, segment[0])-1
                exonIds = []
                while self.exons[exonId] < segment[1]:
                    exonIds += [exonId]
                    exonId += 1

                pair.exonIdsB += exonIds

                spliceSites = [self.exons[e] for e in exonIds[1:]]
                if len(spliceSites) > 0:
                    newExons.append([segment[0], spliceSites[0]])
                    for i in xrange(1, len(spliceSites)):
                        newExons.append([spliceSites[i-1], spliceSites[i]])
                    newExons.append([spliceSites[-1], segment[1]])
                else:
                    newExons.append(segment)
            exonsB = newExons

            n = 0
            while n < len(pair.exonIdsB) and pair.exonIdsB[n] <= pair.exonIdsA[-1]:
                n += 1

            ## TODO: Fix this problem with exons not matching
            if n > 0 and not pair.exonIdsA[-n:] == pair.exonIdsB[:n]:
                continue

            # TODO: Replace gap with end lengths
            if n > 0:
                newRead = read.Read(pair.chromA, exonsA[:-1] + [[exonsA[-1][0], exonsB[n-1][1]]] + exonsB[n:], pair.xs, pair.NH)
                gap = 0
                for i in xrange(n):
                    gap += exonsB[i][0] - exonsA[i-n][1]
            else:
                newRead = read.Read(pair.chromA, exonsA+exonsB, pair.xs, pair.NH)
                gap = exonsB[0][0] - self.exons[pair.exonIdsB[0]] + self.exons[pair.exonIdsA[-1]+1] - exonsA[-1][1]
            newRead.exonIds = pair.exonIdsA + pair.exonIdsB[n:]

            newRead.lenLeft = readLenLeft
            newRead.lenRight = readLenRight
            newRead.readLen = readLenLeft + readLenRight + gap

            # offset of start into first exon
            newRead.startOffset = exonsA[0][0] - self.exons[pair.exonIdsA[0]]

            # offset of end from last exon
            newRead.endOffset = self.exons[pair.exonIdsB[-1] + 1] - exonsB[-1][1]

            if len(newRead.exonIds) == 1:
                self.unspliced += [newRead]
            else:
                self.spliced += [newRead]


    # Returns the element from list1 and the element from list2 with the smallest distance between them
    def findClosestVals(self, list1, list2):
        if list1[0] == list2[0]:
            return (0,0)
        elif list1[0] < list2[0]:
            temp = list1
            list1 = list2
            list2 = temp
            swapped = True
        else:
            swapped = False

        minRange = 0
        minDist = abs(list1[0] - list2[0])
        minVals = (0,0)

        len1 = len(list1)
        len2 = len(list2)
        for i in xrange(len1):
            while (minRange+1) < len2 and list2[minRange+1] < list1[i]:
                minRange += 1
            dist1 = abs(list2[minRange] - list1[i])
            if dist1 < minDist:
                minDist = dist1
                minVals = (i, minRange)

            if (minRange+1) < len2:
                dist2 = abs(list2[minRange+1] - list1[i])
                if dist2 < minDist:
                    minDist = dist2
                    minVals = (i, minRange+1)

        if swapped:
            return (minVals[1], minVals[0])
        else:
            return minVals

    def findReads(self, readLens, pairedLens, coverage, debug=False):
        ''' Find the set of reads that most closely matches the distribution of readLens and the coverage vector
        '''

        #reads = self.findReadsInCoverage_v1(coverage, readLens, None)
        reads = self.findReadsInCoverage_v2(coverage, readLens, None)
        #reads = self.findReadsInCoverage_v3(coverage, readLens, boundaries)
        if debug:
            print '  Found reads, pairing'

        #print ''
        #print readLens
        #print str(len(reads)) + ' recovered'
        #print pairedLens

        if reads == None:
            print 'Error finding reads!'
            exit()

        # Pair left and right read lengths with long fragments
        # First sort reads and lengths by decreasing length
        reads.sort(key=lambda r: r[0], reverse=False)
        pairedLensSorted = sorted(pairedLens.keys(), reverse=True)

        pairedReads = []

        lenId = 0
        while lenId < len(pairedLensSorted) and len(reads) > 0:
            #print ''
            #print reads
            #print pairedLens
            #print pairedLensSorted

            i,j = self.findBestPair(reads, pairedLensSorted[0])
            #print '%d, %d' % (i,j)

            pairedReads.append( [reads[i], reads[j]] )

            if i < j:
                del reads[j]
                del reads[i]
            else:
                del reads[i]
                del reads[j]

            pairedLens[pairedLensSorted[lenId]] -= 1
            if pairedLens[pairedLensSorted[lenId]] == 0:
                lenId += 1
            elif pairedLens[pairedLensSorted[lenId]] < 0:
                print 'Error! Negative value!'

        return reads, pairedReads

    def findBestPair(self, reads, pairLen):
        closestDiff = None
        closesti = 0
        closestj = 0
        for i in xrange(len(reads)-1):
            for j in xrange(i+1, len(reads)):
                if reads[i][1] > reads[j][0]:
                    lenDiff = abs(pairLen - (reads[i][1]-reads[j][0]))
                    if lenDiff == 0:
                        return j,i
                    elif closestDiff == None or lenDiff < closestDiff:
                        closestDiff = lenDiff
                        closesti = j
                        closestj = i

                if reads[j][1] > reads[i][0]:
                    lenDiff = abs(pairLen - (reads[j][1] - reads[i][0]))
                    if lenDiff == 0:
                        return i,j
                    elif closestDiff == None or lenDiff < closestDiff:
                        closestDiff = lenDiff
                        closesti = i
                        closestj = j
        if closesti == closestj:
            print '%d, %d' % (closesti, closestj)
        return closesti, closestj

    def findReadsInCoverage_v1(self, coverage, readLens, boundaries=None):
        ''' Given a coverage vector, return a set of reads the closely fits the coverage vector as well the distribution of read lengths.
            The algorithm creates new reads greedily, starting from both ends at once to ensure the ends of the vector are well marked.
            This algorithm is guaranteed to return a set of reads that covers every base at least to the corresponding depth of the coverage vector.
            In many cases the algorithm will overcompensate by creating extra reads to make sure every base in the coverage vector is covered.
            In such cases new reads have length equal to the median read length in the input distribution.

            coverage: Coverage vector containing the accumulation of many reads
            readLens: Dictionary containing the distribution of all read lengths in the coverage vector
            boundaries: Indicies in the coverage vector which reads cannot cross (e.g. exon boundaries in the unspliced coverage vector)
        '''

        if boundaries == None:
            boundaries = [0, len(coverage)]
        else:
            boundaries = list(boundaries)
        if not len(boundaries) == 2:
            print '%d boundaries' % len(boundaries)
        boundBottom = 0
        boundTop = len(boundaries)-1

        lens = readLens.keys()

        # Find max and mode read lengths
        maxLen = max(lens)

        # Read lengths sorted by frequency, largest to smallest
        lensSorted = sorted(readLens, key=readLens.get, reverse=True)


        reads = []

        # start and end of coverage window
        # Keep finding reads from both ends until they meet in the middle
        start = 0
        while coverage[start] <= 0:
            start += 1
        end = len(coverage)
        while coverage[end-1] <= 0:
            end -= 1

        while end > start:
            while boundaries[boundBottom] <= start:
                boundBottom += 1

            # find a read from the beginning
            readStart = start
            readEnd = start

            currMaxLen = maxLen

            closestEndpoint = None
            for length in xrange(1, min(maxLen+1, boundaries[boundBottom]-start)):
                if (readStart+length == end) or (readStart+length < end and coverage[readStart + length] < coverage[readStart + length - 1]):
                    if length in readLens:
                        readEnd = readStart + length
                        reads.append([readStart, readEnd])

                        readLens[length] -= 1

                        # reorder sorted lengths
                        for i in xrange(len(lensSorted)):
                            if lensSorted[i] == length:
                                break
                        j = i+1
                        while j < len(readLens) and readLens[lensSorted[j]] > readLens[lensSorted[i]]:
                            j += 1
                        if j > i+1:
                            lensSorted = lensSorted[:i] + lensSorted[i+1:j] + [lensSorted[i]] + lensSorted[j:]

                        if readLens[length] == 0:
                            del readLens[length]

                        break
                    else:
                        if closestEndpoint == None:
                            closestEndpoint = readStart + length

                # Don't extend into section where coverage is 0
                if coverage[readStart+length] == 0:
                    currMaxLen = length
                    break
            if readEnd == readStart:
                if closestEndpoint == None:
                    lenId = 0
                    length = lensSorted[lenId]
                    readEnd = readStart + length

                    while readEnd > len(coverage) and lenId < (len(lensSorted)-1):
                        lenId += 1
                        length = lensSorted[lenId]
                        readEnd = readStart + length
                    if readEnd > len(coverage):
                        # No read lengths fit within the end of the vector
                        readEnd = len(coverage)

                    reads.append([readStart, readEnd])

                    if length in readLens:
                        readLens[length] -= 1

                        # reorder sorted lengths
                        for i in xrange(len(lensSorted)):
                            if lensSorted[i] == length:
                                break
                        j = i+1
                        while j < len(readLens) and readLens[lensSorted[j]] > readLens[lensSorted[i]]:
                            j += 1
                        if j > i+1:
                            lensSorted = lensSorted[:i] + lensSorted[i+1:j] + [lensSorted[i]] + lensSorted[j:]

                        if readLens[length] == 0:
                            del readLens[length]
                else:
                    readEnd = closestEndpoint
                    reads.append([readStart, readEnd])

            # Update coverage vector
            for i in xrange(readStart, readEnd):
                coverage[i] -= 1

            # update start
            while start < end and coverage[start] <= 0:
                start += 1
            while end > start and coverage[end-1] <= 0:
                end -= 1


            if end > start:
                while boundaries[boundTop] >= end:
                    boundTop -= 1

                # find a read from the end
                readEnd = end
                readStart = end

                closestEndpoint = None
                for length in xrange(1, min(maxLen+1, end-boundaries[boundTop])):
                    if (end-length == start) or (end-length > start and coverage[end - length] > coverage[end - length - 1]):
                        if length in readLens:
                            readStart = readEnd - length
                            reads.append([readStart, readEnd])

                            readLens[length] -= 1
                                

                            # reorder sorted lengths
                            for i in xrange(len(lensSorted)):
                                if lensSorted[i] == length:
                                    break
                            j = i+1

                            while j < len(readLens) and readLens[lensSorted[j]] > readLens[lensSorted[i]]:
                                j += 1
                            if j > i+1:
                                lensSorted = lensSorted[:i] + lensSorted[i+1:j] + [lensSorted[i]] + lensSorted[j:]


                            if readLens[length] == 0:
                                del readLens[length]

                            break
                        else:
                            if closestEndpoint == None:
                                closestEndpoint = readEnd - length

                if readStart == readEnd:
                    if closestEndpoint == None:
                        length = lensSorted[0]
                        readStart = readEnd - length
                        reads.append([readStart, readEnd])

                        if length in readLens:
                            readLens[length] -= 1

                            # reorder sorted lengths
                            for i in xrange(len(lensSorted)):
                                if lensSorted[i] == length:
                                    break

                            j = i+1
                            while j < len(readLens) and readLens[lensSorted[j]] > readLens[lensSorted[i]]:
                                j += 1
                            if j > i+1:
                                lensSorted = lensSorted[:i] + lensSorted[i+1:j] + [lensSorted[i]] + lensSorted[j:]

                            if readLens[length] == 0:
                                del readLens[length]
                    else:
                        readStart = closestEndpoint
                        reads.append([readStart, readEnd])

                for i in xrange(readStart, readEnd):
                    coverage[i] -= 1

                # update end
                while coverage[end-1] <= 0 and end > start:
                    end -= 1
                while coverage[start] <= 0 and start < end:
                    start += 1
        return reads

    def findReadsInCoverage_v2(self, coverage, readLens, boundaries=None):
        ''' Given a coverage vector, return a set of reads the closely fits the coverage vector as well the distribution of read lengths.
            The algorithm creates new reads greedily, starting from both ends at once to ensure the ends of the vector are well marked.
            This algorithm is guaranteed to return a set of reads that matches the input read length distribution exactly. However these reads
              may not replicate the input coverage vector exactly.

            coverage: Coverage vector containing the accumulation of many reads
            readLens: Dictionary containing the distribution of all read lengths in the coverage vector
            boundaries: Indicies in the coverage vector which reads cannot cross (e.g. exon boundaries in the unspliced coverage vector)
        '''

        countReads = 0
        for length, freq in readLens.items():
            countReads += freq

        if boundaries == None:
            boundaries = [0, len(coverage)]
        else:
            boundaries = list(boundaries)
        boundBottom = 0
        boundTop = len(boundaries)-1


        # Read lengths sorted first by frequency, largest to smallest, then by length (smallest to largest)
        lensSorted = [v[0] for v in sorted(readLens.iteritems(), key = lambda(k,v): (v, -k), reverse=True)]


        reads = []

        # start and end of coverage window
        # Keep finding reads from both ends until they meet in the middle
        start = 0
        while start < len(coverage) and coverage[start] <= 0:
            start += 1
        end = len(coverage)
        while end >= start and coverage[end-1] <= 0:
            end -= 1

        while end > start and len(lensSorted) > 0:
            #print 'Adding from start'
            while boundaries[boundBottom] <= start:
                boundBottom += 1

            # find a read from the beginning
            readStart = start
            readEnd = start

            readFound = False

            for length in lensSorted:
                readEnd = readStart + length

                if readEnd <= boundaries[boundBottom] and readEnd <= end and (readEnd == len(coverage) or coverage[readEnd-1] > coverage[readEnd]):
                    reads.append([readStart, readEnd])

                    readLens[length] -= 1

                    # reorder sorted lengths
                    for i in xrange(len(lensSorted)):
                        if lensSorted[i] == length:
                            break

                    if readLens[length] == 0:
                        del lensSorted[i]
                    else:
                        j = i+1

                        while j < len(lensSorted) and (readLens[lensSorted[j]] > readLens[lensSorted[i]] or (readLens[lensSorted[j]] == readLens[lensSorted[i]] and lensSorted[j] < lensSorted[i])):
                            j += 1
                        if j > i+1:
                            lensSorted = lensSorted[:i] + lensSorted[i+1:j] + [lensSorted[i]] + lensSorted[j:]

                    readFound = True
                    break
            if not readFound:
                # No good end point found; add the most common read from from the readLen distribution
                
                i = 0
                while i+1 < len(lensSorted) and lensSorted[i] > len(coverage)-readStart and readLens[lensSorted[i]] > 0:
                    i += 1
                if lensSorted[i] > len(coverage)-readStart:
                    i = 0
                    readStart = len(coverage) - lensSorted[0]

                length = lensSorted[i]
                readEnd = readStart + length
                reads.append([readStart, readEnd])

                readLens[length] -= 1

                # reorder sorted lengths
                for i in xrange(len(lensSorted)):
                    if lensSorted[i] == length:
                        break

                if readLens[length] == 0:
                    del lensSorted[i]
                else:
                    j = i+1

                    while j < len(lensSorted) and (readLens[lensSorted[j]] > readLens[lensSorted[i]] or (readLens[lensSorted[j]] == readLens[lensSorted[i]] and lensSorted[j] < lensSorted[i])):
                        j += 1
                    if j > 1:
                        lensSorted = lensSorted[:i] + lensSorted[i+1:j] + [lensSorted[i]] + lensSorted[j:]

            # Update coverage vector
            for i in xrange(readStart, readEnd):
                coverage[i] -= 1

            # update start
            while start < end and coverage[start] <= 0:
                start += 1
            while end > start and coverage[end-1] <= 0:
                end -= 1


            if end > start and len(lensSorted) > 0:
                #print 'Adding from end'

                while boundaries[boundTop] >= end:
                    boundTop -= 1

                # find a read from the end
                readEnd = end
                readStart = end

                readFound = False

                for length in lensSorted:
                    readStart = readEnd - length
                    if readStart >= boundaries[boundTop] and (readStart == 0 or coverage[readStart] > coverage[readStart-1]):
                        readStart = readEnd - length
                        reads.append([readStart, readEnd])

                        readLens[length] -= 1
                            

                        # reorder sorted lengths
                        for i in xrange(len(lensSorted)):
                            if lensSorted[i] == length:
                                break

                        if readLens[length] == 0:
                            del lensSorted[i]
                        else:
                            j = i+1

                            while j < len(lensSorted) and (readLens[lensSorted[j]] > readLens[lensSorted[i]] or (readLens[lensSorted[j]] == readLens[lensSorted[i]] and lensSorted[j] < lensSorted[i])):
                                j += 1
                            if j > i+1:
                                lensSorted = lensSorted[:i] + lensSorted[i+1:j] + [lensSorted[i]] + lensSorted[j:]

                        readFound = True
                        break

                if not readFound:
                    # No good end point found; add the most common read from from the readLen distribution
                
                    i = 0
                    while i+1 < len(lensSorted) and lensSorted[i] > len(coverage)-readStart and readLens[lensSorted[i]] > 0:
                        i += 1
                    if lensSorted[i] > len(coverage)-readStart:
                        i = 0
                        readStart = len(coverage) - lensSorted[0]

                    length = lensSorted[i]
                    readEnd = readStart + length
                    reads.append([readStart, readEnd])

                    readLens[length] -= 1

                    # reorder sorted lengths
                    for i in xrange(len(lensSorted)):
                        if lensSorted[i] == length:
                            break

                    if readLens[length] == 0:
                        del lensSorted[i]
                    else:
                        j = i+1

                        while j < len(lensSorted) and (readLens[lensSorted[j]] > readLens[lensSorted[i]] or (readLens[lensSorted[j]] == readLens[lensSorted[i]] and lensSorted[j] < lensSorted[i])):
                            j += 1
                        if j > 1:
                            lensSorted = lensSorted[:i] + lensSorted[i+1:j] + [lensSorted[i]] + lensSorted[j:]

                for i in xrange(readStart, readEnd):
                    coverage[i] -= 1

                # update end
                while coverage[end-1] <= 0 and end > start:
                    end -= 1
                while coverage[start] <= 0 and start < end:
                    start += 1

        if not len(reads) == countReads:
            print 'Error! %d =/= %d!' % (countReads, len(reads))
            exit()
        #print ''
        return reads

    
    def findReadsInCoverage_v3(self, coverage, readLens, boundaries=None):
        ''' Given a coverage vector, return a set of reads that exactly fits the coverage vector as well the distribution of read lengths.
            Uses brute force to find the an set of reads that fits the read length distribution and matches the coverage vector

            coverage: Coverage vector containing the accumulation of many reads
            readLens: Dictionary containing the distribution of all read lengths in the coverage vector
            boundaries: Indicies in the coverage vector which reads cannot cross (e.g. exon boundaries in the unspliced coverage vector)
        '''

        start = 0
        while coverage[start] == 0 and start < len(coverage):
            start += 1
        end = len(coverage)
        while coverage[end-1] == 0:
            end -= 1
        if start == end:
            return []

        print 'Coverage:'
        self.RLE(coverage, sys.stdout)
        print 'Read Lengths: ' + str(readLens)
        reads = self.findReadsInCoverage_BruteForce2(coverage, start, end, readLens)
        print 'Reads: ' + str(reads)
        print ''
        return reads

    
    def findReadsInCoverage_BruteForce(self, coverage, start, end, readLens, level=0):

        length = end-start

        firstLength = readLens.keys()[0]
        # Only 1 read left: return that read iff coverage[start:end] is all 1s of the correct length
        if len(readLens) == 1 and readLens[firstLength] == 1:
            if length == firstLength:
                for i in xrange(start,end):
                    if not coverage[i] == 1:
                        return None
                return [[start,end]]

        # More than 1 read left: pick 1 recursively
        iteration = 0
        uniqueReads = readLens.keys()
        for r in uniqueReads:
            iteration += 1
            if level < 2:
                print str(level) + '  '*(level+1) + str(iteration) + ' / ' + str(len(uniqueReads))
            if iteration > len(readLens):
                print 'Error!'
                exit()

            # make sure read is short enough for available space
            if length >= r:
                readFits = True

                # Subtract read from coverage
                for i in xrange(start, start+r):
                    if coverage[i] == 0:
                        readFits = False
                        break
                    coverage[i] -= 1

                if not readFits:
                    #print '  '*level + 'Read doesn\'t fit'
                    for j in xrange(start, i):
                        coverage[j] += 1
                else:
                    # Update bounds of nonzero portion of coverage
                    #print '  '*level + 'New Coverage: ' + ','.join([str(c) for c in coverage])
                    newStart = start
                    while coverage[newStart] == 0 and newStart < end:
                        newStart += 1
                    newEnd = end
                    while coverage[newEnd-1] == 0 and newEnd > newStart:
                        newEnd -= 1
                    #print 'Start = %d, End = %d' % (newStart, newEnd)

                    readLens[r] -= 1
                    if readLens[r] == 0:
                        del readLens[r]

                    reads = self.findReadsInCoverage_BruteForce(coverage, newStart, newEnd, readLens, level+1)
                    # Add read back to coverage
                    for i in xrange(start, start+r):
                        coverage[i] += 1

                    if r in readLens:
                        readLens[r] += 1
                    else:
                        readLens[r] = 1
        
                    if not reads == None:
                        #print '  '*level + 'Works!'
                        return [[start, start+r]] + reads
        return None

    def findReadsInCoverage_v4(self, coverage, readLens, boundaries=None):
        ''' Given a coverage vector, return a set of reads that exactly fits the coverage vector as well the distribution of read lengths.
            Uses brute force to find the an set of reads that fits the read length distribution and matches the coverage vector

            coverage: Coverage vector containing the accumulation of many reads
            readLens: Dictionary containing the distribution of all read lengths in the coverage vector
            boundaries: Indicies in the coverage vector which reads cannot cross (e.g. exon boundaries in the unspliced coverage vector)
        '''

        start = 0
        end = len(coverage)
        while coverage[start] == 0 and start < end:
            start += 1
        while coverage[end-1] == 0 and end > start:
            end -= 1
        if start == end:
            return []

        print 'Coverage:'
        self.RLE(coverage, sys.stdout)
        print 'Read Lengths: ' + str(readLens)
        
        reads = []
        while start < end:
            newReads = self.findReadsInCoverage_BruteForce2(coverage, start, end, readLens)
            if newReads == None:
                print 'Error! No reads found!'
                exit()
            reads.append(newReads)

            while coverage[start] == 0 and start < end:
                start += 1
            while coverage[end-1] == 0 and end > start:
                end -= 1
        
        print 'Reads: ' + str(reads)
        print ''
        return reads

    def findReadsInCoverage_BruteForce2(self, coverage, start, end, readLens, level=0):

        
        totalLen = 0
        for r,num in readLens.items():
            totalLen += r*num
        covLen = sum(coverage)
        if not totalLen == sum(coverage):
            print 'Error! %d =/= %d' % (covLen, totalLen)
            exit()
        

        length = end-start

        firstLength = readLens.keys()[0]
        # Only 1 read left: return that read iff coverage[start:end] is all 1s of the correct length
        if len(readLens) == 1 and readLens[firstLength] == 1:
            if length == firstLength:
                for i in xrange(start,end):
                    if not coverage[i] == 1:
                        return None
                return [[start,end]]

        # More than 1 read left: pick 1 recursively
        iteration = 0

        # separate reads into likely and unlike for efficiency
        likely = []
        unlikely = []

        for r in readLens.keys():
            if coverage[r] > coverage[r+1]:
                likely += [r]
            else:
                unlikely += [r]

        for r in likely:
            iteration += 1
            if level < 2:
                print str(level) + '  '*(level+1) + str(iteration) + ' / ' + str(len(likely) + len(unlikely))
            if iteration > len(readLens):
                print 'Error!'
                exit()

            # make sure read is short enough for available space
            if length >= r:
                readFits = True

                # Subtract read from coverage
                for i in xrange(start, start+r):
                    if coverage[i] == 0:
                        readFits = False
                        break
                    coverage[i] -= 1

                if not readFits:
                    #print '  '*level + 'Read doesn\'t fit'
                    for j in xrange(start, i):
                        coverage[j] += 1
                else:
                    # Update bounds of nonzero portion of coverage
                    #print '  '*level + 'New Coverage: ' + ','.join([str(c) for c in coverage])
                    newStart = start
                    while coverage[newStart] == 0 and newStart < end:
                        newStart += 1
                    newEnd = end
                    while coverage[newEnd-1] == 0 and newEnd > newStart:
                        newEnd -= 1
                    #print 'Start = %d, End = %d' % (newStart, newEnd)

                    readLens[r] -= 1
                    if readLens[r] == 0:
                        del readLens[r]

                    reads = self.findReadsInCoverage_BruteForce2(coverage, newStart, newEnd, readLens, level+1)
                    # Add read back to coverage
                    for i in xrange(start, start+r):
                        coverage[i] += 1

                    if r in readLens:
                        readLens[r] += 1
                    else:
                        readLens[r] = 1
        
                    if not reads == None:
                        #print '  '*level + 'Works!'
                        return [[start, start+r]] + reads
        for r in unlikely:
            iteration += 1
            if level < 2:
                print str(level) + '  '*(level+1) + str(iteration) + ' / ' + str(len(likely) + len(unlikely))
            if iteration > len(readLens):
                print 'Error!'
                exit()

            # make sure read is short enough for available space
            if length >= r:
                readFits = True

                # Subtract read from coverage
                for i in xrange(start, start+r):
                    if coverage[i] == 0:
                        readFits = False
                        break
                    coverage[i] -= 1

                if not readFits:
                    #print '  '*level + 'Read doesn\'t fit'
                    for j in xrange(start, i):
                        coverage[j] += 1
                else:
                    # Update bounds of nonzero portion of coverage
                    #print '  '*level + 'New Coverage: ' + ','.join([str(c) for c in coverage])
                    newStart = start
                    while coverage[newStart] == 0 and newStart < end:
                        newStart += 1
                    newEnd = end
                    while coverage[newEnd-1] == 0 and newEnd > newStart:
                        newEnd -= 1
                    #print 'Start = %d, End = %d' % (newStart, newEnd)

                    readLens[r] -= 1
                    if readLens[r] == 0:
                        del readLens[r]

                    reads = self.findReadsInCoverage_BruteForce2(coverage, newStart, newEnd, readLens, level+1)
                    # Add read back to coverage
                    for i in xrange(start, start+r):
                        coverage[i] += 1

                    if r in readLens:
                        readLens[r] += 1
                    else:
                        readLens[r] = 1
        
                    if not reads == None:
                        #print '  '*level + 'Works!'
                        return [[start, start+r]] + reads
        return None

    def  findSpliceSites(self, start, end):
        ''' Look for any exon boundaries crossed by the given interval. 
            Return a list of all splice sites in the interval [start, end)
        '''

        startExon = bisect.bisect_right(self.exons, start)
        endExon = startExon
        while endExon < len(self.exons) and self.exons[endExon] <= end:
            endExon += 1

        return self.exons[startExon:endExon]

    def getChromosome(self, index):
        ''' Return chromosome name containing the given index from the whole-genome vector
        '''

        for c in self.chromosomes:
            if self.chromosomes[c] > index:
                return c
            else:
                index -= self.chromosomes[c]

    def getChromosomeAndIndex(self, index):
        ''' Return chromosome name containing the given index from the whole-genome vector
        '''

        for c in self.chromosomes:
            if self.chromosomes[c] > index:
                return c, index
            else:
                index -= self.chromosomes[c]

    def insertInOrder(sortedList, a):
        ''' Insert a in the correct place in a sorted list in increasing order
        '''
        i = 0
        while i < len(sortedList) and a > sortedList[i]:
            i += 1
        return sortedList[:i] + [a] + sortedList[i:]

    def processRead(self, read, pair_chrom=None, pair_index=None):
        ''' If read is unpaired, add it to the correct spliced or unspliced list of reads.
            If read is paired, find its pair or add it to a list to be found later. Once a pair of reads is found, add the combined read to the appropriate list of reads
        '''
        if read.exons == None:
            print 'Error! Read exons = None'
            exit()
        if len(read.exons) == 0:
            print 'Error! Exon length = 0'

        # Update read index based on chromosome
        offset = self.chromOffsets[read.chrom]
        for i in xrange(len(read.exons)):
            read.exons[i] = [read.exons[i][0]+offset, read.exons[i][1]+offset]

        # update list of exons
        alignment = read.exons
        if len(alignment) > 1:
            for i in xrange(len(alignment)-1):
                self.exons.add(alignment[i][1])
                self.exons.add(alignment[i+1][0])

        if pair_chrom == None:
            # unpaired read
            if len(read.exons) == 1:
                self.addUnspliced(read)
            else:
                self.addSpliced(read)

        else:
            pair_index += self.chromOffsets[pair_chrom]

            # TODO: Use binary search here for speed
            #self.unmatched = dict()

            if (pair_index, pair_chrom, read.exons[0][0], read.chrom, read.NH) in self.unmatched:
                match = self.unmatched[(pair_index, pair_chrom, read.exons[0][0], read.chrom, read.NH)][0]

                if len(self.unmatched[(pair_index, pair_chrom, read.exons[0][0], read.chrom, read.NH)]) == 1:
                    del self.unmatched[(pair_index, pair_chrom, read.exons[0][0], read.chrom, read.NH)]
                else:
                    del self.unmatched[(pair_index, pair_chrom, read.exons[0][0], read.chrom, read.NH)][0]

                xs = read.xs or match.xs
                NH = read.NH


                self.paired.append(pairedread.PairedRead(pair_chrom, match.exons, read.chrom, read.exons, xs, NH))
            else:
                if not (read.exons[0][0], read.chrom, pair_index, pair_chrom, read.NH) in self.unmatched:
                    self.unmatched[(read.exons[0][0], read.chrom, pair_index, pair_chrom, read.NH)] = [read]
                else:
                    self.unmatched[(read.exons[0][0], read.chrom, pair_index, pair_chrom, read.NH)] += [read]

    def RLE(self, vector, filehandle):
        val = vector[0]
        length = 0

        for v in vector:
            if v == val:
                length += 1
            else:
                if length == 1:
                    filehandle.write(str(val) + '\n')
                else:
                    filehandle.write(str(val) + '\t' + str(length) + '\n')
                val = v
                length = 1
        if length == 1:
            filehandle.write(str(val) + '\n')
        else:
            filehandle.write(str(val) + '\t' + str(length) + '\n')

    def writeSAM(self, filehandle):
        ''' Write all alignments to a SAM file
        '''
        readId = 0
        for read in self.unspliced:
            exons = read.exons
            chrom = read.chrom
            offset = self.chromOffsets[chrom]
            filehandle.write(read.chrom+':'+str(readId) + '\t0\t' + read.chrom + '\t' + str(exons[0][0]-offset) + '\t50\t' + str(exons[0][1]-exons[0][0]) + 'M\t*\t0\t0\t*\t*\tNH:i:' + str(read.NH) + '\n')
            readId += 1

        for read in self.spliced:
            exons = read.exons
            cigar = [str(exons[0][1] - exons[0][0]) + 'M']
            for i in xrange(1, len(exons)):
                if exons[i][0] - exons[i-1][1] == 0:
                    prevLen = int(cigar[-1][:-1])
                    cigar[-1] = str(prevLen + exons[i][1] - exons[i][0]) + 'M'
                else:

                    ####
                    if exons[i][1] == exons[i][0]:
                        print exons
                        exit()

                    cigar += [str(exons[i][0] - exons[i-1][1]) + 'N']
                    cigar += [str(exons[i][1] - exons[i][0]) + 'M']

                    if exons[i][0] < exons[i-1][1]:
                        print exons
                        exit()
            cigar = ''.join(cigar)

            chrom = read.chrom
            offset = self.chromOffsets[chrom]

            if 'N' in cigar:
                filehandle.write(chrom+':'+str(readId) + '\t0\t' + chrom + '\t' + str(exons[0][0]-offset) + '\t50\t' + cigar + '\t*\t0\t0\t*\t*\tXS:A:' + read.xs + '\tNH:i:' + str(read.NH) + '\n')
            else:
                filehandle.write(chrom+':'+str(readId) + '\t0\t' + chrom + '\t' + str(exons[0][0]-offset) + '\t50\t' + cigar + '\t*\t0\t0\t*\t*\tNH:i:' + str(read.NH) + '\n')
            readId += 1
        
        for pair in self.paired:
            exonsA = pair.exonsA
            cigarA = [str(exonsA[0][1] - exonsA[0][0]) + 'M']
            for i in xrange(1, len(exonsA)):
                if exonsA[i][0] - exonsA[i-1][1] == 0:
                    prevLen = int(cigarA[-1][:-1])
                    cigarA[-1] = str(prevLen + exonsA[i][1] - exonsA[i][0]) + 'M'
                else:
                    if exonsA[i][1] == exonsA[i][0]:
                        print exons
                        exit()

                    cigarA += [str(exonsA[i][0] - exonsA[i-1][1]) + 'N']
                    cigarA += [str(exonsA[i][1] - exonsA[i][0]) + 'M']
            cigarA = ''.join(cigarA)

            exonsB = pair.exonsB
            cigarB = [str(exonsB[0][1] - exonsB[0][0]) + 'M']
            for i in xrange(1, len(exonsB)):
                if exonsB[i][0] - exonsB[i-1][1] == 0:
                    prevLen = int(cigarB[-1][:-1])
                    cigarB[-1] = str(prevLen + exonsB[i][1] - exonsB[i][0]) + 'M'
                else:
                    ####
                    if exonsB[i][1] == exonsB[i][0]:
                        print exonsB
                        exit()

                    cigarB += [str(exonsB[i][0] - exonsB[i-1][1]) + 'N']
                    cigarB += [str(exonsB[i][1] - exonsB[i][0]) + 'M']
            cigarB = ''.join(cigarB)

            # Distance from start of first read to end of second read
            totalLen = exonsB[-1][1] - exonsA[0][0]

            chromA = pair.chromA
            chromB = pair.chromB
            offsetA = self.chromOffsets[chromA]
            if chromA == chromB:
                filehandle.write(chromA+':'+str(readId) + '\t161\t' + chromA + '\t' + str(exonsA[0][0]-offsetA) + '\t50\t' + cigarA + '\t=\t' + str(exonsB[0][0]-offsetA) + '\t' + str(totalLen) + '\t*\t*\tNH:i:' + str(pair.NH))
                if 'N' in cigarA:
                    filehandle.write('\tXS:A:' + pair.xs)
                filehandle.write('\n')

                filehandle.write(chromA+':'+str(readId) + '\t81\t' + chromA + '\t' + str(exonsB[0][0]-offsetA) + '\t50\t' + cigarB + '\t=\t' + str(exonsA[0][0]-offsetA) + '\t' + str(-totalLen) + '\t*\t*\tNH:i:' + str(pair.NH))
                if 'N' in cigarB:
                    filehandle.write('\tXS:A:' + pair.xs)
                filehandle.write('\n')
            else:
                offsetB = self.chromOffsets[chromB]
                filehandle.write(chromA+':'+str(readId) + '\t161\t' + chromA + '\t' + str(exonsA[0][0]-offsetA) + '\t50\t' + cigarA + '\t' + chromB + '\t' + str(exonsB[0][0]-offsetB) + '\t' + str(totalLen) + '\t*\t*\tNH:i:' + str(pair.NH))
                if 'N' in cigarA:
                    filehandle.write('\tXS:A:' + pair.xs)
                filehandle.write('\n')

                filehandle.write(chromA+':'+str(readId) + '\t81\t' + chromB + '\t' + str(exonsB[0][0]-offsetB) + '\t50\t' + cigarB + '\t' + chromA + '\t' + str(exonsA[0][0]-offsetA) + '\t' + str(-totalLen) + '\t*\t*\tNH:i:' + str(pair.NH))
                if 'N' in cigarB:
                    filehandle.write('\tXS:A:' + pair.xs)
                filehandle.write('\n')
            readId += 1