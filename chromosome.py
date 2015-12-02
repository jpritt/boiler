#! /usr/bin/env python

import bisect
import bucket
import read
import pairedread
import os.path
import sys

class Chromosome:
    ''' A set of reads aligned to a chromosome '''

    def __init__(self, name, length):
        self.name = name
        self.length = length

        self.exons = set([0, length])

        self.unspliced = []
        self.spliced = []

        # paired reads for which the mate still needs to be found
        self.unmatched = []
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

        if len(self.spliced) > 0:
            self.compressSpliced(file_prefix + '.spliced.txt')

        if len(self.unspliced) > 0:
            self.compressUnspliced(file_prefix + '.unspliced.txt')

    def compressSpliced(self, filename):
        ''' Compress the spliced alignments to a single file

            filename: Name of file to compress to
        '''
            
        # Compute coverage levels across every exon junction
        junctions = dict()

        for r in self.spliced:
            exonIds = r.exonIds

            '''
            if self.name == '3L':
                for n in xrange(1,len(exonIds)):
                    if self.exons[exonIds[n]] == 7256049:
                        print 'Exons: ' + str(r.exons)
                        print 'Exon ids: ' + str(exonIds)
                        print 'Junction exons: ' + str([[self.exons[i], self.exons[i+1]] for i in r.exonIds])
                        print 'Left, right = %d, %d' % (r.lenLeft, r.lenRight)
            '''

            if not r.xs == None:
                # XS is defined for this read
                key = '\t'.join([str(e) for e in exonIds]) + '\t' + r.xs
                if not key in junctions:
                    covLength = 0
                    for e in exonIds:
                        covLength += self.exons[e+1] - self.exons[e]
                    junctions[key] = bucket.Junction(exonIds, covLength)
                j = junctions[key]
            else:
                # XS not defined for this read, so see which one (+/-) is in the dict already or add + by default
                key1 = '\t'.join([str(e) for e in exonIds]) + '\t+'
                key2 = '\t'.join([str(e) for e in exonIds]) + '\t-'
                if key1 in junctions:
                    read.xs = '+'
                    j = junctions[key1]
                elif key2 in junctions:
                    read.xs = '-'
                    j = junctions[key2]
                else:
                    read.xs = '+'
                    covLength = 0
                    for e in exonIds:
                        covLength += self.exons[e+1] - self.exons[e]
                    junctions[key1] = bucket.Junction(exonIds, covLength)
                    j = junctions[key1]

            # update junction coverage vector in dictionary
            for i in xrange(r.startOffset, len(j.coverage)-r.endOffset):
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

        with open(filename, 'w') as f:

            # Write exons
            f.write('\t'.join([str(e) for e in self.exons]) + '\n')


            # Write junction information
            for key, junc in junctions.items():
                '''
                if self.name == '3L':
                    printDetails = False
                    for e in junc.exons:
                        if self.exons[e] > 7365300 and self.exons[e] < 7365900:
                            printDetails = True

                    if printDetails:
                        juncBounds = []
                        for e in junc.exons:
                            juncBounds.append([self.exons[e], self.exons[e+1]])
                        print 'Junction: ' + str(juncBounds)
                        print 'Coverage: ' + str(junc.coverage)
                        print 'Read lens: ' + str(junc.readLens)
                        print 'Lens left: ' + str(junc.lensLeft)
                        print 'Lens right: ' + str(junc.lensRight)
                        print ''
                '''

                # write junction information
                f.write('>' + key + '\n')

                # write read lengths
                f.write('\t'.join( [str(k)+','+str(v) for k,v in junc.readLens.items()] ) + '\n')

                # Write left and right read lengths
                f.write('\t'.join([str(k)+','+str(v) for k,v in junc.lensLeft.items()]) + '\n')
                f.write('\t'.join([str(k)+','+str(v) for k,v in junc.lensRight.items()]) + '\n')

                # write coverage
                self.RLE(junc.coverage, f)

    def compressUnspliced(self, filename):
        ''' Compress the unspliced alignments as a run-length-encoded coverage vector

            filename: Name of file to compress to
        '''

        # sort reads into exons
        readExons = []
        for i in xrange(len(self.exons)-1):
            readExons += [[]]
        for r in self.unspliced:
            j = 0
            while r.exons[0][0] >= self.exons[j+1]:
                j += 1
            readExons[j] += [r]

        with open(filename, 'w') as f:
            for i in xrange(len(readExons)):
                length = self.exons[i+1] - self.exons[i]

                exonStart = self.exons[i]

                cov = [0] * length

                # Distribution of all read lengths
                readLens = dict()

                # Distribution of all gap lengths in paired-end reads
                lensLeft = dict()
                lensRight = dict()

                for read in readExons[i]:
                    alignment = read.exons

                    start = alignment[0][0] - exonStart
                    end = alignment[0][1] - exonStart

                    # update coverage vector
                    for base in xrange(start, end):
                        cov[base] += 1

                    # update read lengths distribution
                    length = end - start
                    if length in readLens:
                        readLens[length] += 1
                    else:
                        readLens[length] = 1

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
                if self.name == '3L' and exonStart == 7256049:
                    print 'Segment Num: ' + str(i)
                    print 'Exon start: ' + str(self.exons[i])
                    #print 'Exons: ' + self.exons
                    print 'Reads: ' + '\t'.join([str(r.exons) for r in readExons[i]])
                    print 'Coverage: '
                    self.RLE(cov, sys.stdout)
                    print 'Read lengths: ' + str(readLens)
                    #exit()
                '''

                # Write read lengths
                f.write('>' + '\t'.join([str(k)+','+str(v) for k,v in readLens.items()]) + '\n')

                # Write left and right read lengths
                f.write('\t'.join([str(k)+','+str(v) for k,v in lensLeft.items()]) + '\n')
                f.write('\t'.join([str(k)+','+str(v) for k,v in lensRight.items()]) + '\n')

                # Write coverage vector
                self.RLE(cov, f)

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

                            '''
                            if self.name == '3L':
                                printDetails = False
                                for e in junc.exons:
                                    if self.exons[e] > 7365300 and self.exons[e] < 7365900:
                                        printDetails = True

                                if printDetails:
                                    juncBounds = []
                                    for e in junc.exons:
                                        juncBounds.append([self.exons[e], self.exons[e+1]])
                                    print 'Junction: ' + str(juncBounds)
                                    print 'Coverage: ' + str(junc.coverage)
                                    print 'Read lens: ' + str(junc.readLens)
                                    print 'Lens left: ' + str(junc.lensLeft)
                                    print 'Lens right: ' + str(junc.lensRight)
                                    print ''
                            '''



                            #print junc.exons
                            unpaired, paired = self.findReads(junc.readLens, junc.lensLeft, junc.lensRight, junc.coverage)

                            '''
                            if self.name == '3L':
                                if len(junc.exons) == 4 and junc.exons[-1] == 2081:
                                    print 'Junction exons: ' + str(junc.exons)
                                    print 'Read lengths: ' + str(junc.readLens)

                                    for r in unpaired:
                                        print 'Unpaired read: ' + str(r)
                                    for p in paired:
                                        print 'Paired read: ' + str(p[0]) + ', ' + str(p[1])
                            '''

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
                                #offsets.append(mapping[-1] + juncBounds[j][0] - juncBounds[j-1][1])
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

                                self.spliced.append(read.Read(readExons, xs))

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
                                
                                '''
                                if self.name == '3L':
                                    if len(junc.exons) == 4 and junc.exons[-1] == 2081:
                                        print 'Exons A: ' + str(readExonsA)
                                        print 'Exons B: ' + str(readExonsB)
                                '''

                                self.paired.append(pairedread.PairedRead(readExonsA, readExonsB))

                        # Start of a new junction
                        key = line[1:].rstrip().split('\t')
                        junctionExons = [int(e) for e in key[:-1]]
                        xs = key[-1]
                        length = 0
                        for e in junctionExons:
                            length += self.exons[e+1] - self.exons[e]

                        junc = bucket.Junction(junctionExons, length)
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

            '''
            if self.name == '3L':
                printDetails = False
                for e in junc.exons:
                    if self.exons[e] > 7365300 and self.exons[e] < 7365900:
                        printDetails = True

                if printDetails:
                    juncBounds = []
                    for e in junc.exons:
                        juncBounds.append([self.exons[e], self.exons[e+1]])
                    print 'Junction: ' + str(juncBounds)
                    print 'Coverage: ' + str(junc.coverage)
                    print 'Read lens: ' + str(junc.readLens)
                    print 'Lens left: ' + str(junc.lensLeft)
                    print 'Lens right: ' + str(junc.lensRight)
                    print ''
            '''

            # process the final junction
            unpaired, paired = self.findReads(junc.readLens, junc.lensLeft, junc.lensRight, junc.coverage)

            '''
            if self.name == '3L':
                if len(junc.exons) == 4 and junc.exons[-1] == 2081:
                    print 'Junction exons: ' + str(junc.exons)
                    print 'Read lengths: ' + str(junc.readLens)

                    for r in unpaired:
                        print 'Unpaired read: ' + r.exons
                    for p in paired:
                        print 'Paired read: ' + p.exonsA + ', ' + p.exonsB

            #print 'Expanded:'
            #print paired
            #print ''
            '''

            juncBounds = []
            for j in junctionExons:
                juncBounds.append([self.exons[j], self.exons[j+1]])

            #print 'Junction: ' + str(junc.exons)
            #print 'Bounds: ' + str(juncBounds)

            # Offset of start of junction from beginning of chromosome
            juncOffset = juncBounds[0][0]
            
            # marks indices in junction coverage vector where exons are
            mapping = [0]
            # offsets of start of each exon in junction
            offsets = [0]
            for j in xrange(1,len(juncBounds)):
                mapping.append(mapping[-1] + juncBounds[j-1][1] - juncBounds[j-1][0])
                #offsets.append(mapping[-1] + juncBounds[j][0] - juncBounds[j-1][1])
                offsets.append(juncBounds[j][0] - juncBounds[0][0])
            #print 'Mapping: ' + str(mapping)
            #print 'Offsets: ' + str(offsets)
            #print ''

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

                self.spliced.append(read.Read(readExons, xs))

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

                self.paired.append(pairedread.PairedRead(readExonsA, readExonsB))

    def expandUnspliced(self, filename):
        ''' Expand a file containing compressed unspliced alignments
        '''
        #print self.exons[0]
        
        '''
        with open(filename, 'r') as f:
            firstLine = True
            for line in f:
                row = [int(r) for r in line.rstrip().split('\t')]
                if firstLine:
                    start = row[0]
                    end = row[1]
                else:
                    lenLeft = int(row[0])
                    lenRight = int(row[1])

                    if lenLeft > 0 or lenRight > 0:
                        self.paired.append(pairedread.PairedRead([[start, start+lenLeft]], [[end-lenRight, end]]))
                    else:
                        self.unspliced.append(read.Read([[start, end]]))

                firstLine = not firstLine

        return
        '''

        with open(filename, "r") as f:
            lineNum = 0
            segmentNum = 0
            readLens = None
            lensLeft = None
            lensRight = None
            coverage = []

            # Read read lengths and coverage vector
            for line in f:
                if line[0] == '>':

                    '''
                    #print 'Segment num: ' + str(segmentNum)
                    if self.name == '3L' and self.exons[segmentNum] == 7256049:
                        print 'Segment num: ' + str(segmentNum)
                        print 'Coverage: '
                        self.RLE(coverage, sys.stdout)
                        print 'Read lengths: ' + str(readLens)
                        exit()
                    elif self.name == '3L' and segmentNum == 2081:
                        print 'Start: ' + str(start)
                        print 'Coverage: '
                        self.RLE(coverage, sys.stdout)
                        print 'Read lengths: ' + str(readLens)
                        exit()
                    '''

                    if not readLens == None and len(readLens) > 0:
                        '''
                        if self.name == '3L' and self.exons[segmentNum-1] == 7256049:
                            print 'Segment Num: ' + str(segmentNum-1)
                            print 'Exon start: ' + str(self.exons[segmentNum-1])
                            print 'Coverage: '
                            self.RLE(coverage, sys.stdout)
                            print 'Read lengths: ' + str(readLens)
                            exit()
                        '''

                        unpaired, paired = self.findReads(readLens, lensLeft, lensRight, coverage)
                        start = self.exons[segmentNum-1]

                        

                        for r in unpaired:
                            self.unspliced.append(read.Read([[r[0]+start, r[1]+start]]))

                        for p in paired:
                            self.paired.append(pairedread.PairedRead([[p[0][0]+start, p[0][1]+start]], [[p[1][0]+start, p[1][1]+start]]))

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


            unpaired, paired = self.findReads(readLens, lensLeft, lensRight, coverage)
            start = self.exons[segmentNum]

            for r in unpaired:
                self.unspliced.append(read.Read([[r[0]+start, r[1]+start]]))

            for p in paired:
                self.paired.append(pairedread.PairedRead([[p[0][0]+start, p[0][1]+start]], [[p[1][0]+start, p[1][1]+start]]))

    def finalizeExons(self):
        ''' Convert the set of exon boundaries to a list
        '''
        self.exons = list(sorted(self.exons))
        #self.exons = list(self.exons)

    def finalizeReads(self):
        ''' Now that all exon boundaries are known, fix unspliced regions that cross exon boundaries
            and finalized paired-end reads
        '''

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

            ## Fix this problem with exons not matching
            if n > 0 and not pair.exonIdsA[-n:] == pair.exonIdsB[:n]:
                continue

            # TODO: Replace gap with end lengths
            if n > 0:
                newRead = read.Read(exonsA[:-1] + [[exonsA[-1][0], exonsB[n-1][1]]] + exonsB[n:], pair.xs)
                gap = 0
                for i in xrange(n):
                    gap += exonsB[i][0] - exonsA[i-n][1]
                #gap = exonsB[n-1][0] - exonsA[-1][1]
            else:
                newRead = read.Read(exonsA+exonsB, pair.xs)
                gap = exonsB[0][0] - self.exons[pair.exonIdsB[0]] + self.exons[pair.exonIdsA[-1]+1] - exonsA[-1][1]
            newRead.exonIds = pair.exonIdsA + pair.exonIdsB[n:]

            newRead.lenLeft = readLenLeft
            newRead.lenRight = readLenRight
            newRead.readLen = readLenLeft + readLenRight + gap

            # offset of start into first exon
            newRead.startOffset = exonsA[0][0] - self.exons[pair.exonIdsA[0]]

            # offset of end from last exon
            newRead.endOffset = self.exons[pair.exonIdsB[-1] + 1] - exonsB[-1][1]

            '''
            if newRead.startOffset == 131 and newRead.endOffset == 16079:
                print ''
                print 'ExonsA: ' + str(exonsA)
                print 'Exon ids A: ' + str(pair.exonIdsA)
                juncBounds = []
                for i in pair.exonIdsA:
                    juncBounds.append([self.exons[i], self.exons[i+1]])
                print 'Junction bounds: ' + str(juncBounds)
                print 'Length left = ' + str(readLenLeft)
                print ''

                print 'ExonsB: ' + str(exonsB)
                print 'Exon ids B: ' + str(pair.exonIdsB)
                juncBounds = []
                for i in pair.exonIdsB:
                    juncBounds.append([self.exons[i], self.exons[i+1]])
                print 'Junction bounds: ' + str(juncBounds)
                print 'Length right = ' + str(readLenRight)
                print ''

                print 'Gap = ' + str(gap)
                print 'Read length = ' + str(newRead.readLen)
                print ''
                exit()
            '''

            if len(newRead.exonIds) == 1:
                self.unspliced += [newRead]
            else:
                self.spliced += [newRead]
        #print ''
        #print ''



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

    def findReads(self, readLens, lensLeft, lensRight, coverage, boundaries=None):
        ''' Find the set of reads that most closely matches the distribution of readLens and the coverage vector
        '''

        countLeft = 0
        for value in lensLeft.values():
            countLeft += value
        countRight = 0
        for value in lensRight.values():
            countRight += value
        if not countLeft == countRight:
            print 'Error! Left and right read dists not the same size!'

        '''
        rle = []
        val = coverage[0]
        length = 0

        for v in coverage:
            if v == val:
                length += 1
            else:
                if length == 1:
                    rle += [str(val)]
                else:
                    rle += [str(val) + 'x' + str(length)]
                val = v
                length = 1
        if length == 1:
            rle += [str(val)]
        else:
            rle += [str(val) + 'x' + str(length)]
        '''


        #countReads = 0
        #for r in readLens.values():
        #    countReads += r

        totalReads = 0
        totalLeft = 0
        totalRight = 0
        for v in readLens.values():
            totalReads += v
        for v in lensLeft.values():
            totalLeft += v
        for v in lensRight.values():
            totalRight += v

        #print 'Read lengths (' + str(totalReads) + '): ' + str(readLens)
        #self.RLE(coverage, sys.stdout)
        #print 'Lengths left (' + str(totalLeft) + ') : ' + str(lensLeft)
        #print 'Lengths right (' + str(totalRight) + '): ' + str(lensRight)


        #reads = self.findReadsInCoverage_v1(coverage, readLens, boundaries)
        reads = self.findReadsInCoverage_v2(coverage, readLens, boundaries)
        #reads = self.findReadsInCoverage_v3(coverage, readLens, boundaries)


        #print '%d reads found: ' % len(reads)
        #print ''

        if not len(reads) == totalReads:
            print '%d / %d found' % (len(reads), totalReads)

        if reads == None:
            print 'Error finding reads!'
            exit()

        # Pair left and right read lengths with long fragments
        # First sort reads and lengths by decreasing length
        reads.sort(key=lambda r: r[1]-r[0], reverse=True)
        lensLeftSorted = sorted(lensLeft.keys(), reverse=True)
        lensRightSorted = sorted(lensRight.keys(), reverse=True)

        pairedReads = []

        readId = 0

        if len(lensLeftSorted) > len(reads):
            print 

        while len(lensLeftSorted) > 0:
            if readId >= len(reads):
                print 'Error! More left/right lengths than fragments!'
                exit()

            start = reads[readId][0]
            end = reads[readId][1]
            fragmentLen = end-start

            lenLeftId = 0
            lenLeft = lensLeftSorted[lenLeftId]
            lenRightId = 0
            lenRight = lensRightSorted[lenRightId]


            pairedReads.append( [[start, start+lenLeft], [end-lenRight, end]] )

            lensLeft[lenLeft] -= 1
            if lensLeft[lenLeft] == 0:
                del lensLeftSorted[lenLeftId]

            lensRight[lenRight] -= 1
            if lensRight[lenRight] == 0:
                del lensRightSorted[lenRightId]

            readId += 1

        return reads[readId:], pairedReads

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


            closestEndpoint = None
            for length in xrange(1, min(maxLen+1, boundaries[boundBottom]-start)):
            #for length in xrange(1, maxLen+1):
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
            if readEnd == readStart:
                if closestEndpoint == None:
                    length = lensSorted[0]
                    readEnd = readStart + length
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
                #for length in xrange(1, maxLen+1):
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
        #lensSorted = sorted(readLens, key=readLens.get, reverse=True)
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

        #print 'Coverage: ' + ','.join([str(c) for c in coverage[start:end]])
        #print 'Read lengths: ' + str(readLens)
        #print 'Lengths sorted: ' + str(lensSorted)

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

                    #if readLens[length] == 0:
                    #    del readLens[length]

                    readFound = True
                    break
            if not readFound:
                #print '  No good endpoint found'
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

                #if readLens[length] == 0:
                #    del readLens[length]
            #print '  Adding read (%d, %d)' % (readStart, readEnd)
            #print '  Lengths sorted: ' + str(lensSorted)
            #print '  Read Lengths: ' + str(readLens)

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

                        #if readLens[length] == 0:
                        #    del readLens[length]

                        readFound = True
                        break

                if not readFound:
                    #print '  No good endpoint found'
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

                    #if readLens[length] == 0:
                    #    del readLens[length]
                #print '  Adding read (%d, %d)' % (readStart, readEnd)
                #print '  Lengths sorted: ' + str(lensSorted)
                #print '  Read Lengths: ' + str(readLens)

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
        #self.RLE(coverage, sys.stdout)

        '''
        totalLen = 0
        for r,num in readLens.items():
            totalLen += r*num
        covLen = sum(coverage)
        if not totalLen == sum(coverage):
            print 'Error! %d =/= %d' % (covLen, totalLen)
            exit()
        '''

        length = end-start

        firstLength = readLens.keys()[0]
        # Only 1 read left: return that read iff coverage[start:end] is all 1s of the correct length
        if len(readLens) == 1 and readLens[firstLength] == 1:
            #print '  '*level + 'Coverage: ' + ','.join([str(c) for c in coverage])
            #print '  '*level + 'Read lengths: ' + str(readLens)
            if length == firstLength:
                for i in xrange(start,end):
                    if not coverage[i] == 1:
                        return None
                return [[start,end]]

        # More than 1 read left: pick 1 recursively
        iteration = 0
        uniqueReads = readLens.keys()
        for r in uniqueReads:
            #print '  '*level + 'Testing ' + str(r)
            #print '  '*level + 'Coverage: ' + ','.join([str(c) for c in coverage])
            #print '  '*level + 'Read lengths: ' + str(readLens)
            #print '  '*level + 'Trying read ' + str(r)


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
        
                    '''
                    if reads == None:
                        #print '  '*level + 'Doesn\'t work'
                        # Undo all the changes above

                        # Add read back to coverage
                        for i in xrange(r):
                            coverage[i] += 1

                        if r in readLens:
                            readLens[r] += 1
                        else:
                            readLens[r] = 1
                    '''
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
        

        
        #if level == 5:
        #    for i in xrange(start,end):
        #        if coverage[i] == 0:
        #            return None
        #    return [[start,end]]
        

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


    def insertInOrder(sortedList, a):
        ''' Insert a in the correct place in a sorted list in increasing order
        '''
        i = 0
        while i < len(sortedList) and a > sortedList[i]:
            i += 1
        return sortedList[:i] + [a] + sortedList[i:]

    def processRead(self, read, pair, pair_offset):
        ''' If read is unpaired, add it to the correct spliced or unspliced list of reads.
            If read is paired, find its pair or add it to a list to be found later. Once a pair of reads is found, add the combined read to the appropriate list of reads
        '''
        #print 'Processing read [' + '  '.join([str(e[0])+','+str(e[1]) for e in read.exons]) + '], ' + str(pair) + ', ' + str(pair_offset)

        if pair == 0 and pair_offset == 0:
            # unpaired read
            if len(read.exons) == 1:
                #print '  Adding unspliced'
                self.addUnspliced(read)
            else:
                #print '  Adding spliced'
                self.addSpliced(read)

                # update list of exons
                alignment = read.exons
                if len(alignment) > 1:
                    for i in xrange(len(alignment)-1):
                        self.exons.add(alignment[i][1])
                        self.exons.add(alignment[i+1][0])
        else:
            # TODO: Use binary search here for speed
            matched = False

            debug = False
            #for r in read.exons:
            #    if r[0] == 

            i = 0
            while i < len(self.unmatched):
                r = self.unmatched[i]
                if r[0] == pair and r[1] == read.exons[0][0]:
                    matched = True

                    # create new read
                    match = r[2]
                    del self.unmatched[i]

                    # Tophat alignments should be in increasing order of index
                    if match.exons[0][0] > read.exons[0][0]:
                        print 'Error! Tophat file not sorted!'
                        exit()

                    #newExons = match.exons[:-1] + [[match.exons[-1][0], read.exons[0][1]]] + read.exons[1:]

                    xs = read.xs or match.xs
                    #print '  Matched to ' + '  '.join([str(e[0])+','+str(e[1]) for e in match.exons])


                    # update list of exons
                    alignment = match.exons
                    if len(alignment) > 1:
                        for i in xrange(len(alignment)-1):
                            self.exons.add(alignment[i][1])
                            self.exons.add(alignment[i+1][0])

                    alignment = read.exons
                    if len(alignment) > 1:
                        for i in xrange(len(alignment)-1):
                            self.exons.add(alignment[i][1])
                            self.exons.add(alignment[i+1][0])

                    self.paired.append(pairedread.PairedRead(match.exons, read.exons, xs))

                    break
                else:
                    i += 1
            if not matched:
                #print '  No match found'
                self.unmatched += [(read.exons[0][0], pair, read)]
        #print ''

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

        for read in self.unspliced:
            exons = read.exons
            filehandle.write(self.name + str(exons[0][0]) + '\t0\t' + self.name + '\t' + str(exons[0][0]) + '\t50\t' + str(exons[0][1]-exons[0][0]) + 'M\t*\t0\t0\t*\t*\n')
            #filehandle.write('Name\t0\t' + self.name + '\t' + str(exons[0][0]) + '\t50\t' + str(exons[0][1]-exons[0][0]) + 'M\t*\t0\t0\t*\t*\n')

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

            filehandle.write(self.name + str(exons[0][0]) + '\t0\t' + self.name + '\t' + str(exons[0][0]) + '\t50\t' + cigar + '\t*\t0\t0\t*\t*\tXS:A:' + read.xs + '\n')
        
        for pair in self.paired:
            exonsA = pair.exonsA
            cigarA = [str(exonsA[0][1] - exonsA[0][0]) + 'M']
            for i in xrange(1, len(exonsA)):
                if exonsA[i][0] - exonsA[i-1][1] == 0:
                    prevLen = int(cigarA[-1][:-1])
                    cigarA[-1] = str(prevLen + exonsA[i][1] - exonsA[i][0]) + 'M'
                else:
                    ####
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

            filehandle.write(self.name + str(exonsA[0][0]) + '\t0\t' + self.name + '\t' + str(exonsA[0][0]) + '\t50\t' + cigarA + '\t=\t' + str(exonsB[0][0]) + '\t' + str(totalLen) + '\t*\t*\tXS:A:' + read.xs + '\n')
            filehandle.write(self.name + str(exonsB[0][0]) + '\t0\t' + self.name + '\t' + str(exonsB[0][0]) + '\t50\t' + cigarB + '\t=\t' + str(exonsA[0][0]) + '\t' + str(-totalLen) + '\t*\t*\tXS:A:' + read.xs + '\n')
            
            '''
            if not pair.xs == None:
                filehandle.write(self.name + str(exonsA[0][0]) + '\t0\t' + self.name + '\t' + str(exonsA[0][0]) + '\t50\t' + cigarA + '\t=\t' + str(exonsB[0][0]) + '\t' + str(totalLen) + '\t*\t*\tXS:A:' + pair.xs + '\n')
                filehandle.write(self.name + str(exonsB[0][0]) + '\t0\t' + self.name + '\t' + str(exonsB[0][0]) + '\t50\t' + cigarB + '\t=\t' + str(exonsA[0][0]) + '\t' + str(-totalLen) + '\t*\t*\tXS:A:' + pair.xs + '\n')
            else:
                if len(exonsA) > 1:
                    filehandle.write(self.name + str(exonsA[0][0]) + '\t0\t' + self.name + '\t' + str(exonsA[0][0]) + '\t50\t' + cigarA + '\t=\t' + str(exonsB[0][0]) + '\t' + str(totalLen) + '\t*\t*\tXS:A:+\n')
                else:
                    filehandle.write(self.name + str(exonsA[0][0]) + '\t0\t' + self.name + '\t' + str(exonsA[0][0]) + '\t50\t' + cigarA + '\t=\t' + str(exonsB[0][0]) + '\t' + str(totalLen) + '\t*\t*\n')

                if len(exonsB) > 1:
                    filehandle.write(self.name + str(exonsB[0][0]) + '\t0\t' + self.name + '\t' + str(exonsB[0][0]) + '\t50\t' + cigarB + '\t=\t' + str(exonsA[0][0]) + '\t' + str(-totalLen) + '\t*\t*\tXS:A:+\n')
                else:
                    filehandle.write(self.name + str(exonsB[0][0]) + '\t0\t' + self.name + '\t' + str(exonsB[0][0]) + '\t50\t' + cigarB + '\t=\t' + str(exonsA[0][0]) + '\t' + str(-totalLen) + '\t*\t*\n')
            '''