#! /usr/bin/env python

import bisect
import junction
import read
import os.path

class Chromosome:
    ''' A set of reads aligned to a chromosome '''

    def __init__(self, name, length):
        self.name = name
        self.length = length

        self.unspliced = []
        self.spliced = []

        self.exons = set([0, length])

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
        # Run-length encode each junction
        with open(filename, 'w') as f:

            # Compute coverage levels across every exon junction
            junctions = dict()

            readExons = []
            for read in self.spliced:
                alignment = read.exons

                # compute the length of the read
                readLen = 0

                # Find exons included in this read
                readExons = []

                for segment in alignment:
                    readLen += segment[1] - segment[0]

                    exonId = bisect.bisect_right(self.exons, segment[0])-1
                    while self.exons[exonId] < segment[1]:
                        readExons += [exonId]
                        exonId += 1

                # offset of start into first exon
                startOffset = alignment[0][0] - self.exons[readExons[0]]

                # offset of end from last exon
                endOffset =  self.exons[readExons[-1] + 1] - alignment[-1][1]

                key = '\t'.join([str(e) for e in readExons]) + '\t' + read.xs
                if not key in junctions:
                    covLength = 0
                    for e in readExons:
                        covLength += self.exons[e+1] - self.exons[e]
                    junctions[key] = junction.Junction(readExons, covLength)
                j = junctions[key]

                # update junction coverage vector in dictionary
                for i in xrange(startOffset, len(j.coverage)-endOffset):
                    j.coverage[i] += 1

                # update readLens
                if readLen in j.readLens:
                    j.readLens[readLen] += 1
                else:
                    j.readLens[readLen] = 1

            # Write exons
            f.write('\t'.join([str(e) for e in self.exons]) + '\n')


            # Write junction information
            for key, junc in junctions.items():
                # write junction information
                f.write('>' + key + '\n')

                # write read lengths
                f.write('\t'.join( [str(k)+','+str(v) for k,v in junc.readLens.items()] ) + '\n')

                # write coverage
                self.RLE(junc.coverage, f)

    def compressUnspliced(self, filename):
        ''' Compress the unspliced alignments as a run-length-encoded coverage vector

            filename: Name of file to compress to
        '''
        cov = [0] * self.length

        readLens = dict()

        for read in self.unspliced:
            alignment = read.exons
            start = alignment[0][0]
            end = alignment[0][1]

            # update coverage vector
            for base in xrange(start, end):
                cov[base] += 1

            # update read lengths distribution
            length = end - start
            if length in readLens:
                readLens[length] += 1
            else:
                readLens[length] = 1

        with open(filename, 'w') as f:
            # Write read lengths
            f.write('\t'.join([str(k)+','+str(v) for k,v in readLens.items()]) + '\n')

            # Write coverage vector
            self.RLE(cov, f)

    def expand(self, file_prefix):
        ''' Expand both spliced and unspliced alignments
        '''

        splicedName = file_prefix + '.spliced.txt'
        if os.path.isfile(splicedName):
            self.expandSpliced(splicedName)

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
                            reads = self.findReads(junc.readLens, junc.coverage)
                            for r in reads:
                                readExons = []
                                start = r[0] + self.exons[junctionExons[0]]
                                readExons.append( [start, self.exons[junctionExons[0]+1]] )

                                for i in xrange(1, len(junctionExons)-1):
                                    if junctionExons[i] == junctionExons[i-1]+1:
                                        readExons[-1][1] = self.exons[junctionExons[i]+1]
                                    else:
                                        readExons.append( [self.exons[junctionExons[i]], self.exons[junctionExons[i]+1]] )

                                end = self.exons[junctionExons[-1]+1] + r[1] - len(junc.coverage) 
                                readExons.append( [self.exons[junctionExons[-1]], end] )
                                self.spliced.append(read.Read(readExons, xs))

                        # Start of a new junction
                        key = line[1:].rstrip().split('\t')
                        junctionExons = [int(e) for e in key[:-1]]
                        xs = key[-1]
                        length = 0
                        for e in junctionExons:
                            length += self.exons[e+1] - self.exons[e]

                        junc = junction.Junction(junctionExons, length)
                        junc.readLens = None

                    else:
                        if junc.readLens == None:
                            # First line after '>' contains read length distribution
                            junc.readLens = dict()
                            for readLen in line.rstrip().split('\t'):
                                readLen = readLen.split(',')
                                junc.readLens[int(readLen[0])] = int(readLen[1])
                            junc.coverage = []
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
            reads = self.findReads(junc.readLens, junc.coverage)
            for r in reads:
                readExons = []
                start = r[0] + self.exons[junctionExons[0]]
                readExons.append( [start, self.exons[junctionExons[0]+1]] )

                for i in xrange(1, len(junctionExons)-1):
                    if junctionExons[i] == junctionExons[i-1]+1:
                        readExons[-1][1] = self.exons[junctionExons[i]+1]
                    else:
                        readExons.append( [self.exons[junctionExons[i]], self.exons[junctionExons[i]+1]] )

                end = self.exons[junctionExons[-1]+1] + r[1] - len(junc.coverage) 
                readExons.append( [self.exons[junctionExons[-1]], end] )
                self.spliced.append(read.Read(readExons, xs))



    def expandUnspliced(self, filename):
        ''' Expand a file containing compressed unspliced alignments
        '''
        
        with open(filename, "r") as f:
            readLens = None
            coverage = []

            # Read read lengths and coverage vector
            for line in f:
                # First line contains read lengths
                if readLens == None:
                    readLens = dict()
                    for row in line.rstrip().split('\t'):
                        readLen = row.split(',')
                        readLens[int(readLen[0])] = int(readLen[1])
                else:
                    row = line.rstrip().split("\t")

                    if len(row) == 1:
                        length = 1
                    else:
                        length = int(row[1])

                    coverage += [int(row[0])] * length

        reads = self.findReads(readLens, coverage)
        self.unspliced = []
        for r in reads:
            self.unspliced.append(read.Read([r]))
        

    def finalizeExons(self):
        self.exons = list(sorted(self.exons))

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

    def findReads(self, readLens, coverage):
        ''' Find the set of reads that most closely matches the distribution of readLens and the coverage vector
        '''

        debug = False

        lens = readLens.keys()

        # Find max and mode read lengths
        maxLen = max(lens)

        # Read lengths sorted by frequency, largest to smallest
        lensSorted = sorted(readLens, key=readLens.get, reverse=True)

        reads = []

        # start and end of coverage window
        # Keep finding reads from both ends until they meet in the middle
        start = 0
        while coverage[start] == 0:
            start += 1
        end = len(coverage)
        while coverage[end-1] == 0:
            end -= 1



        while end > start:
            if debug:
                print 'Coverage range: (%d, %d)' % (start, end)
                print lensSorted
                print readLens

            # find a read from the beginning
            readStart = start
            readEnd = start

            closestEndpoint = None
            for length in xrange(1, maxLen+1):
                if readStart+length < end and coverage[readStart + length] < coverage[readStart + length - 1]:

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
            while start < end and coverage[start] == 0:
                start += 1
            while end > start and coverage[end-1] == 0:
                end -= 1


            if end > start:
                if debug:
                    print 'Coverage range: (%d, %d)' % (start, end)
                    print lensSorted
                    print readLens

                # find a read from the end
                readEnd = end
                readStart = end

                closestEndpoint = None
                for length in xrange(1, maxLen+1):
                    if end-length >= start and coverage[end - length] > coverage[end - length - 1]:
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
                                #print '  j = %d' % j
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
                while coverage[end-1] == 0 and end > start:
                    end -= 1
                while coverage[start] == 0 and start < end:
                    start += 1

                
        return reads

    # def findReads(self, readLens, coverage):
    #     ''' Find the set of reads that most closely matches the distribution of readLens and the coverage vector
    #     '''

    #     # Find most common read length
    #     mode = 0
    #     modeVal = 0
    #     for k,v in readLens.items():
    #         if v > modeVal:
    #             mode = k
    #             modeVal = v

    #     # Find min and max read lengths
    #     minLen = min(readLens)
    #     maxLen = max(readLens)

    #     reads = []
    #     for i in xrange(len(coverage)):
    #         while coverage[i] > 0:
    #             # find next endpoints of reads
    #             j = i+1
    #             endpoints = []
    #             #while j < len(coverage) and j < exons[currExon+1] and (j-i+1) < (maxLen+3):
    #             while j < (len(coverage)-1) and (j-i+1) < (maxLen+3):
    #                 if coverage[j] > coverage[j+1]:
    #                     endpoints += [j-i+1] 
    #                 j += 1

    #             if len(readLens) == 0:
    #                 if len(endpoints) > 0:
    #                     length = endpoints[0]
    #                 else:
    #                     length = mode
    #             elif len(endpoints) == 0:
    #                 if mode in readLens:
    #                     length = mode
    #                 else:
    #                     length = readLens[len(readLens)/2]
    #             elif endpoints[0] < minLen:
    #                 length = minLen
    #             else:
    #                 (id1, id2) = self.findClosestVals(readLens.keys(), endpoints)
    #                 length = readLens.keys()[id1]

    #             if length in readLens:
    #                 readLens[length] -= 1
    #                 if readLens[length] == 0:
    #                     del readLens[length]
    #                     if len(readLens) > 0:
    #                         minLen = min(readLens)

    #             for x in xrange(length):
    #                 coverage[i+x] -= 1
                
    #             reads.append([i, i+length])
    #     return reads

    def insertInOrder(sortedList, a):
        ''' Insert a in the correct place in a sorted list in increasing order
        '''
        i = 0
        while i < len(sortedList) and a > sortedList[i]:
            i += 1
        return sortedList[:i] + [a] + sortedList[i:]

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
                    cigar += [str(exons[i][0] - exons[i-1][1]) + 'N']
                    cigar += [str(exons[i][1] - exons[i][0]) + 'M']
            cigar = ''.join(cigar)

            filehandle.write(self.name + str(exons[0][0]) + '\t0\t' + self.name + '\t' + str(exons[0][0]) + '\t50\t' + cigar + '\t*\t0\t0\t*\t*\tXS:A:' + read.xs + '\n')
            #filehandle.write('Name\t0\t' + self.name + '\t' + str(exons[0][0]) + '\t50\t' + cigar + '\t*\t0\t0\t*\t*\tXS:A:' + read.xs + '\n')
