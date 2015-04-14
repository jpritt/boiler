#! /usr/bin/env python

import bisect
import junction
import read
import pairedread
import os.path
import sys

import objgraph
import random

import time

class Alignments:
    ''' A set of reads aligned to a genome '''

    def __init__(self, chromosomes):
        ''' Initialize a genome for alignments

            chromosomes: A dictionary with keys corresponding to chromosome names and values corresponding to lengths
        '''
        self.chromosomes = chromosomes
        self.chromosomeNames = sorted(chromosomes.keys())

        # Initialize exon breaks between all chromosomes
        self.exons = [0]

        # Offset of each chromosome from the start of the genome
        self.chromOffsets = dict()

        nextOffset = 0
        for c in self.chromosomeNames:
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
        for i in range(len(alignment)-1):
            self.exons.add(alignment[i][1])
            self.exons.add(alignment[i+1][0])

    def finalizeExons(self):
        ''' Convert the set of exon boundaries to a list
        '''

        self.exons = list(sorted(self.exons))

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
                for i in range(1, len(spliceSites)):
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
                    for i in range(1, len(spliceSites)):
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
                    for i in range(1, len(spliceSites)):
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
                for i in range(n):
                    gap += exonsB[i][0] - exonsA[i-n][1]
                #gap = exonsB[n-1][0] - exonsA[-1][1]
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
        del self.paired

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
        for i in range(len1):
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

    def findPairsSemiRecursive(self, starts, ends, reads, readLens, readLensSorted, i, currentScore=0, depthCutoff=None):
        maxDepth = 10

        unpaired = []
        paired = []
        score = 0

        while i < len(starts):
            u, p, s = self.findPairsRecursive(starts, ends, reads, readLens, readLensSorted, i, currentScore, 0, depthCutoff, maxDepth)
            unpaired += u
            paired += p
            score += s

            while i < len(starts) and starts[i] == 0:
                i += 1

        return unpaired, paired, score


    def findPairsRecursive(self, starts, ends, reads, readLens, readLensSorted, i, currentScore=0, currentDepth=0, depthCutoff=None, maxDepth=None):
        '''
            starts: array with ones wherever reads start
            ends: array with ones wherever reads end
            reads: fragments extracted from coverage vector
            readLens: lengths of paired and unpaired reads
        '''

        if depthCutoff == None:
            depthCutoff = len(reads)
        
        if not (maxDepth == None) and currentDepth >= maxDepth:
            return [],[],0

        if i == len(starts):
            return [], [], 0

        if starts[i] == 0:
            while i < len(starts) and starts[i] == 0:
                i += 1
            return self.findPairsRecursive(starts, ends, reads, readLens, readLensSorted, i, currentScore, currentDepth, depthCutoff, maxDepth)

        id1 = 0
        while not reads[id1][0] == i:
            id1 += 1

        bestScore = None
        bestRead = None
        bestPaired = None
        bestUnpaired = None
        bestSingle = None
        for k in range(len(readLensSorted)):
            l = readLensSorted[k]
            if ends[i+l] > 0:
                for id2 in range(len(reads)):
                    if reads[id2][1] == i+l:
                        # update reads, readLens, starts
                        if id1 == id2:
                            read = reads[id1]
                            single = True

                            del reads[id1]

                            # update starts and ends
                            starts[read[0]] -= 1
                            ends[read[1]] -= 1
                        else:
                            read = [reads[id1], reads[id2]]
                            single = False

                            if id1 > id2:
                                del reads[id1]
                                del reads[id2]
                            else:
                                del reads[id2]
                                del reads[id1]

                            starts[read[0][0]] -= 1
                            ends[read[0][1]] -= 1
                            starts[read[1][0]] -= 1
                            ends[read[1][1]] -= 1

                        # update readLens
                        readLens[l] -= 1

                        readLensSortedNew = readLensSorted[:]
                        if readLens[l] <= 0:
                            del readLensSortedNew[k]

                        # recurse
                        newi = i
                        while newi < len(starts) and starts[newi] == 0:
                            newi += 1
                        unpaired, paired, score = self.findPairsRecursive(starts, ends, reads, readLens, readLensSortedNew, newi, currentScore, currentDepth+1, depthCutoff, maxDepth)

                        # add read back
                        if single:
                            reads.append(read)

                            starts[read[0]] += 1
                            ends[read[1]] += 1
                        else:
                            reads.append(read[0])
                            reads.append(read[1])

                            starts[read[0][0]] += 1
                            ends[read[0][1]] += 1
                            starts[read[1][0]] += 1
                            ends[read[1][1]] += 1

                        readLens[l] += 1

                        if score == 0:
                            if single:
                                return unpaired+[read], paired, 0
                            else:
                                return unpaired, paired+[read], 0
                        elif score < bestScore:
                            bestRead = read
                            bestScore = score
                            bestPaired = paired
                            bestUnpaired = unpaired
                            bestSingle = single

                        if score < depthCutoff:
                            depthCutoff = score
        if bestScore == None:
            if currentScore >= depthCutoff:
                return [],[],1

            read = reads[id1]
            del reads[id1]

            # update starts and ends
            starts[read[0]] -= 1
            ends[read[1]] -= 1

            newi = i
            while newi < len(starts) and starts[newi] == 0:
                newi += 1
            unpaired, paired, score = self.findPairsRecursive(starts, ends, reads, readLens, readLensSorted, newi, currentScore+1, currentDepth+1, depthCutoff, maxDepth)


            reads.append(read)
            starts[read[0]] += 1
            ends[read[1]] += 1
            return unpaired, paired, score+1
        else:
            if bestSingle:
                return bestUnpaired+[bestRead], bestPaired, bestScore
            else:
                return bestUnpaired, bestPaired+[bestRead], bestScore



    def findPairs2(self, length, reads, readLens):#pairedLens, unpairedLens):
        # Sort read lengths from largest to smallest
        readLensSorted = sorted(readLens, reverse=True)

        starts = [0] * (length+1)
        ends = [0] * (length+1)

        for r in reads:
            starts[r[0]] += 1
            ends[r[1]] += 1

        paired = []
        unpaired = []

        unmatched = []

        i = 0
        while i < length:
            if starts[i] > 0:
                id1 = 0
                while not reads[id1][0] == i:
                    id1 += 1

                startRead = reads[id1]

                starts[i] -= 1
                ends[reads[id1][1]] -= 1
                del reads[id1]

                foundRead = False
                for j in range(len(readLensSorted)):
                    l = readLensSorted[j]
                    if i+l <= length and ends[i+l] > 0:
                        readLens[l] -= 1
                        if readLens[l] == 0:
                            del readLensSorted[j]
                        foundRead = True

                        # Add paired read to list
                        id2 = 0
                        while not reads[id2][1] == i+l:
                            id2 += 1

                        starts[reads[id2][0]] -= 1
                        ends[i+l] -= 1

                        paired += [[startRead, reads[id2]]]
                        del reads[id2]

                        break
                
                if not foundRead:
                    readLen = startRead[1]-startRead[0]
                    if readLen in readLensSorted:
                        unpaired += [startRead]

                        readLens[readLen] -= 1
                        if readLens[readLen] == 0:
                            for j in range(len(readLensSorted)):
                                if readLensSorted[j] == readLen:
                                    del readLensSorted[j]
                                    break
                    else:
                        unmatched += [startRead]
                
            if starts[i] <= 0:
                i += 1

        if len(unmatched) > 0:
            countReadLens = 0
            for k,v in readLens.items():
                countReadLens += v

            while len(unmatched) > max(1, countReadLens):
                #paired += [[unmatched[0], unmatched[len(unmatched)/2]]]
                paired += [[unmatched[0], unmatched[-1]]]
                #del unmatched[len(unmatched)/2]
                del unmatched[-1]
                del unmatched[0]
                countReadLens -= 1

            for r in unmatched:
                unpaired += [r]

        return unpaired, paired

    def findPairs(self, length, reads, pairedLens, unpairedLens):
        paired = []
        unpaired = []

        if len(pairedLens) > 0:
            # Sort read pairedLensSortedlengths from largest to smallest
            pairedLensSorted = sorted(pairedLens, reverse=True)

            starts = [0] * (length+1)
            ends = [0] * (length+1)


            for r in reads:
                starts[r[0]] += 1
                ends[r[1]] += 1

            unmatched = []

            i = 0
            j = length

            while i < j and starts[i] <= 0:
                i += 1
            while j > i and ends[j] <= 0:
                j -= 1 

            while i < j and len(pairedLensSorted) > 0:
                id1 = 0
                while not reads[id1][0] == i:
                    id1 += 1

                startRead = reads[id1]

                starts[i] -= 1
                ends[reads[id1][1]] -= 1
                del reads[id1]

                foundRead = False
                for r in range(len(pairedLensSorted)):
                    l = pairedLensSorted[r]
                    if i+l <= length and ends[i+l] > 0:
                        pairedLens[l] -= 1
                        if pairedLens[l] == 0:
                            del pairedLensSorted[r]
                        foundRead = True

                        # Add paired read to list
                        id2 = 0
                        while not reads[id2][1] == i+l:
                            id2 += 1

                        starts[reads[id2][0]] -= 1
                        ends[i+l] -= 1

                        paired += [[startRead, reads[id2]]]
                        del reads[id2]

                        break
                
                if not foundRead:
                    readLen = startRead[1]-startRead[0]
                    if readLen in unpairedLens:
                        unpaired += [startRead]

                        if unpairedLens[readLen] == 1:
                            del unpairedLens[readLen]
                        else:
                            unpairedLens[readLen] -= 1
                    else:
                        unmatched += [startRead]
                    
                while i < j and starts[i] <= 0:
                    i += 1
                while j > i and ends[j] <= 0:
                    j -= 1

                if j > i and len(pairedLensSorted) > 0:
                    id1 = 0
                    while not reads[id1][1] == j:
                        id1 += 1

                    startRead = reads[id1]

                    starts[reads[id1][0]] -= 1
                    ends[j] -= 1
                    del reads[id1]

                    foundRead = False
                    for r in range(len(pairedLensSorted)):
                        l = pairedLensSorted[r]
                        if j-l >= 0 and starts[j-l] > 0:
                            pairedLens[l] -= 1
                            if pairedLens[l] == 0:
                                del pairedLensSorted[r]
                            foundRead = True

                            # Add paired read to list
                            id2 = 0
                            while not reads[id2][0] == j-l:
                                id2 += 1

                            starts[j-l] -= 1
                            ends[reads[id2][1]] -= 1

                            paired += [[reads[id2], startRead]]
                            del reads[id2]

                            break
                    
                    if not foundRead:
                        readLen = startRead[1]-startRead[0]
                        if readLen in unpairedLens:
                            unpaired += [startRead]

                            if unpairedLens[readLen] == 1:
                                del unpairedLens[readLen]
                            else:
                                unpairedLens[readLen] -= 1
                        else:
                            unmatched += [startRead]

                    while i < j and starts[i] <= 0:
                        i += 1
                    while j > i and ends[j] <= 0:
                        j -= 1

            # Pair up remaining reads until we meet the quota of paired-end reads
            if len(pairedLensSorted) > 0:
                numPairedReads = 0
                for l in pairedLensSorted:
                    numPairedReads += pairedLens[l]

                for _ in range(numPairedReads):
                    paired += [[unmatched[0], unmatched[-1]]]
                    del unmatched[-1]
                    del unmatched[0]

            # Add remaining unmatched reads as unpaired reads
            for r in unmatched:
                unpaired += [r]

        # Add remaining reads as unpaired reads
        for r in reads:
            unpaired += [r]

        return unpaired, paired

    def findReads(self, readLens, lensLeft, lensRight, coverage):
        ''' Find the set of reads that most closely matches the distribution of readLens and the coverage vector
        '''

        countReads = 0
        for length,count in readLens.items():
            countReads += count

        fragmentLens = dict()
        countPaired = 0
        for length,count in lensLeft.items():
            countPaired += count
            fragmentLens[length] = count
        for length,count in lensRight.items():
            if length in fragmentLens:
                fragmentLens[length] += count
            else:
                fragmentLens[length] = count

        if countPaired == 0:
            pairedLens = dict()
            unpairedLens = readLens
            fragmentLens = readLens
        elif countPaired == countReads:
            pairedLens = readLens
            unpairedLens = dict()
        else:
            # Assign longest read lengths to be paired-end reads
            pairedLens = dict()
            unpairedLens = dict()
            sortedLens = sorted(readLens, reverse=True)
            i = 0
            while countPaired > 0:
                l = sortedLens[i]
                if countPaired >= readLens[l]:
                    pairedLens[l] = readLens[l]
                else:
                    pairedLens[l] = countPaired
                    unpairedLens[l] = readLens[l] - countPaired
                    if l in fragmentLens:
                        fragmentLens[l] += readLens[l] - countPaired
                    else:
                        fragmentLens[l] = readLens[l] - countPaired
                i += 1
            while i < len(sortedLens):
                l = sortedLens[i]
                unpairedLens[l] = readLens[l]
                if l in fragmentLens:
                    fragmentLens[l] += readLens[l] - countPaired
                else:
                    fragmentLens[l] = readLens[l] - countPaired

        reads = self.findReadsInCoverage_v1(coverage, fragmentLens)

        s = 'Reads (' + str(len(reads)) + '): ' + str(reads)
        s2 = 'Read lengths: ' + str(readLens)
        unpaired, paired = self.findPairs(len(coverage), reads, pairedLens, unpairedLens)

        if len(paired) > 0:
            print(s)
            print(s2)
            print('Unpaired ('+str(len(unpaired))+'): ' + str(unpaired))
            print('Paired: ('+str(len(paired))+')' + str(paired))
            exit()


        return unpaired, paired
        

    def findBestLeftLen(self, coverage, start, end, lensLeftSorted):
        # Find best left length for paired or unpaired reads
        bestScore = 0
        leftScores = [0]*len(lensLeftSorted)
        scoreWeight = 1
        for i in range(lensLeftSorted[0]):
            # penalize scores for fragments containing zeros
            if start+i > end or coverage[start+i] == 0:
                if coverages[start+i-1] > 0:
                    scoreWeight *= 0.5
            else:
                j = 0
                while j < len(lensLeftSorted) and lensLeftSorted[j] > i:
                    j += 1
                if j < len(lensLeftSorted) and lensLeftSorted[j] == i:
                    if start+i == end or coverage[start+i-1] > coverage[start+i]:
                        leftScores[j] = scoreWeight
                    else:
                        leftScores[j] = 0.8 * scoreWeight
                    if leftScores[j] > bestScore:
                        bestScore = leftScores[j]
        for i in range(len(leftScores)):
            if leftScores[i] == bestScore:
                return i, bestScore

    def findBestRightLen(self, coverage, start, end, lensRightSorted):
        bestScore = 0
        rightScores = [0]*len(lensRightSorted)
        scoreWeight = 1
        for i in range(min(end-start, lensRightSorted[0])):
            if coverage[start+l-i] == 0:
                if coverages[start+l-i+1] > 0:
                    scoreWeight *= 0.5
            else:
                j = 0
                while j < len(lensRightSorted) and lensRightSorted[j] > i:
                    j += 1
                if j < len(lensRightSorted) and lensRightSorted[j] == i:
                    if coverage[start+l-i] > coverage[start+l-j-1]:
                        rightScores[j] = scoreWeight
                    else:
                        rightScores[j] = 0.8 * scoreWeight
                    if rightScores[j] > bestScore:
                        bestScore = rightScores[j]
        for i in range(len(rightScores)):
            if rightScores[i] == bestScore:
                return i, bestScore

    def findBestFragmentLen(self, coverage, start, end, fragmentLensSorted, lensRightSorted):
        bestScore = 0
        bestFragmentLen = 0
        bestLenRight = 0
        for i in range(len(fragmentLensSorted)):
            l = fragmentLensSorted[i]

            if start+l > end:
                fragmentScore = 0
            elif start+l == end or coverage[start+l-1] > coverage[start+l]:
                fragmentScore = 1
            else:
                fragmentScore = 0.9

            if score > 0:
                lenRightId, rightScore = self.findBestRightLen(coverage, start, end, lensRightSorted)

                score = fragmentScore * rightScore
                if score == 1:
                    return i, lenRightId, score
                elif score > bestScore:
                    bestScore = score
                    bestFragmentLen = i
                    bestLenRight = lenRightId
        return bestFragmentLen, bestLenRight, bestScore

    def scoreFragment(self, coverage, start, end, overlapLen=0):
        ''' Return the score for the continuous fragment between start and end, weighted by the given value
            coverage should be run-length-encoded

            OverlapLen is the length of overlap between paired ends at the beginning of this read
        '''
        # reward if ends are a step above outside region
        startReward = False
        endReward = False

        numZeros = 0

        i = 0
        offset = 0
        while offset+coverage[i][1] <= start:
            offset += coverage[i][1]
            i += 1
        if offset == start and (i == 0 or coverage[i][0] > coverage[i-1][0]):
            startReward = True
        while offset < end:
            if coverage[i][0] <= 0:
                numZeros += min(end, offset+coverage[i][1]) - max(start, offset)
            offset += coverage[i][1]
            i += 1
        if offset == end and (i == len(coverage) or coverage[i][0] < coverage[i-1][0]):
            endReward = True

        score = 1
        if not startReward:
            score *= 0.95
        if not endReward:
            score *= 0.95
        if numZeros > 0:
            length = end-start
            score *= float(length - numZeros) / (2*length)

        return score

    def findReadsInCoverage_v1_new(self, coverage, fragmentLens, lensLeft, lensRight):
        ''' Given a coverage vector, return a set of reads the closely fits the coverage vector as well the distribution of read lengths.
            The algorithm creates new reads greedily, starting from both ends at once to ensure the ends of the vector are well marked.
            This algorithm is guaranteed to return a set of reads that covers every base at least to the corresponding depth of the coverage vector.
            In many cases the algorithm will overcompensate by creating extra reads to make sure every base in the coverage vector is covered.
            In such cases new reads have length equal to the median read length in the input distribution.

            coverage: Coverage vector containing the accumulation of many reads
            fragmentLens: Dictionary containing the distribution of all read lengths in the coverage vector
            boundaries: Indicies in the coverage vector which reads cannot cross (e.g. exon boundaries in the unspliced coverage vector)
        '''
        countFragments = 0
        for _,count in fragmentLens.items():
            countFragments += count
        totalFragments = countFragments

        countPaired = 0
        for _,count in lensLeft.items():
            countPaired += count

        countUnpaired = 0

        scoreThreshold = 0.01

        #lens = fragmentLens.keys()

        # Find max and mode read lengths
        #maxLen = max(lens)

        # Read lengths sorted by frequency, largest to smallest
        #fragmentLensSorted = sorted(fragmentLens, key=fragmentLens.get, reverse=True)
        #lensLeftSorted = sorted(lensLeft, key=lensLeft.get, reverse=True)
        #lensRightSorted = sorted(lensRight, key=lensRight.get, reverse=True)

        # Read lengths sorted largest to smallest
        fragmentLensSorted = sorted(fragmentLens, reverse=True)
        lensLeftSorted = sorted(lensLeft, reverse=True)
        lensRightSorted = sorted(lensRight, reverse=True)

        # create separate list with unpaired lengths
        unpairedLens = dict()
        unpairedLensSorted = []
        while countUnpaired < countFragments-countPaired:
            diff = countFragments - countPaired - countUnpaired
            smallestLen = fragmentLensSorted[-1]

            unpairedLensSorted = [smallestLen] + unpairedLensSorted
            if diff < fragmentLens[smallestLen]:
                unpairedLens[smallestLen] = diff
                countUnpaired += diff
                fragmentLens[smallestLen] -= diff
            else:
                unpairedLens[smallestLen] = fragmentLens[smallestLen]
                countUnpaired += fragmentLens[smallestLen]
                del fragmentLens[smallestLen]
                del fragmentLensSorted[-1]

        unpaired = []
        paired = []

        # start and end of coverage window
        # Keep finding reads from both ends until they meet in the middle
        start = 0
        startOffset = 0
        end = len(coverage)

        covLength = 0
        for segment in coverage:
            covLength += segment[1]
        endOffset = covLength

        while start < end and coverage[start][0] <= 0:
            startOffset += coverage[start][1]
            start += 1
        while end > start and coverage[end-1][0] <= 0:
            endOffset -= coverage[end-1][1]
            end -= 1

        unpairedWeight = 0.5

        while end > start and (len(fragmentLensSorted) + len(unpairedLensSorted) > 0):
            # find a read from the beginning
            # find best unpaired read
            bestUnpairedId = 0
            bestUnpairedScore = 0

            # penalty for unpaired reads to encourage paired-end reads
            for i in range(len(unpairedLensSorted)):
                score = unpairedWeight * self.scoreFragment(coverage, startOffset, startOffset+unpairedLensSorted[i])**2
                if score > bestUnpairedScore:
                    bestUnpairedScore = score
                    bestUnpairedId = i

            
            # find best left length for paired reads
            bestLeftId = 0
            bestLeftScore = 0
            for i in range(len(lensLeftSorted)):
                score = self.scoreFragment(coverage, startOffset, startOffset+lensLeftSorted[i])
                if score > bestLeftScore:
                    bestLeftScore = score
                    bestLeftId = i

            if bestLeftScore > 0:
                # find best fragment length and right length combination for paired reads
                bestFragmentId = 0
                bestRightId = 0
                bestFragmentScore = 0
                for i in range(len(fragmentLensSorted)):
                    for j in range(len(lensRightSorted)):
                        score = self.scoreFragment(coverage, startOffset+fragmentLensSorted[i]-lensRightSorted[j], startOffset+fragmentLensSorted[i])
                        if score > bestFragmentScore:
                            bestFragmentScore = score
                            bestFragmentId = i
                            bestRightId = j

                bestPairedScore = bestLeftScore * bestFragmentScore
            else:
                bestPairedScore = 0

            if bestPairedScore > scoreThreshold and bestPairedScore >= bestUnpairedScore:
                # Add paired read
                lenLeft = lensLeftSorted[bestLeftId]
                fragmentLen = fragmentLensSorted[bestFragmentId]
                lenRight = lensRightSorted[bestRightId]
                paired.append([[startOffset, startOffset+lenLeft], [startOffset+fragmentLen-lenRight, startOffset+fragmentLen]])

                # Update coverage vector
                if lenLeft+lenRight > fragmentLen:
                    coverage = self.updateRLE(coverage, startOffset, fragmentLen, -1)
                else:
                    coverage = self.updateRLE(coverage, startOffset, lenLeft, -1)
                    coverage = self.updateRLE(coverage, startOffset+fragmentLen-lenRight, lenRight, -1)

                # remove from length distributions
                if fragmentLens[fragmentLen] == 1:
                    del fragmentLens[fragmentLen]
                    del fragmentLensSorted[bestFragmentId]
                else:
                    fragmentLens[fragmentLen] -= 1

                if lensRight[lenRight] == 1:
                    del lensRight[lenRight]
                    del lensRightSorted[bestRightId]
                else:
                    lensRight[lenRight] -= 1

                if lensLeft[lenLeft] == 1:
                    del lensLeft[lenLeft]
                    del lensLeftSorted[bestLeftId]
                else:
                    lensLeft[lenLeft] -= 1
            elif bestUnpairedScore > scoreThreshold:
                # Add unpaired read
                length = unpairedLensSorted[bestUnpairedId]
                unpaired.append([startOffset, startOffset+length])

                # Update coverage vector
                coverage = self.updateRLE(coverage, startOffset, min(length, len(coverage)-startOffset), -1)

                # remove from length distribution
                if unpairedLens[length] == 1:
                    del unpairedLens[length]
                    del unpairedLensSorted[bestUnpairedId]
                else:
                    unpairedLens[length] -= 1
            else:
                if bestLeftScore > scoreThreshold:
                    # Add right as unpaired
                    length = lensLeftSorted[bestLeftId]
                    unpaired.append([startOffset, startOffset+length])

                    # Update coverage vector
                    coverage = self.updateRLE(coverage, startOffset, min(length, len(coverage)-startOffset), -1)

                else:
                    # No decent matches
                    coverage[start][0] = 0

            # update start and end
            start = 0
            startOffset = 0
            end = len(coverage)
            endOffset = covLength

            while start < end and coverage[start][0] <= 0:
                startOffset += coverage[start][1]
                start += 1
            while end > start and coverage[end-1][0] <= 0:
                endOffset -= coverage[end-1][1]
                end -= 1

            if end > start and (len(fragmentLensSorted) + len(unpairedLensSorted) > 0):
                # find a read from the end
                # find best unpaired read
                bestUnpairedId = 0
                bestUnpairedScore = 0

                # penalty for unpaired reads to encourage paired-end reads
                for i in range(len(unpairedLensSorted)):
                    score = unpairedWeight * self.scoreFragment(coverage, endOffset-unpairedLensSorted[i], endOffset)**2
                    if score > bestUnpairedScore:
                        bestUnpairedScore = score
                        bestUnpairedId = i

                
                # find best left length for paired reads
                bestRightId = 0
                bestRightScore = 0
                for i in range(len(lensRightSorted)):
                    score = self.scoreFragment(coverage, endOffset-lensRightSorted[i], endOffset)
                    if score > bestRightScore:
                        bestRightScore = score
                        bestRightId = i

                if bestRightScore > 0:
                    # find best fragment length and right length combination for paired reads
                    bestFragmentId = 0
                    bestLeftId = 0
                    bestFragmentScore = 0
                    for i in range(len(fragmentLensSorted)):
                        for j in range(len(lensLeftSorted)):
                            score = self.scoreFragment(coverage, endOffset-fragmentLensSorted[i], endOffset-fragmentLensSorted[i]+lensLeftSorted[j])
                            if score > bestFragmentScore:
                                bestFragmentScore = score
                                bestFragmentId = i
                                bestLeftId = j

                    bestPairedScore = bestRightScore * bestFragmentScore
                else:
                    bestPairedScore = 0

                if bestPairedScore > scoreThreshold and bestPairedScore >= bestUnpairedScore:
                    # Add paired read
                    lenLeft = lensLeftSorted[bestLeftId]
                    fragmentLen = fragmentLensSorted[bestFragmentId]
                    lenRight = lensRightSorted[bestRightId]
                    paired.append([[endOffset-fragmentLen, endOffset-fragmentLen+lenLeft], [endOffset-lenRight, endOffset]])

                    # Update coverage vector
                    if lenLeft+lenRight > fragmentLen:
                        coverage = self.updateRLE(coverage, endOffset-fragmentLen, fragmentLen, -1)
                    else:
                        coverage = self.updateRLE(coverage, endOffset-fragmentLen, lenLeft, -1)
                        coverage = self.updateRLE(coverage, endOffset-lenRight, lenRight, -1)

                    # remove from length distributions
                    if fragmentLens[fragmentLen] == 1:
                        del fragmentLens[fragmentLen]
                        del fragmentLensSorted[bestFragmentId]
                    else:
                        fragmentLens[fragmentLen] -= 1

                    if lensRight[lenRight] == 1:
                        del lensRight[lenRight]
                        del lensRightSorted[bestRightId]
                    else:
                        lensRight[lenRight] -= 1

                    if lensLeft[lenLeft] == 1:
                        del lensLeft[lenLeft]
                        del lensLeftSorted[bestLeftId]
                    else:
                        lensLeft[lenLeft] -= 1
                elif bestUnpairedScore > scoreThreshold:
                    # Add unpaired read
                    length = unpairedLensSorted[bestUnpairedId]
                    unpaired.append([endOffset-length, endOffset])

                    # Update coverage vector
                    coverage = self.updateRLE(coverage, max(0,endOffset-length), min(length, endOffset), -1)

                    # remove from length distribution
                    if unpairedLens[length] == 1:
                        del unpairedLens[length]
                        del unpairedLensSorted[bestUnpairedId]
                    else:
                        unpairedLens[length] -= 1
                else:
                    if bestRightScore > scoreThreshold:
                        # Add right as unpaired
                        length = lensRightSorted[bestRightId]
                        unpaired.append([endOffset-length, endOffset])

                        # Update coverage vector
                        coverage = self.updateRLE(coverage, max(0,endOffset-length), min(length, endOffset), -1)

                    else:
                        # No decent matches
                        coverage[end-1][0] = 0

                # update start and end
                start = 0
                startOffset = 0
                end = len(coverage)
                endOffset = covLength

                while start < end and coverage[start][0] <= 0:
                    startOffset += coverage[start][1]
                    start += 1
                while end > start and coverage[end-1][0] <= 0:
                    endOffset -= coverage[end-1][1]
                    end -= 1
        if (len(fragmentLensSorted) + len(unpairedLensSorted) > 0) or (len(unpaired)+len(paired) < totalFragments):
            print('Found %d of %d, %d left' % (len(unpaired) + len(paired), totalFragments, len(fragmentLensSorted) + len(unpairedLensSorted)))
            exit()
        return unpaired, paired

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
            print('%d boundaries' % len(boundaries))
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
            for length in range(1, min(maxLen+1, boundaries[boundBottom]-start)):
                if (readStart+length == end) or (readStart+length < end and coverage[readStart + length] < coverage[readStart + length - 1]):
                    if length in readLens:
                        readEnd = readStart + length
                        reads.append([readStart, readEnd])

                        readLens[length] -= 1

                        # reorder sorted lengths
                        for i in range(len(lensSorted)):
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
                        for i in range(len(lensSorted)):
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
            for i in range(readStart, readEnd):
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
                for length in range(1, min(maxLen+1, end-boundaries[boundTop])):
                #for length in range(1, maxLen+1):
                    if (end-length == start) or (end-length > start and coverage[end - length] > coverage[end - length - 1]):
                        if length in readLens:
                            readStart = readEnd - length
                            reads.append([readStart, readEnd])

                            readLens[length] -= 1
                                

                            # reorder sorted lengths
                            for i in range(len(lensSorted)):
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
                            for i in range(len(lensSorted)):
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

                for i in range(readStart, readEnd):
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
        lensSorted = [v[0] for v in sorted(readLens.iteritems(), key = lambda k,v: (v, -k), reverse=True)]


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
                    for i in range(len(lensSorted)):
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
                for i in range(len(lensSorted)):
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
            for i in range(readStart, readEnd):
                coverage[i] -= 1

            # update start
            while start < end and coverage[start] <= 0:
                start += 1
            while end > start and coverage[end-1] <= 0:
                end -= 1


            if end > start and len(lensSorted) > 0:
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
                        for i in range(len(lensSorted)):
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
                    for i in range(len(lensSorted)):
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

                for i in range(readStart, readEnd):
                    coverage[i] -= 1

                # update end
                while coverage[end-1] <= 0 and end > start:
                    end -= 1
                while coverage[start] <= 0 and start < end:
                    start += 1

        if not len(reads) == countReads:
            print('Error! %d =/= %d!' % (countReads, len(reads)))
            exit()
        return reads

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

        for c in self.chromosomeNames:
            if self.chromosomes[c] > index:
                return c
            else:
                index -= self.chromosomes[c]

    def getChromosomeAndIndex(self, index):
        ''' Return chromosome name containing the given index from the whole-genome vector
        '''

        for c in self.chromosomeNames:
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
            print('Error! Read exons = None')
            exit()
        if len(read.exons) == 0:
            print('Error! Exon length = 0')

        # Update read index based on chromosome
        offset = self.chromOffsets[read.chrom]
        for i in range(len(read.exons)):
            read.exons[i] = [read.exons[i][0]+offset, read.exons[i][1]+offset]

        # update list of exons
        alignment = read.exons
        if len(alignment) > 1:
            for i in range(len(alignment)-1):
                self.exons.add(offset + alignment[i][1])
                self.exons.add(offset + alignment[i+1][0])

        
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

            key = (pair_index, pair_chrom, read.exons[0][0], read.chrom)

            if key in self.unmatched:
                foundMatch = False
                for i in range(len(self.unmatched[key])):
                    match = self.unmatched[key][i]

                    if read.NH == match.NH:
                        if len(self.unmatched[key]) == 1:
                            del self.unmatched[key]
                        else:
                            del self.unmatched[key][i]

                        xs = read.xs or match.xs
                        NH = read.NH

                        if self.conflicts(read.exons, match.exons):
                            if len(read.exons) == 1:
                                self.addUnspliced(read)
                            else:
                                self.addSpliced(read)

                            if len(match.exons) == 1:
                                self.addUnspliced(match)
                            else:
                                self.addSpliced(match)
                        else:
                            self.paired.append(pairedread.PairedRead(pair_chrom, match.exons, read.chrom, read.exons, xs, NH))


                        foundMatch = True
                        break
                if not foundMatch:
                    if (read.exons[0][0], read.chrom, pair_index, pair_chrom) in self.unmatched:
                        self.unmatched[(read.exons[0][0], read.chrom, pair_index, pair_chrom)] += [read]
                    else:
                        self.unmatched[(read.exons[0][0], read.chrom, pair_index, pair_chrom)] = [read]
            else:
                if (read.exons[0][0], read.chrom, pair_index, pair_chrom) in self.unmatched:
                    self.unmatched[(read.exons[0][0], read.chrom, pair_index, pair_chrom)] += [read]
                else:
                    self.unmatched[(read.exons[0][0], read.chrom, pair_index, pair_chrom)] = [read]

    def conflicts(self, exonsA, exonsB):
        '''
            Returns true if any of the exons from A or B overlaps one of the introns from the other set of exons
        '''
        for e in exonsB:
            if e[0] > exonsA[-1][0]:
                break

            for i in range(len(exonsA)-1):
                if e[0] >= exonsA[-i-1][0]:
                    break
                elif e[1] > exonsA[-i-2][1]:
                    return True

        countA = len(exonsA)
        for i in range(countA):
            e = exonsA[countA-i-1]
            if e[1] < exonsB[0][0]:
                break

            for i in range(len(exonsB)-1):
                if e[1] <= exonsB[i][1]:
                    break
                elif e[1] > exonsB[i][1] and e[0] < exonsB[i+1][0]:
                    return True

        return False

    def writeSAM(self, filehandle):
        ''' Write all alignments to a SAM file
        '''
        
        # write header
        filehandle.write('@HD\tVN:1.0\tSO:unsorted\n')
        for k,v in self.chromosomes.items():
            filehandle.write('@SQ\tSN:' + str(k) + '\tLN:' + str(v) + '\n')

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
            for i in range(1, len(exons)):
                if exons[i][0] - exons[i-1][1] == 0:
                    prevLen = int(cigar[-1][:-1])
                    cigar[-1] = str(prevLen + exons[i][1] - exons[i][0]) + 'M'
                else:
                    cigar += [str(exons[i][0] - exons[i-1][1]) + 'N']
                    cigar += [str(exons[i][1] - exons[i][0]) + 'M']
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
            for i in range(1, len(exonsA)):
                if exonsA[i][0] - exonsA[i-1][1] == 0:
                    prevLen = int(cigarA[-1][:-1])
                    cigarA[-1] = str(prevLen + exonsA[i][1] - exonsA[i][0]) + 'M'
                else:
                    if exonsA[i][1] == exonsA[i][0]:
                        print(exonsA)
                        exit()

                    cigarA += [str(exonsA[i][0] - exonsA[i-1][1]) + 'N']
                    cigarA += [str(exonsA[i][1] - exonsA[i][0]) + 'M']
            cigarA = ''.join(cigarA)

            exonsB = pair.exonsB
            cigarB = [str(exonsB[0][1] - exonsB[0][0]) + 'M']
            for i in range(1, len(exonsB)):
                if exonsB[i][0] - exonsB[i-1][1] == 0:
                    prevLen = int(cigarB[-1][:-1])
                    cigarB[-1] = str(prevLen + exonsB[i][1] - exonsB[i][0]) + 'M'
                else:
                    ####
                    if exonsB[i][1] == exonsB[i][0]:
                        print(exonsB)
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