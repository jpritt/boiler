#! /usr/bin/env python

import bisect
import junction
import read
import pairedread
import os.path
import sys
import copy
import random

from random import shuffle

import time

class Alignments:
    ''' A set of reads aligned to a genome '''

    def __init__(self, chromosomes, debug=False):
        ''' Initialize a genome for alignments

            chromosomes: A dictionary with keys corresponding to chromosome names and values corresponding to lengths
        '''

        self.debug = debug

        self.origPairedDepth = 0
        self.calcPairedDepth = 0
        self.countDense = 0

        self.pairedCounts = [0] * 14
        self.unpairedCounts = [0] * 14


        self.totalPaired = 0
        self.totalUnpaired = 0

        self.countUnpaired = 0
        self.badUnpaired = 0
        self.countPaired = 0
        self.badPaired = 0

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

        # List of potential gene boundaries as tuples
        self.gene_bounds = []
        # If 2 reads are less than this far apart, combine them into a single gene
        self.overlap_radius = 50

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

    def finalizeUnmatched(self):
        # Finalize unmatched (discordant) reads
        count = 0

        for name,reads in self.unmatched.items():
            if reads:
                count += len(reads)
                for r in reads:
                    self.updateGeneBounds([r.exons[0][0], r.exons[-1][1]])

                    if len(r.exons) == 1:
                        self.addUnspliced(r)
                    else:
                        self.addSpliced(r)

        if self.debug:
            print('%d unmatched' % count)

    def finalizeUnmatched_profile(self):
        import cProfile
        import io
        import pstats
        pr = cProfile.Profile()
        pr.enable()

        start = time.time()
        self.finalizeUnmatched()
        end = time.time()
        print('Total time: %0.3f s' % (end-start))
        print('')

        pr.disable()
        s = io.StringIO()
        sortby = 'tottime'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats(30)
        print(s.getvalue())

        exit()

    def finalizeExons(self):
        ''' Convert the set of exon boundaries to a list
        '''


        if self.debug:
            print('Sorting spliced subexons')
        splice_exons = sorted(list(self.exons))

        # Convert gene boundaries to exon boundaries
        if self.debug:
            print('Sorting gene subexons')
        gene_exons = [0] * (2*len(self.gene_bounds))
        for i in range(len(self.gene_bounds)):
            gene_exons[2*i] = self.gene_bounds[i][0]
            gene_exons[2*i+1] = self.gene_bounds[i][1]
        if self.debug:
            print('%d gene subexons' % len(gene_exons))
            print('%d spliced subexons' % len(self.exons))

        gene_exons.sort()

        # Merge the 2 exon lists
        self.exons = [0] * (len(splice_exons) + len(gene_exons))

        splice_i = 0
        gene_i = 0
        combined_i = 0
        while splice_i < len(splice_exons) and gene_i < len(gene_exons):
            if splice_exons[splice_i] < gene_exons[gene_i]:
                self.exons[combined_i] = splice_exons[splice_i]
                splice_i += 1
            elif splice_exons[splice_i] > gene_exons[gene_i]:
                self.exons[combined_i] = gene_exons[gene_i]
                gene_i += 1
            else:
                self.exons[combined_i] = splice_exons[splice_i]
                splice_i += 1
                gene_i += 1
            combined_i += 1

        if splice_i == len(splice_exons):
            while gene_i < len(gene_exons):
                self.exons[combined_i] = gene_exons[gene_i]
                gene_i += 1
            combined_i += 1
        elif gene_i == len(gene_exons):
            while splice_i < len(splice_exons):
                self.exons[combined_i] = splice_exons[splice_i]
                splice_i += 1
            combined_i += 1

        if combined_i < len(self.exons):
            self.exons = self.exons[:combined_i]

    def finalizeReads(self):
        ''' Now that all exon boundaries are known, fix unspliced regions that cross exon boundaries
            and finalized paired-end reads
        '''

        # TODO: Find splice sites in spliced reads, paired reads, and discordant

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
                self.spliced += [r]

                del self.unspliced[x]
            else:
                x += 1

        for r in self.unspliced:
            r.readLen = r.exons[0][1] - r.exons[0][0]

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
            r.endOffset = self.exons[r.exonIds[-1] + 1] - exons[-1][1]

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

            # TODO: Replace gap with end lengths
            if n > 0:
                newRead = read.Read(pair.chromA, exonsA[:-1] + [[exonsA[-1][0], exonsB[n-1][1]]] + exonsB[n:], pair.xs, pair.NH)
                gap = 0
                for i in range(n):
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

    def findPairedDist(self, dists, assigned, numReads, value):
        for i in range(numReads):
            if assigned[i]:
                continue
            for j in range(i):
                if assigned[j]:
                    continue
                if dists[i][j] == value:
                    return i,j

    def findSingleDist(self, dists, assigned, numReads, value):
        for i in range(numReads):
            if assigned[i]:
                continue
            if dists[i][i] == value:
                return i

    def findReadsWithPairedCov(self, unpaired_lens, paired_lens, lens_left, lens_right, coverage_fragments):
        fragmentLens = dict()
        for k,v in unpaired_lens.items():
            fragmentLens[k] = v
        for k,v in paired_lens.items():
            if k in fragmentLens:
                fragmentLens[k] += v
            else:
                fragmentLens[k] = v

        countPaired = 0
        for k,v in lens_left:
            countPaired += v

        reads = self.findReadsInCoverage_v1(coverage_fragments, fragmentLens)
        reads.sort(key=lambda x: x[1]-x[0], reverse=True)

        lens_left_sorted = sorted(lens_left, key=lens_left.get, reverse=True)
        lens_right_sorted = sorted(lens_left, key=lens_right.get, reverse=True)

        unpaired = []
        paired = []
        for n in range(len(reads)):
            if countPaired == 0:
                break

            r = reads[n]
            l = r[1]-r[0]
            if l in unpaired_lens and unpaired_lens[l] > 0:
                unpaired.append(r)
                unpaired_lens[l] -= 1
            else:
                countPaired -= 1
                if lens_left_sorted[-1] > l or lens_right_sorted[-1] > l:
                    paired.append([r, r])
                else:
                    i = 0
                    while lens_left_sorted[i] > l:
                        i += 1
                    ll = lens_left_sorted[i]
                    lens_left[ll] -= 1
                    if lens_left[ll] <= 0:
                        del lens_left_sorted[i]

                    i = 0
                    while lens_right_sorted[i] > l:
                        i += 1
                    lr = lens_right_sorted[i]
                    lens_right[lr] -= 1
                    if lens_right[lr] <= 0:
                        del lens_right_sorted[i]

                    paired.append([[r[0], r[0]+ll], [r[1]-lr, r[1]]])
        unpaired += reads[n:]

        return unpaired, paired


    def findPairsGreedy(self, reads, pairedLens):
        #for k,v in pairedLens.items():
        #    self.origPairedDepth += k * v

        pairedLensSorted = sorted(pairedLens.keys(), reverse=True)

        reads.sort()
        numReads = len(reads)

        # Keep track of which reads we have already assigned
        assigned = [0] * numReads

        # Map each distance to the paired reads that match it
        pairDists = dict()
        singleDists = dict()

        # Create a distance matrix between all pairs of reads
        dists = [0] * numReads
        for i in range(numReads):
            dists[i] = [0] * (i+1)
            for j in range(i):
                d = reads[i][1] - reads[j][0]
                dists[i][j] = d
                if d in pairDists:
                    pairDists[d] += 1
                else:
                    pairDists[d] = 1

            d = reads[i][1] - reads[i][0]
            dists[i][i] = d
            if d in singleDists:
                singleDists[d] += 1
            else:
                singleDists[d] = 1

        paired = []
        countPaired = 0
        for k,v in pairedLens.items():
            countPaired += v
        while countPaired > 0:
            bestFreq = None
            bestL = None
            foundPair = False

            # Look for unique paired read lengths
            for l in pairedLensSorted:
                if not (l in pairDists and pairDists[l] > 0 and pairedLens[l] > 0):
                    continue

                expected = pairedLens[l]
                freq = pairDists[l]
                if freq == 0:
                    continue
                if freq <= expected:
                    pairedLens[l] -= 1
                    i,j = self.findPairedDist(dists, assigned, numReads, l)
                    paired.append([reads[j], reads[i]])

                    assigned[i] = 1
                    assigned[j] = 1

                    if dists[i][i] > 0:
                        singleDists[dists[i][i]] -= 1
                        dists[i][i] = 0
                    if dists[j][j] > 0:
                        singleDists[dists[j][j]] -= 1
                        dists[j][j] = 0

                    for x in range(i):
                        if dists[i][x] > 0:
                            pairDists[dists[i][x]] -= 1
                            dists[i][x] = 0
                    for x in range(i+1,numReads):
                        if dists[x][i] > 0:
                            pairDists[dists[x][i]] -= 1
                            dists[x][i] = 0
                    for x in range(j):
                        if dists[j][x] > 0:
                            pairDists[dists[j][x]] -= 1
                            dists[j][x] = 0
                    for x in range(j+1,numReads):
                        if dists[x][j] > 0:
                            pairDists[dists[x][j]] -= 1
                            dists[x][j] = 0

                    foundPair = True
                    countPaired -= 1
                    break
                elif bestFreq == None or (freq-expected) < bestFreq:
                    bestFreq = freq - expected
                    bestL = l

            # No unique paired lengths, so choose one from the least frequent
            if not foundPair:
                if bestFreq == None:
                    break
                else:
                    pairedLens[bestL] -= 1
                    i,j = self.findPairedDist(dists, assigned, numReads, bestL)
                    paired.append([reads[j], reads[i]])

                    assigned[i] = 1
                    assigned[j] = 1

                    if dists[i][i] > 0:
                        singleDists[dists[i][i]] -= 1
                        dists[i][i] = 0
                    if dists[j][j] > 0:
                        singleDists[dists[j][j]] -= 1
                        dists[j][j] = 0

                    for x in range(i):
                        if dists[i][x] > 0:
                            pairDists[dists[i][x]] -= 1
                            dists[i][x] = 0
                    for x in range(i+1,numReads):
                        if dists[x][i] > 0:
                            pairDists[dists[x][i]] -= 1
                            dists[x][i] = 0
                    for x in range(j):
                        if dists[j][x] > 0:
                            pairDists[dists[j][x]] -= 1
                            dists[j][x] = 0
                    for x in range(j+1,numReads):
                        if dists[x][j] > 0:
                            pairDists[dists[x][j]] -= 1
                            dists[x][j] = 0

                    countPaired -= 1

        remaining = [0] * (numReads - sum(assigned))
        i = 0
        for j in range(numReads):
            if not assigned[j]:
                remaining[i] = reads[j]
                i += 1

        remainingPairedLens = dict()
        for k,v in pairedLens.items():
            if v > 0:
                remainingPairedLens[k] = v

        if countPaired > 0:
            newUnpaired, newPaired = self.findPairsGreedy2(remaining, remainingPairedLens)
            unpaired = newUnpaired
            paired += newPaired
        else:
            unpaired = remaining

        #for p in paired:
        #    self.calcPairedDepth += p[1][1] - p[0][0]

        return unpaired, paired

    def findPairsGreedy2(self, reads, pairedLens):
        origPairedDepth = 0
        for k,v in pairedLens.items():
            origPairedDepth += k * v
            self.origPairedDepth += k * v

        numReads = len(reads)
        reads.sort()
        pairedLensSorted = sorted(pairedLens.keys(), reverse=True)

        paired = []

        countPaired = 0
        for k,v in pairedLens.items():
            countPaired += v
        countUnpaired = numReads - 2 * countPaired

        # Create a distance matrix between all pairs of reads
        dists = [0] * numReads
        for i in range(numReads):
            dists[i] = [0] * i
            for j in range(i):
                d = reads[i][1] - reads[j][0]
                dists[i][j] = d
        assigned = [0] * numReads

        lenIndex = 0
        while countPaired > 0:
            targetL = pairedLensSorted[lenIndex]

            bestL = None
            bestDiff = None
            bestPos = None

            for i in range(numReads):
                for j in range(i,0,-1):
                    l = dists[i][j-1]
                    diff = abs(l - targetL)
                    if l > 0 and (bestDiff == None or diff < bestDiff):
                        bestDiff = diff
                        bestL = dists[i][j-1]
                        bestPos = (j-1, i)
                    elif l > targetL:
                        break

            if bestL == None:
                break
            else:
                pairedLens[targetL] -= 1
                if pairedLens[targetL] == 0:
                    lenIndex += 1

                i = bestPos[0]
                j = bestPos[1]
                paired.append([reads[i], reads[j]])
                assigned[i] = 1
                assigned[j] = 1

                for x in range(i):
                    dists[i][x] = 0
                for x in range(i+1,numReads):
                    dists[x][i] = 0
                for x in range(j):
                    dists[j][x] = 0
                for x in range(j+1,numReads):
                    dists[x][j] = 0

                countPaired -= 1

        unpaired = [0] * countUnpaired
        i = 0
        for j in range(numReads):
            if not assigned[j]:
                unpaired[i] = reads[j]
                i += 1

        calcPairedDepth = 0
        for p in paired:
            calcPairedDepth += p[1][1] - p[0][0]
            self.calcPairedDepth += p[1][1] - p[0][0]

        return unpaired, paired

    def findPairsWithBoundaries(self, reads, unpaired_lens, paired_lens, boundaries):
        '''
        Use the fact that all reads must span all the subexons to improve our pairing

        :param reads:
        :param unpaired_lens:
        :param paired_lens:
        :param lens_left:
        :param lens_right:
        :param boundaries: Boundaries between subexons
        :return:
        '''

        numUnpaired = 0
        for k,v in unpaired_lens.items():
            numUnpaired += v
        numPaired = 0
        for k,v in paired_lens.items():
            numPaired += v
        numReads = len(reads)

        left_reads = []
        spanning_reads = []
        right_reads = []

        for r in reads:
            if r[0] < boundaries[0]:
                if r[1] > boundaries[-1]:
                    spanning_reads.append(r)
                else:
                    left_reads.append(r)
            elif r[1] > boundaries[-1]:
                right_reads.append(r)
            else:
                print('Read does not overlap left or right...?')
                exit()

        left_reads.sort()
        spanning_reads.sort()
        right_reads.sort()
        paired_lens_sorted = sorted(paired_lens, reverse=True)

        countUnpaired = 0
        for k,v in unpaired_lens.items():
            countUnpaired += v
        if len(spanning_reads) < countUnpaired and self.debug:
            print('%d < %d' % (len(spanning_reads), countUnpaired))

        countPaired = 0
        for k,v in paired_lens.items():
            countPaired += v

        paired = []
        unpaired = None
        while countPaired > 0:
            if not unpaired and len(spanning_reads) <= countUnpaired:
                unpaired = spanning_reads
                spanning_reads = []

            if (not left_reads and not right_reads) or (not left_reads and not spanning_reads) or (not right_reads and not spanning_reads):
                break

            if len(left_reads) >= len(right_reads):
                p = self.findLeftPair(left_reads, spanning_reads, right_reads, paired_lens_sorted)
                l = p[1][1] - p[0][0]
                if l in paired_lens_sorted:
                    paired_lens[l] -= 1
                    if paired_lens[l] <= 0:
                        i = 0
                        while not paired_lens_sorted[i] == l:
                            i += 1
                        del paired_lens_sorted[i]
                paired.append(p)
                countPaired -= 1
            else:
                p = self.findRightPair(left_reads, spanning_reads, right_reads, paired_lens_sorted)
                l = p[1][1] - p[0][0]
                if l in paired_lens_sorted:
                    paired_lens[l] -= 1
                    if paired_lens[l] <= 0:
                        i = 0
                        while not paired_lens_sorted[i] == l:
                            i += 1
                        del paired_lens_sorted[i]
                paired.append(p)
                countPaired -= 1

        s1 = 'After pairing: %d paired remaining, %d left, %d right' % (countPaired, len(left_reads), len(right_reads))
        remaining_reads = left_reads + spanning_reads + right_reads
        i = 0
        j = len(remaining_reads)
        for _ in range(countPaired):
            paired.append([remaining_reads[i], remaining_reads[j-1]])
            i += 1
            j -= 1
        s2 = str(i) + ', ' + str(j)
        if unpaired:
            unpaired += remaining_reads[i:j]
        else:
            unpaired = remaining_reads[i:j]
        unpaired.sort()

        if len(unpaired) != countUnpaired and self.debug:
            print(s1)
            print(s2)
            print('%d reads found, looking for %d unpaired and %d paired' % (numReads, numUnpaired, numPaired))
            print('Found %d unpaired reads (instead of %d) and %d paired' % (len(unpaired), countUnpaired, len(paired)))
            print('')

        return unpaired, paired


    def findLeftPair(self, left_reads, spanning_reads, right_reads, paired_lens_sorted):
        '''
        :param left_reads:
        :param spanning_reads:
        :param right_reads:
        :param paired_lens:
        :return:
        '''

        if not left_reads:
            return None

        r = left_reads[0]

        # Index of closest non-exact pair
        closest_i = None
        # Error distance for this pair
        closest_d = None
        # True if closest pair is in right_reads, false if in spanning_reads
        closest_right = True
        for l in paired_lens_sorted:
            end = r[0] + l

            for i in range(len(right_reads)):
                s = right_reads[i]
                if s[1] == end:
                    del left_reads[0]
                    del right_reads[i]
                    return [r, s]
                elif not closest_d or abs(s[1]-end) < closest_d:
                    closest_d = abs(s[1]-end)
                    closest_i = i

        for l in paired_lens_sorted:
            end = r[0] + l

            for i in range(len(spanning_reads)):
                s = spanning_reads[i]
                if s[1] == end:
                    del left_reads[0]
                    del spanning_reads[i]
                    return [r, s]
                elif not closest_d or abs(s[1]-end) < closest_d:
                    closest_d = abs(s[1]-end)
                    closest_i = i
                    closest_right = False

        if closest_i == None:
            print('Error! No closest pair found!')
            exit()
        else:
            del left_reads[0]
            if closest_right:
                s = right_reads[closest_i]
                del right_reads[closest_i]
            else:
                s = spanning_reads[closest_i]
                del spanning_reads[closest_i]
            return [r, s]

    def findRightPair(self, left_reads, spanning_reads, right_reads, paired_lens_sorted):
        '''
        :param left_reads:
        :param spanning_reads:
        :param right_reads:
        :param paired_lens:
        :return:
        '''

        if not right_reads:
            return None

        r = right_reads[-1]

        # Index of closest non-exact pair
        closest_i = None
        # Error distance for this pair
        closest_d = None
        # True if closest pair is in left_reads, false if in spanning_reads
        closest_left = True
        for l in paired_lens_sorted:
            start = r[1] - l

            for i in range(len(left_reads)):
                s = left_reads[i]
                if s[0] == start:
                    del right_reads[-1]
                    del left_reads[i]
                    return [s, r]
                elif not closest_d or abs(s[0]-start) < closest_d:
                    closest_d = abs(s[0]-start)
                    closest_i = i

        for l in paired_lens_sorted:
            start = r[1] - l

            for i in range(len(spanning_reads)):
                s = spanning_reads[i]
                if s[0] == start:
                    del right_reads[-1]
                    del spanning_reads[i]
                    return [s, r]
                elif not closest_d or abs(s[0]-start) < closest_d:
                    closest_d = abs(s[0]-start)
                    closest_i = i
                    closest_left = False

        if closest_i == None:
            print('Error! No closest pair found!')
            exit()
        else:
            del right_reads[-1]
            if closest_left:
                s = left_reads[closest_i]
                del left_reads[closest_i]
            else:
                s = spanning_reads[closest_i]
                del spanning_reads[closest_i]
            return [s, r]





    def findPairsRandom(self, reads, paired_lens):
        countPairs = 0
        for k,v in paired_lens.items():
            countPairs += v

        reads.sort()

        unmatched = []
        paired = []
        while countPairs > 0 and reads:
            p = self.findLeftPairRandom(reads, paired_lens)
            if p:
                paired.append(p)
                countPairs -= 1
            else:
                unmatched.append(reads[0])
                del reads[0]

            if countPairs == 0 or not reads:
                break

            p = self.findRightPairRandom(reads, paired_lens)
            if p:
                paired.append(p)
                countPairs -= 1
            else:
                unmatched.append(reads[-1])
                del reads[-1]

        reads += unmatched
        reads.sort()

        while countPairs > 0:
            p = self.findClosestLeftPair(reads, paired_lens)
            paired.append(p)
            countPairs -= 1

            if countPairs == 0:
                break

            p = self.findClosestRightPair(reads, paired_lens)
            paired.append(p)
            countPairs -= 1

        return reads, paired


    def findLeftPairRandom(self, reads, paired_lens):
        start = reads[0][0]

        matches = []

        for i in range(1, len(reads)):
            end = reads[i][1]
            if end-start in paired_lens:
                matches.append(i)

        if matches:
            # If there is at least one exact match, return one of them at random
            i = random.choice(matches)

            pair = [reads[0], reads[i]]

            l = reads[i][1] - start
            if paired_lens[l] > 1:
                paired_lens[l] -= 1
            else:
                del paired_lens[l]
            del reads[i]
            del reads[0]

            return pair
        else:
            return None

    def findRightPairRandom(self, reads, paired_lens):
        end = reads[-1][1]

        matches = []

        for i in range(len(reads)-1):
            start = reads[i][0]
            if end-start in paired_lens:
                matches.append(i)

        if matches:
            # If there is at least one exact match, return one of them at random
            i = random.choice(matches)

            pair = [reads[i], reads[-1]]

            l = end - reads[i][0]
            if paired_lens[l] > 1:
                paired_lens[l] -= 1
            else:
                del paired_lens[l]
            del reads[-1]
            del reads[i]

            return pair
        else:
            return None

    def findClosestLeftPair(self, reads, paired_lens):
        start = reads[0][0]

        closestD = None
        closestI = None
        closestL = None
        for i in range(1, len(reads)):
            l = reads[i][1] - start
            for pl in paired_lens:
                if closestD == None or abs(l-pl) < closestD:
                    closestD = abs(l-pl)
                    closestI = i
                    closestL = pl

        if closestD == None:
            print('Error! Trying to pair only 1 read?')
            exit()

        pair = [reads[0], reads[closestI]]
        del reads[closestI]
        del reads[0]
        return pair

    def findClosestRightPair(self, reads, paired_lens):
        end = reads[-1][1]

        closestD = None
        closestI = None
        closestL = None
        for i in range(len(reads)-1):
            l = end - reads[i][0]
            for pl in paired_lens:
                if closestD == None or abs(l-pl) < closestD:
                    closestD = abs(l-pl)
                    closestI = i
                    closestL = pl

        if closestD == None:
            print('Error! Trying to pair only 1 read?')
            exit()

        pair = [reads[closestI], reads[-1]]
        del reads[-1]
        del reads[closestI]
        return pair





    def findPairs(self, reads, pairedLens):
        if len(pairedLens) == 0:
            return reads, []

        for k,v in pairedLens.items():
            self.origPairedDepth += k * v

        length = 0
        for r in reads:
            if r[1] > length:
                length = r[1]

        paired = []
        unpaired = []

        reads.sort()

        # Sort read pairedLensSorted lengths
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


                    if reads[id2][0] < startRead[0]:
                        paired += [[reads[id2], startRead]]
                    else:
                        paired += [[startRead, reads[id2]]]

                    #paired += [[startRead, reads[id2]]]
                    del reads[id2]

                    break

            if not foundRead:
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


                        if reads[id2][0] < startRead[0]:
                            paired += [[reads[id2], startRead]]
                        else:
                            paired += [[startRead, reads[id2]]]

                        #paired += [[reads[id2], startRead]]
                        del reads[id2]
                        break

                if not foundRead:
                    unmatched += [startRead]

                while i < j and starts[i] <= 0:
                    i += 1
                while j > i and ends[j] <= 0:
                        j -= 1

        # Pair up remaining reads until we meet the quota of paired-end reads

        paired_lens_new = dict()
        for k,v in pairedLens.items():
            if v > 0:
                paired_lens_new[k] = v

        u, p = self.findPairsGreedy(unmatched+reads, paired_lens_new)

        unpaired = u
        paired += p


        for p in paired:
            self.calcPairedDepth += p[1][1] - p[0][0]

        return unpaired, paired

    def findClosestRead(self, reads, s):
        '''

        :param reads:
        :param start:
        :return: The read in reads that starts closest to s
        '''

        min_dist = None
        best_id = None

        for i in range(len(reads)):
            if reads[i][0] == s:
                return i
            else:
                d = abs(reads[i][0] - s)
                if not min_dist or d < min_dist:
                    min_dist = d
                    best_id = i
        if not best_id:
            print('No reads?!')
            exit()
        return best_id

    def findReads(self, unpairedLens, pairedLens, lensLeft, lensRight, coverage, boundaries=None, debug=False):
        ''' Find the set of reads that most closely matches the distribution of readLens and the coverage vector
        '''

        if debug:
            print('Unpaired:')
            print(unpairedLens)
            print('Paired:')
            print(pairedLens)
            print('Left:')
            print(lensLeft)
            print('Right:')
            print(lensRight)
            print('  -->')

        fragmentLens = copy.copy(unpairedLens)
        for k,v in lensLeft.items():
            if k > 0:
                if k in fragmentLens:
                    fragmentLens[k] += v
                else:
                    fragmentLens[k] = v
        for k,v in lensRight.items():
            if k > 0:
                if k in fragmentLens:
                    fragmentLens[k] += v
                else:
                    fragmentLens[k] = v


        if debug:
            countUnpaired = 0
            for k,v in unpairedLens.items():
                countUnpaired += v
            countPaired = 0
            for k,v in pairedLens.items():
                countPaired += v

        reads = self.findReadsInCoverage_v1(coverage, fragmentLens)

        if debug:
            print('Reads:')
            print(reads)

        if boundaries:
            unpaired, paired = self.findPairsWithBoundaries(reads, unpairedLens, pairedLens, boundaries)
        else:
            #unpaired, paired = self.findPairs(reads, pairedLens)
            unpaired, paired = self.findPairsRandom(reads, pairedLens)

        if debug:
            print('  Unpaired: %d\t->\t%d' % (countUnpaired, len(unpaired)))
            print('  Paired:   %d\t->\t%d' % (countPaired, len(paired)))
            print('')


        if boundaries:
            start = boundaries[0]
            end = boundaries[-1]
            for r in unpaired:
                if r[0] >= start or r[1] <= end:
                    self.badUnpaired += 1
                self.countUnpaired += 1
            for p in paired:
                if p[0][0] >= start or p[1][1] <= end:
                    self.badPaired += 1
                self.countPaired += 1

        return unpaired, paired

    def findReadsWithPairs(self, unpairedLens, pairs, lensLeft, lensRight, coverage, debug=False):
        ''' Find the set of reads that most closely matches the distribution of readLens and the coverage vector
        '''

        countPaired = 0
        for k,v in lensLeft.items():
            countPaired += v

        fragmentLens = copy.copy(unpairedLens)
        for k,v in lensLeft.items():
            if k > 0:
                if k in fragmentLens:
                    fragmentLens[k] += v
                else:
                    fragmentLens[k] = v
        for k,v in lensRight.items():
            if k > 0:
                if k in fragmentLens:
                    fragmentLens[k] += v
                else:
                    fragmentLens[k] = v

        unpaired, paired = self.findReadsInCoverageWithPairs(coverage, unpairedLens, lensLeft, lensRight, pairs)
        unpaired.sort()
        paired.sort()
        return unpaired, paired

    def findReadsInCoverageWithPairs(self, coverage, unpairedLens, lensLeft, lensRight, pairs):
        if len(pairs) > 0:
            paired = []

            # Reads sorted by frequency, most to least common
            lensSortedLeft = sorted(lensLeft, key=lensLeft.get, reverse=True)
            lensSortedRight = sorted(lensRight, key=lensRight.get, reverse=True)

            # Longest read length
            longestLeft = max(lensSortedLeft)
            longestRight = max(lensSortedRight)

            pairs = sorted(pairs, key = lambda p: p[1]-p[0])

            for (start,end) in pairs:
                pairLength = end - start

                if coverage[start] == 0 or coverage[end-1] == 0:
                    if self.debug:
                        print('Endpoint is 0!')
                    continue

                ''' Find starting read '''
                maxLen = 1
                while maxLen < longestLeft and start+maxLen < len(coverage) and coverage[start+maxLen] > 0:
                    maxLen += 1
                if maxLen > pairLength:
                    maxLen = pairLength

                startLen = None
                # Find a length that ends on a step down
                for l in lensSortedLeft:
                    if l <= maxLen and (start+l == len(coverage) or coverage[start+l-1] > coverage[start+l]):
                        startLen = l
                        break
                # Find any length that fits in coverage
                if not startLen:
                    for l in lensSortedLeft:
                        if l <= maxLen and start+l < len(coverage):
                            startLen = l
                            break
                # Find any remaining length
                if not startLen:
                    startLen = maxLen
                else:
                    lensLeft[startLen] -= 1
                    if lensLeft[startLen] == 0:
                        i = 0
                        while not lensSortedLeft[i] == startLen:
                            i += 1
                        del lensSortedLeft[i]

                for i in range(start, start+startLen):
                    coverage[i] -= 1


                ''' Find ending read'''
                maxLen = 0
                while maxLen < longestRight and end > maxLen and coverage[end-maxLen-1] > 0:
                    maxLen += 1
                if maxLen > pairLength:
                    maxLen = pairLength

                endLen = None
                # Find a length that starts on a step up
                for l in lensSortedRight:
                    if l <= maxLen and (end == l+1 or coverage[end-l] > coverage[end-l-1]):
                        endLen = l
                        break
                # Find any length that fits in coverage
                if not endLen:
                    for l in lensSortedRight:
                        if l <= maxLen and end-l >= 0:
                            endLen = l
                            break

                if not endLen:
                    endLen = maxLen
                else:
                    lensRight[endLen] -= 1
                    if lensRight[endLen] == 0:
                        i = 0
                        while not lensSortedRight[i] == endLen:
                            i += 1
                        del lensSortedRight[i]

                for i in range(end-endLen, end):
                    coverage[i] -= 1


                paired.append([[start, start+startLen], [end-endLen, end]])
        else:
            paired = []

        if len(unpairedLens) > 0:
            unpaired = self.findReadsInCoverage_v1(coverage, unpairedLens)
        else:
            unpaired = []

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
            #while boundaries[boundBottom] <= start:
            #    boundBottom += 1

            # find a read from the beginning
            readStart = start
            readEnd = start

            closestEndpoint = None
            for length in range(1, maxLen+1):
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
                #while boundaries[boundTop] >= end:
                #    boundTop -= 1

                # find a read from the end
                readEnd = end
                readStart = end

                closestEndpoint = None
                for length in range(1, maxLen+1):
                #for length in range(1, maxLen+1):
                    if (readEnd-length == start) or (readEnd-length > start and coverage[readEnd - length] > coverage[readEnd - length - 1]):
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

    def findSpliceSites(self, start, end):
        ''' Look for any exon boundaries crossed by the given interval. 
            Return a list of all splice sites in the interval [start, end)
        '''

        startExon = bisect.bisect_right(self.exons, start)
        endExon = startExon
        while endExon < len(self.exons) and self.exons[endExon] < end:
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

    def processRead(self, read, name, paired=False):
        ''' If read is unpaired, add it to the correct spliced or unspliced list of reads.
            If read is paired, find its pair or add it to a list to be found later. Once a pair of reads is found, add the combined read to the appropriate list of reads
        '''

        # Update read index based on chromosome
        offset = self.chromOffsets[read.chrom]
        for i in range(len(read.exons)):
            read.exons[i] = [read.exons[i][0]+offset, read.exons[i][1]+offset]

        # update list of exons
        alignment = read.exons
        if len(alignment) > 1:
            for i in range(len(alignment)-1):
                self.exons.add(alignment[i][1])
                self.exons.add(alignment[i+1][0])

        #self.updateGeneBounds([read.exons[0][0], read.exons[-1][1]])

        if not paired:
            self.updateGeneBounds([read.exons[0][0], read.exons[-1][1]])

            # unpaired read
            if len(read.exons) == 1:
                self.addUnspliced(read)
            else:
                self.addSpliced(read)
        else:
            read.pairOffset += offset
            if name in self.unmatched:
                foundMatch = False
                for i in range(len(self.unmatched[name])):
                    match = self.unmatched[name][i]
                    if read.pairOffset == match.exons[0][0] and match.pairOffset == read.exons[0][0] and not self.conflicts(read.exons, match.exons):
                        xs = read.xs or match.xs
                        NH = min(read.NH, match.NH)

                        self.paired.append(pairedread.PairedRead(match.chrom, match.exons, read.chrom, read.exons, xs, NH))

                        start = min(match.exons[0][0], read.exons[0][0])
                        end = max(match.exons[-1][1], read.exons[-1][1])
                        if end-start < 100000:
                            self.updateGeneBounds([start, end])

                        # Update NH values
                        if match.NH <= read.NH:
                            del self.unmatched[name][i]

                            if match.NH < read.NH:
                                read.NH -= match.NH
                                self.unmatched[name].append(read)
                        else:
                            match.NH -= read.NH

                        foundMatch = True
                        break

                if not foundMatch:
                    self.unmatched[name].append(read)
            else:
                self.unmatched[name] = [read]

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
                    # Exon in B overlaps an intron in A
                    return 1

        countA = len(exonsA)
        for i in range(countA):
            e = exonsA[countA-i-1]
            if e[1] < exonsB[0][1]:
                break

            for i in range(len(exonsB)-1):
                if e[1] <= exonsB[i][1]:
                    break
                elif e[1] > exonsB[i][1] and e[0] < exonsB[i+1][0]:
                    # Exon in A overlaps an intron in B
                    return 2

        return 0


    def updateGeneBounds(self, segment):

        '''
        if segment[1] - segment[0] > 1000000:
            print(segment)

        if len(self.gene_bounds) == 0:
            self.gene_bounds = [segment]
            return

        j = len(self.gene_bounds)
        while j > 0 and self.gene_bounds[j-1][0] - segment[1] > self.overlap_radius:
            j -= 1
        i = j
        while i > 0 and segment[0] - self.gene_bounds[i-1][1] <= self.overlap_radius:
            i -= 1

        if i < len(self.gene_bounds):
            segment[0] = min(segment[0], self.gene_bounds[i][0])
        if j > 0:
            segment[1] = max(segment[1], self.gene_bounds[j-1][1])

        if j == i+1:
            self.gene_bounds[i] = segment
        else:
            self.gene_bounds = self.gene_bounds[:i] + [segment] + self.gene_bounds[j:]
        '''


        i = bisect.bisect_left(self.gene_bounds, segment)

        j = i
        while j > 0 and segment[0] - self.gene_bounds[j-1][1] < self.overlap_radius:
            j -= 1
        k = i
        while k < len(self.gene_bounds) and self.gene_bounds[k][0] - segment[1] < self.overlap_radius:
            k += 1


        new_gene = segment[:]
        if j < i:
            new_gene[0] = min(self.gene_bounds[j][0], new_gene[0])
        if k > i:
            new_gene[1] = max(self.gene_bounds[k-1][1], segment[1])

        self.gene_bounds = self.gene_bounds[:j] + [new_gene] + self.gene_bounds[k:]


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

    def countReads(self, d):
        count = 0
        for k,v in d.items():
            count += v
        return count












    '''
        Outdated/unused functions
    '''

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






    def findPairsDumb(self, reads, unpairedLens, pairedLens):

        for k,v in pairedLens.items():
            self.origPairedDepth += k * v

        countPaired = 0
        for k,v in pairedLens.items():
            countPaired += v

        paired = []
        i = 0
        j = len(reads)-1
        for _ in range(countPaired):
            if i >= j:
                break

            paired += [[reads[i], reads[j]]]
            i += 1
            j -= 1

        for p in paired:
            self.calcPairedDepth += p[1][1] - p[0][0]

        return reads[i:j+1], paired

    def findPairs2(self, length, reads, unpairedLens, pairedLens):
        '''
        starts = [0] * (length+1)
        ends = [0] * (length+1)

        for r in reads:
            starts[r[0]] += 1
            ends[r[1]] += 1

        print('Finding pairs')
        startTime = time.time()
        unpaired, paired, score = self.findPairsRecursive(starts, ends, 0, reads, unpairedLens, pairedLens)
        endTime = time.time()
        print('Finished in %0.2fs, score = %d' % (endTime-startTime, score))

        return unpaired, paired
        '''

        reads.sort()

        dists = [0] * len(reads)
        for i in range(len(reads)):
            dists[i] = [0] * len(reads)
            for j in range(i+1):
                dists[i][j] = abs(reads[j][1] - reads[i][0])

        print('Finding pairs from %d reads' % len(reads))
        startTime = time.time()
        unpaired, paired, score = self.findPairsRecursive(reads, dists, unpairedLens, pairedLens, currentScore=0, maxScore=len(reads))
        endTime = time.time()
        print('Finished in %0.2fs, score = %d' % (endTime-startTime, score))

        return unpaired, paired

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

    def findPairsRecursive(self, reads, dists, unpairedLens, pairedLens, currentScore=0, maxScore=None, numUnpaired=None, numPaired=None):
        if numUnpaired == None:
            numUnpaired = 0
            for l,n in  unpairedLens.items():
                numUnpaired += n
        if numPaired == None:
            numPaired = 0
            for l,n in pairedLens.items():
                numPaired += n

        numReads = len(reads)
        if maxScore == None:
            maxScore = numReads
        if currentScore >= maxScore:
            return [], [], maxScore

        if numPaired == 0 and numUnpaired == 0:
            return [], [], currentScore

        bestScore = None
        bestUnpaired = []
        bestPaired = []

        if numPaired > 0:
            for l in pairedLens:
                if pairedLens[l] == 0:
                    continue

                for i in range(numReads):
                    for j in range(i):
                        if dists[i][j] == l:
                            # Update distance matrix
                            #newDists = [r[:j]+r[j+1:i]+r[i+1:] for r in dists[:j]+dists[j+1:i]+dists[i+1:]]
                            newDists = dists[:][:]
                            for x in range(numReads):
                                newDists[i][x] = -1
                                newDists[x][j] = -1

                            pairedLens[l] -= 1

                            #unpaired, paired, score = self.findPairsRecursive(reads[:j] + reads[j+1:i] + reads[i+1:], newDists, unpairedLens, pairedLens, currentScore, maxScore)
                            unpaired, paired, score = self.findPairsRecursive(reads, newDists, unpairedLens, pairedLens, currentScore, maxScore, numUnpaired, numPaired-1)

                            if score == 0:
                                return unpaired, paired + [reads[j], reads[i]], score
                            elif bestScore == None or score < bestScore:
                                bestScore = score
                                bestUnpaired = unpaired
                                bestPaired = paired + [reads[j], reads[i]]
                                maxScore = bestScore

                            pairedLens[l] += 1

        if numUnpaired > 0:
            for l in unpairedLens:
                if unpairedLens[l] == 0:
                    continue

                for i in range(len(reads)):
                    if dists[i][i] == l:
                        #newDists = [r[:i] + r[i+1:] for r in dists[:i]+dists[i+1:]]
                        newDists = dists[:][:]
                        for x in range(numReads):
                            newDists[i][x] = -1
                            newDists[x][i] = -1

                        unpairedLens[l] -= 1

                        #unpaired, paired, score = self.findPairsRecursive(reads[:i] + reads[i+1:], newDists, unpairedLens, pairedLens, currentScore, maxScore)
                        unpaired, paired, score = self.findPairsRecursive(reads, newDists, unpairedLens, pairedLens, currentScore, maxScore, numUnpaired-1, numPaired)

                        if score == 0:
                            return unpaired + [reads[i]], paired, score
                        elif bestScore == None or score < bestScore:
                            bestScore = score
                            bestUnpaired = unpaired + [reads[i]]
                            bestPaired = paired
                            maxScore = bestScore

                        unpairedLens[l] += 1

        # No lengths match
        if numPaired > 0:
            if bestScore == None:
                for l in pairedLens:
                    if pairedLens[l] == 0:
                        continue

                    for i in range(numReads):
                        for j in range(i):
                            # Update distance matrix
                            #newDists = [r[:j]+r[j+1:i]+r[i+1:] for r in dists[:j]+dists[j+1:i]+dists[i+1:]]
                            newDists = dists[:][:]
                            for x in range(numReads):
                                newDists[i][x] = -1
                                newDists[x][j] = -1

                            pairedLens[l] -= 1

                            #unpaired, paired, score = self.findPairsRecursive(reads[:j] + reads[j+1:i] + reads[i+1:], newDists, unpairedLens, pairedLens, currentScore+1, maxScore)
                            unpaired, paired, score = self.findPairsRecursive(reads, newDists, unpairedLens, pairedLens, currentScore+1, maxScore, numUnpaired, numPaired-1)

                            if score == 0:
                                return unpaired, paired + [reads[j], reads[i]], score
                            elif bestScore == None or score < bestScore:
                                bestScore = score
                                bestUnpaired = unpaired
                                bestPaired = paired + [reads[j], reads[i]]
                                maxScore = bestScore

                            pairedLens[l] += 1

            if numUnpaired > 0:
                for l in unpairedLens:
                    if unpairedLens[l] == 0:
                        continue

                    for i in range(len(reads)):
                        #newDists = [r[:i] + r[i+1:] for r in dists[:i]+dists[i+1:]]
                        newDists = dists[:][:]
                        for x in range(numReads):
                            newDists[i][x] = -1
                            newDists[x][i] = -1

                        unpairedLens[l] -= 1

                        #unpaired, paired, score = self.findPairsRecursive(reads[:i] + reads[i+1:], newDists, unpairedLens, pairedLens, currentScore+1, maxScore)
                        unpaired, paired, score = self.findPairsRecursive(reads, newDists, unpairedLens, pairedLens, currentScore+1, maxScore, numUnpaired-1, numPaired)

                        if score == 0:
                            return unpaired + [reads[i]], paired, score
                        elif bestScore == None or score < bestScore:
                            bestScore = score
                            bestUnpaired = unpaired + [reads[i]]
                            bestPaired = paired
                            maxScore = bestScore

                        unpairedLens[l] += 1

        if (not bestUnpaired or len(bestUnpaired) == 0) and (not bestPaired or len(bestPaired) == 0):
            print('Error!')
            print(unpairedLens)
            print(pairedLens)
            print(reads)
            print(bestUnpaired)
            print(bestPaired)
            print(bestScore)
            exit()
        return bestUnpaired, bestPaired, bestScore


    '''
    def findPairsRecursive(self, starts, ends, i, reads, unpairedLens, pairedLens):
        # Test if reads list is empty
        if not reads:
            return [], [], 0

        while i < len(starts) and starts[i] == 0:
            i += 1

        #print('')
        #print(reads)
        #print(i)
        #print(starts[i])

        starts[i] -= 1

        bestScore = None
        bestUnpaired = None
        bestPaired = None

        # Test all unpaired lengths
        testedLengths = []
        for rid in range(len(reads)):
            r = reads[rid]
            if r[0] == i:
                l = r[1] - r[0]
                if l in unpairedLens and unpairedLens[l] > 0 and not l in testedLengths:
                    #print(r)
                    ends[r[1]] -= 1
                    unpairedLens[l] -= 1

                    # Recurse
                    unpaired, paired, score = self.findPairsRecursive(starts, ends, i, reads[:rid]+reads[rid+1:], unpairedLens, pairedLens)

                    ends[r[1]] += 1
                    unpairedLens[l] += 1
                    testedLengths.append(l)

                    if score == 0:
                        return unpaired, paired, score
                    elif bestScore == None or score < bestScore:
                        bestScore = score
                        bestUnpaired = unpaired + [r]
                        bestPaired = paired

        # Test all paired lengths
        id1 = 0
        while not reads[id1][0] == i:
            id1 += 1
        startRead = reads[id1]
        ends[startRead[1]] -= 1
        newReads = reads[:id1] + reads[id1+1:]

        for l in pairedLens:
            if pairedLens[l] == 0 or i+l >= len(starts):
                continue

            if ends[i+l] > 0:
                id2 = 0
                while not newReads[id2][1] == i+l:
                    id2 += 1
                endRead = newReads[id2]
                starts[endRead[0]] -= 1
                ends[i+l] -= 1
                pairedLens[l] -= 1

                #print(str(startRead) + ', ' + str(endRead))
                unpaired, paired, score = self.findPairsRecursive(starts, ends, i, newReads[:id2]+newReads[id2+1:], unpairedLens, pairedLens)

                starts[endRead[0]] += 1
                ends[i+l] += 1
                pairedLens[l] += 1

                if score == 0:
                    return unpaired, paired, score
                elif bestScore == None or score < bestScore:
                    bestScore = score
                    bestUnpaired = unpaired
                    bestPaired = paired + [[startRead, endRead]]

        if bestScore == None:
            for rid in range(len(newReads)):
                endRead = newReads[rid]
                starts[endRead[0]] -= 1
                ends[endRead[1]] -= 1

                unpaired, paired, score = self.findPairsRecursive(starts, ends, i, newReads[:rid]+newReads[rid+1:], unpairedLens, pairedLens)

                starts[endRead[0]] += 1
                ends[endRead[1]] += 1

                # increment score since the current length is a mismatch
                score += 1

                if score == 1:
                    return unpaired, paired, score
                elif bestScore == None or score < bestScore:
                    bestScore = score
                    bestUnpaired = unpaired
                    bestPaired = paired + [[startRead, endRead]]

        return bestUnpaired, bestPaired, bestScore
    '''

    '''
    def findPairsRecursive(self, starts, ends, reads, unpairedLens, pairedLens, pairedLensSorted, i, currentScore=0, scoreCutoff=None, currentDepth=0, maxDepth=None):

        if depthCutoff == None:
            depthCutoff = len(reads)

        if not (maxDepth == None) and currentDepth >= maxDepth:
            return [],[],0

        if i == len(starts):
            return [], [], 0

        if starts[i] == 0:
            while i < len(starts) and starts[i] == 0:
                i += 1
            return self.findPairsRecursive(starts, ends, reads, readLens, readLensSorted, i, currentScore, scoreCutoff, currentDepth, maxDepth)


        bestScore = None
        bestRead = None
        bestPaired = None
        bestUnpaired = None
        bestSingle = None

        # Test all unpaired lengths
        testedLengths = []
        for rid in range(len(reads)):
            read = reads[rid]
            if read[0] == i:
                l = read[1] - read[0]
                if not l in testedLengths:
                    starts[i] -= 1
                    ends[i+l] -= 1
                    del reads[rid]
                    unpairedLens[l] -= 1

                    # recurse
                    newi = i
                    while newi < len(starts) and starts[newi] == 0:
                        newi += 1
                    unpaired, paired, score = self.findPairsRecursive(starts, ends, reads, readLens, readLensSortedNew, newi, currentScore, scoreCutoff, currentDepth+1, maxDepth)

                    reads.append(read)
                    starts[read[0]] += 1
                    ends[read[1]] += 1
                    unpairedLens[l] += 1

                    if score == 0:
                        return unpaired+[start+[read], paired, 0
                    elif score < bestScore:
                        bestRead = read
                        bestScore = score
                        bestPaired = paired
                        bestUnpaired = unpaired
                        bestSingle = single

                    if score < scoreCutoff:
                        scoreCutoff = score

                    testedLengths.append(l)

        # Test all paired lengths
        id1 = 0
        while not reads[id1][0] == i:
            id1 += 1
            read = reads[id1]
            starts[i] -= 1
            ends[read[1]] -= 1
            del reads[id1]

        for l in range(len(pairedLensSorted)):
            if pairedLens[l] == 0:
                continue

            if ends[i+l] == 0 and

            if ends[i+l] > 0:
                id2 = 0
                while not reads[id2][1] == i+l:
                    id2 += 1
                read = [startRead, reads[id2]]
                starts[read[1][0]] -= 1
                ends[id2] -= 1
                del reads[id2]

                # update readLens
                pairedLens[l] -= 1

                # recurse
                newi = i
                while newi < len(starts) and starts[newi] == 0:
                    newi += 1
                unpaired, paired, score = self.findPairsRecursive(starts, ends, reads, readLens, readLensSortedNew, newi, currentScore, scoreCutoff, currentDepth+1, maxDepth)

                reads.append(read[1])
                starts[read[1][0]] += 1
                ends[read[1][1]] += 1

                pairedLens[l] += 1

                if score == 0:
                    return unpaired, paired+[read], 0
                elif score < bestScore:
                    bestRead = read
                    bestScore = score
                    bestPaired = paired
                    bestUnpaired = unpaired
                    bestSingle = single

                if score < scoreCutoff:
                    scoreCutoff = score
        if bestScore == None:
            if currentScore >= scoreCutoff:
                return [],[],1

            read = reads[id1]
            del reads[id1]

            # update starts and ends
            starts[read[0]] -= 1
            ends[read[1]] -= 1

            newi = i
            while newi < len(starts) and starts[newi] == 0:
                newi += 1
            unpaired, paired, score = self.findPairsRecursive(starts, ends, reads, readLens, readLensSorted, newi, currentScore+1, scoreCutoff, currentDepth+1, maxDepth)


            reads.append(read)
            starts[read[0]] += 1
            ends[read[1]] += 1
            return unpaired, paired, score+1
        else:
            if bestSingle:
                return bestUnpaired+[bestRead], bestPaired, bestScore
            else:
                return bestUnpaired, bestPaired+[bestRead], bestScore
    '''
