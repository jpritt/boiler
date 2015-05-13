#! /usr/bin/env python

import time
import bisect

class Pairs:
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






    def findPairsGreedy(self, reads, unpairedLens, pairedLens):
        origReads = reads[:]
        origPaired = dict()
        for k,v in pairedLens.items():
            origPaired[k] = v
        origUnpaired = dict()
        for k,v in unpairedLens.items():
            origUnpaired[k] = v

        origPairedDepth = 0
        for k,v in pairedLens.items():
            #self.origPairedDepth += k * v
            origPairedDepth += k * v

        pairedLensSorted = sorted(pairedLens.keys(), reverse=True)

        reads.sort()
        numReads = len(reads)

        # Keep track of which reads we have already assigned
        assigned = [0] * numReads

        # Map each distance to the paired reads that match it
        pairDists = dict()
        #singleDists = dict()

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

        #print('\n'.join(['\t'.join([str(a) for a in r]) for r in dists]))
        print(pairedLens)

        paired = []
        unpaired = []
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

                    print('Found unique length %d' % l)
                    print('\t(%d,%d) - (%d,%d)' % (reads[j][0], reads[j][1], reads[i][0], reads[i][1]))

                    assigned[i] = 1
                    assigned[j] = 1

                    if dists[i][i] > 0:
                        #singleDists[dists[i][i]] -= 1
                        dists[i][i] = 0
                    if dists[j][j] > 0:
                        #singleDists[dists[j][j]] -= 1
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
                    #print('Best pair: ' + str(bestL))
                    pairedLens[bestL] -= 1
                    i,j = self.findPairedDist(dists, assigned, numReads, bestL)
                    paired.append([reads[j], reads[i]])

                    print('Found duplicate length %d (%d)' % (bestL, bestFreq))
                    print('\t(%d,%d) - (%d,%d)' % (reads[j][0], reads[j][1], reads[i][0], reads[i][1]))

                    assigned[i] = 1
                    assigned[j] = 1

                    if dists[i][i] > 0:
                        #singleDists[dists[i][i]] -= 1
                        dists[i][i] = 0
                    if dists[j][j] > 0:
                        #singleDists[dists[j][j]] -= 1
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
            print('')

        print('%d unmatched' % countPaired)
        print(pairedLens)
        print('\t'.join([str(r) for r in reads]))
        print('\t'.join([str(r) for r in assigned]))

        remaining = [0] * (numReads - sum(assigned))
        i = 0
        for j in range(numReads):
            if not assigned[j]:
                remaining[i] = reads[j]
                i += 1

        i = 0
        j = len(remaining)-1
        for _ in range(countPaired):
            if i >= j:
                break

            paired.append([remaining[i], remaining[j]])
            i += 1
            j -= 1

        calcPairedDepth = 0
        for p in paired:
            #self.calcPairedDepth += p[1][1] - p[0][0]
            calcPairedDepth += p[1][1] - p[0][0]

        print('%d\t-->\t%d' % (origPairedDepth, calcPairedDepth))

        '''
        if not calcPairedDepth == origPairedDepth:
            print('%d\t-->\t%d' % (origPairedDepth, calcPairedDepth))
            print(origReads)
            print(origPaired)
            print(origUnpaired)
            print('-->')
            print(paired)
            print(unpaired+remaining[i:j+1])
            print('')
        '''

        return unpaired+remaining[i:j+1], paired






    def findPairsGreedy2(self, reads, unpairedLens, pairedLens):
        origReads = reads[:]
        origPaired = dict()
        for k,v in pairedLens.items():
            origPaired[k] = v
        origUnpaired = dict()
        for k,v in unpairedLens.items():
            origUnpaired[k] = v

        origPairedDepth = 0
        for k,v in pairedLens.items():
            origPairedDepth += k * v
            #self.origPairedDepth += k * v

        reads.sort()
        pairedLensSorted = sorted(pairedLens.keys())

        paired = []

        countPaired = 0
        for k,v in pairedLens.items():
            countPaired += v

        # Create a distance matrix between all pairs of reads
        numReads = len(reads)
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
                    if bestL == None or (l > 0 and diff < bestDiff):
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
                paired += [[reads[i], reads[j]]]
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


        countUnpaired = numReads - 2 * countPaired
        unpaired = [0] * countUnpaired
        i = 0
        for j in range(numReads):
            if not assigned[j]:
                unpaired[i] = reads[j]
                i += 1

        calcPairedDepth = 0
        for p in paired:
            calcPairedDepth += p[1][1] - p[0][0]
            #self.calcPairedDepth += p[1][1] - p[0][0]

        if float(calcPairedDepth) / float(origPairedDepth) > 1.15:
            print('%d\t-->\t%d' % (origPairedDepth, calcPairedDepth))
            print(origReads)
            print(origPaired)
            print(origUnpaired)
            print('-->')
            print(paired)
            print(unpaired)

        return unpaired, paired

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

    def findPairs(self, reads, unpairedLens, pairedLens):
        length = 0
        for r in reads:
            if r[1] > length:
                length = r[1]

        origPairedDepth = 0
        for k,v in pairedLens.items():
            origPairedDepth += k * v

        origReads = reads[:]
        origPaired = dict()
        for k,v in pairedLens.items():
            origPaired[k] = v
        origUnpaired = dict()
        for k,v in unpairedLens.items():
            origUnpaired[k] = v

        paired = []
        unpaired = []

        reads.sort()

        if len(pairedLens) > 0:
            # Sort read pairedLensSorted lengths
            pairedLensSorted = sorted(pairedLens)

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


                            if reads[id2][0] < startRead[0]:
                                paired += [[reads[id2], startRead]]
                            else:
                                paired += [[startRead, reads[id2]]]

                            #paired += [[reads[id2], startRead]]
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
                    if len(unmatched) > 1:
                        if unmatched[-1][0] < unmatched[0][0]:
                            paired += [[unmatched[-1], unmatched[0]]]
                        else:
                            paired += [[unmatched[0], unmatched[-1]]]
                        del unmatched[-1]
                        del unmatched[0]

            # Add remaining unmatched reads as unpaired reads
            for r in unmatched:
                unpaired += [r]

        # Add remaining reads as unpaired reads
        for r in reads:
            unpaired += [r]

        calcPairedDepth = 0
        for p in paired:
            calcPairedDepth += p[1][1] - p[0][0]

        print('%d -> %d' % (origPairedDepth, calcPairedDepth))

        return unpaired, paired

reads = [[298, 374], [531, 607], [610, 686], [1250, 1326], [1577, 1653], [2019, 2095], [2101, 2177], [2144, 2220], [2192, 2268], [2429, 2505], [2429, 2505], [2739, 2815], [4358, 4434], [4474, 4550], [5030, 5106], [7254, 7330], [7272, 7348], [7281, 7357], [7323, 7399], [7480, 7556]]
unpairedLens = {76:14}
pairedLens = {328: 1, 277: 1, 134: 1}

p = Pairs()
#startTime = time.time()
unpaired, paired = p.findPairsGreedy2(reads, unpairedLens, pairedLens)
#endTime = time.time()
#print('%0.3f s' % (endTime-startTime))

print('')
print(unpaired)
print(paired)