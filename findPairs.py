import time
import bisect
import random

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
        origCountPaired = countPaired
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
                    #print('Best pair: ' + str(bestL))
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

        #print('%d\t/\t%d' % (countPaired, origCountPaired))


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

        '''
        i = 0
        j = len(remaining)-1
        for _ in range(countPaired):
            if i >= j:
                break

            paired.append([remaining[i], remaining[j]])
            i += 1
            j -= 1
        '''

        if countPaired > 0:
            newUnpaired, newPaired = self.findPairsGreedy2(remaining, unpairedLens, remainingPairedLens)
            unpaired = newUnpaired
            paired += newPaired
        else:
            unpaired = remaining


        calcPairedDepth = 0
        for p in paired:
            #self.calcPairedDepth += p[1][1] - p[0][0]
            calcPairedDepth += p[1][1] - p[0][0]

        '''
        if countPaired >= 10: #not calcPairedDepth == origPairedDepth:
            print('%d\t/\t%d' % (countPaired, origCountPaired))
            print('%d\t-->\t%d' % (origPairedDepth, calcPairedDepth))
            print(origReads)
            print(origPaired)
            print(origUnpaired)
            print('-->')
            print(paired)
            print(unpaired+remaining[i:j+1])
            print('')
        '''


        #return remaining[i:j+1], paired
        return unpaired, paired




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

        numReads = len(reads)
        reads.sort()
        pairedLensSorted = sorted(pairedLens.keys(), reverse=True)

        paired = []

        countPaired = 0
        for k,v in pairedLens.items():
            countPaired += v
        countUnpaired = numReads - 2 * countPaired
        #print(countPaired)

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

        #print(countPaired)
        #print(countUnpaired)
        #print(len(assigned))
        #print(sum(assigned))
        unpaired = [0] * countUnpaired
        i = 0
        for j in range(numReads):
            if not assigned[j]:
                unpaired[i] = reads[j]
                i += 1
        #print('')

        calcPairedDepth = 0
        for p in paired:
            calcPairedDepth += p[1][1] - p[0][0]

        print('%d / %d' % (calcPairedDepth, origPairedDepth))


        '''
        if float(calcPairedDepth) / float(origPairedDepth) < 0.9:
            print('%d\t-->\t%d' % (origPairedDepth, calcPairedDepth))
            print(origReads)
            print(origPaired)
            print(origUnpaired)
            print('-->')
            print(paired)
            print(unpaired)
            print('')
        '''

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

        #if len(reads) > 0 and len(unmatched) > 0:
        #    print('Reads and unmatched both left!')
        #    exit()

        print('Found %d pairs, still looking for %d pairs, %d reads left and %d unmatched' % (len(paired), countPairs, len(reads), len(unmatched)))

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

reads = [[1263, 1339], [6404, 6445], [3595, 3671], [6339, 6415], [3595, 3671], [6312, 6388], [3597, 3673], [5970, 6046], [3797, 3873], [5969, 6045], [3818, 3894], [5859, 5935], [3820, 3896], [5765, 5841], [3869, 3945], [5724, 5800], [3928, 4004], [5577, 5653], [3960, 4036], [5501, 5577], [4076, 4152], [5496, 5572], [4137, 4213], [5481, 5557], [4211, 4287], [5340, 5416], [4219, 4295], [5332, 5408], [4327, 4403], [5318, 5394], [4382, 4458], [5297, 5373], [4447, 4523], [5264, 5340], [4459, 4535], [5242, 5318], [4511, 4587], [5109, 5185], [4520, 4596], [5017, 5093], [4531, 4607], [5013, 5089], [4546, 4622], [5010, 5086], [4576, 4652], [5010, 5086], [4581, 4657], [4940, 5016], [4598, 4674], [4810, 4886], [4641, 4717], [4784, 4860], [4669, 4745], [4744, 4820], [4692, 4768], [4710, 4786], [4707, 4783]]
unpairedLens = {41:1, 76: 12}
pairedLens = {171: 1, 391: 1, 393: 1, 206: 1, 493: 1, 401: 1, 404: 1, 239: 1, 350: 1, 224: 1, 225: 1, 163: 1, 358: 1, 299: 1, 172: 1, 237: 1, 623: 1, 304: 1, 446: 1, 379: 1, 254: 1, 319: 1}
lensLeft = {76: 22}
lensRight = {76: 22}

countUnpaired = 0
for k,v in unpairedLens.items():
    countUnpaired += v
countPaired = 0
for k,v in pairedLens.items():
    countPaired += v
print('Searching for %d unpaired, %d paired' % (countUnpaired, countPaired))

p = Pairs()
#startTime = time.time()
unpaired, paired = p.findPairsRandom(reads, pairedLens)
#endTime = time.time()
#print('%0.3f s' % (endTime-startTime))

#print('')
#print(sorted(unpaired))
#print('')
#print(sorted(paired))


print('Found %d unpaired, %d paired' % (len(unpaired), len(paired)))