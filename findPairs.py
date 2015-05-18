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

reads = [[6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976955, 6977031], [6976956, 6977032], [6976965, 6977041], [6976987, 6977063], [6977054, 6977130], [6977072, 6977148], [6977074, 6977150], [6977075, 6977151], [6977085, 6977161], [6977098, 6977174], [6977101, 6977177], [6977118, 6977194], [6977118, 6977194], [6977120, 6977196], [6977124, 6977200], [6977149, 6977225], [6977153, 6977229], [6977161, 6977237], [6977167, 6977243], [6977170, 6977246], [6977176, 6977252], [6977179, 6977255], [6977184, 6977260], [6977184, 6977260], [6977187, 6977263], [6977190, 6977266], [6977193, 6977269], [6977197, 6977273], [6977200, 6977276], [6977200, 6977276], [6977201, 6977277], [6977206, 6977282], [6977216, 6977292], [6977227, 6977303], [6977230, 6977306], [6977233, 6977309], [6977238, 6977314], [6977242, 6977318], [6977253, 6977329], [6977254, 6977330], [6977255, 6977331], [6977259, 6977335], [6977263, 6977339], [6977278, 6977354], [6977282, 6977358], [6977283, 6977359], [6977283, 6977359], [6977290, 6977366], [6977293, 6977369], [6977299, 6977375], [6977304, 6977380], [6977304, 6977380], [6977309, 6977385], [6977319, 6977395], [6977326, 6977402], [6977328, 6977404], [6977332, 6977408], [6977346, 6977422], [6977348, 6977424], [6977351, 6977427], [6977356, 6977432], [6977365, 6977441], [6977368, 6977444], [6977370, 6977446], [6977382, 6977458], [6977383, 6977459], [6977383, 6977459], [6977384, 6977460], [6977388, 6977464], [6977419, 6977495], [6977423, 6977499], [6977428, 6977504], [6977432, 6977508], [6977433, 6977509], [6977446, 6977522], [6977452, 6977528], [6977487, 6977563], [6977503, 6977579], [6977508, 6977584], [6977509, 6977585], [6977516, 6977592], [6977524, 6977600], [6977527, 6977603], [6977527, 6977603], [6977543, 6977619], [6977554, 6977630], [6977579, 6977655], [6977582, 6977658], [6977584, 6977660], [6977591, 6977667], [6977595, 6977671], [6977596, 6977672], [6977606, 6977682], [6977612, 6977688], [6977641, 6977717], [6977641, 6977717], [6977658, 6977734], [6977661, 6977737], [6977671, 6977747], [6977672, 6977748], [6977672, 6977748], [6977675, 6977751], [6977676, 6977752], [6977681, 6977757], [6977689, 6977765], [6977692, 6977768], [6977701, 6977777], [6977702, 6977778], [6977706, 6977782], [6977708, 6977784], [6977711, 6977787], [6977711, 6977787], [6977716, 6977792], [6977723, 6977799], [6977726, 6977802], [6977727, 6977803], [6977736, 6977812], [6977738, 6977814], [6977743, 6977819], [6977745, 6977821], [6977750, 6977826], [6977753, 6977829], [6977754, 6977830], [6977754, 6977830], [6977755, 6977831], [6977760, 6977836], [6977772, 6977848], [6977779, 6977855], [6977782, 6977858], [6977795, 6977871], [6977811, 6977887], [6977823, 6977899], [6977828, 6977904], [6977836, 6977912], [6977838, 6977914], [6977839, 6977915], [6977846, 6977922], [6977849, 6977925], [6977863, 6977939], [6977870, 6977946], [6977871, 6977947], [6977871, 6977947], [6977874, 6977950], [6977875, 6977951], [6977878, 6977954], [6977893, 6977969], [6977902, 6977978], [6977903, 6977979], [6977911, 6977987], [6977927, 6978003], [6977929, 6978005], [6977941, 6978017], [6977943, 6978019], [6977945, 6978021], [6977951, 6978027], [6977955, 6978031], [6977960, 6978036], [6977965, 6978041], [6977969, 6978045], [6977981, 6978057], [6977983, 6978059], [6977983, 6978059], [6977989, 6978065], [6977991, 6978067], [6977994, 6978070], [6977996, 6978072], [6977999, 6978075], [6978000, 6978076], [6978000, 6978076], [6978007, 6978083], [6978034, 6978110], [6978034, 6978110], [6978039, 6978115], [6978045, 6978121], [6978047, 6978123], [6978053, 6978129], [6978054, 6978130], [6978054, 6978130], [6978056, 6978132], [6978065, 6978141], [6978066, 6978142], [6978072, 6978148], [6978081, 6978157], [6978093, 6978169], [6978096, 6978172], [6978107, 6978183], [6978112, 6978188], [6978114, 6978190], [6978120, 6978196], [6978121, 6978197], [6978121, 6978197], [6978129, 6978205], [6978129, 6978205], [6978133, 6978209], [6978138, 6978214], [6978142, 6978218], [6978144, 6978220], [6978145, 6978221], [6978157, 6978233], [6978157, 6978233], [6978161, 6978237], [6978162, 6978238], [6978169, 6978245], [6978171, 6978247], [6978171, 6978247], [6978188, 6978264], [6978192, 6978268], [6978193, 6978269], [6978194, 6978270], [6978197, 6978273], [6978198, 6978274], [6978198, 6978274], [6978198, 6978274], [6978198, 6978274], [6978198, 6978274], [6978198, 6978274], [6978198, 6978274], [6978198, 6978274], [6978198, 6978274], [6978198, 6978274], [6978198, 6978274], [6978200, 6978276], [6978223, 6978299], [6978224, 6978300], [6978224, 6978300], [6978228, 6978304], [6978229, 6978305], [6978236, 6978312], [6978244, 6978320], [6978247, 6978323], [6978252, 6978328], [6978255, 6978331], [6978258, 6978334], [6978287, 6978363], [6978303, 6978379], [6978314, 6978390], [6978328, 6978404], [6978328, 6978404], [6978341, 6978417], [6978349, 6978425], [6978354, 6978430], [6978358, 6978434], [6978361, 6978437], [6978376, 6978452], [6978376, 6978452], [6978382, 6978458], [6978384, 6978460], [6978389, 6978465], [6978392, 6978468], [6978400, 6978476], [6978405, 6978481], [6978417, 6978493], [6978420, 6978496], [6978420, 6978496], [6978422, 6978498], [6978422, 6978498], [6978423, 6978499], [6978430, 6978506], [6978434, 6978510], [6978443, 6978519], [6978444, 6978520], [6978447, 6978523], [6978448, 6978524], [6978451, 6978527], [6978454, 6978530], [6978459, 6978535], [6978492, 6978568], [6978531, 6978607], [6978535, 6978611], [6978539, 6978615], [6978541, 6978617], [6978546, 6978622]]
unpairedLens = {76:54}
pairedLens = {258: 1, 260: 1, 262: 1, 264: 2, 268: 2, 528: 1, 273: 1, 274: 1, 275: 1, 276: 2, 280: 1, 285: 3, 287: 1, 288: 1, 290: 1, 294: 1, 295: 1, 297: 2, 300: 2, 302: 1, 305: 1, 435: 1, 309: 1, 311: 3, 312: 1, 313: 1, 318: 2, 320: 1, 395: 1, 325: 1, 326: 1, 327: 4, 329: 1, 332: 1, 334: 1, 336: 1, 339: 1, 597: 1, 86: 1, 348: 1, 350: 1, 351: 1, 353: 1, 354: 1, 356: 1, 358: 1, 359: 1, 360: 1, 363: 1, 108: 2, 365: 2, 110: 2, 370: 2, 374: 1, 375: 1, 377: 1, 385: 1, 107: 1, 389: 1, 393: 1, 139: 1, 397: 1, 402: 2, 403: 1, 152: 1, 414: 1, 415: 1, 417: 1, 162: 1, 420: 1, 165: 2, 497: 1, 169: 1, 427: 1, 174: 1, 431: 1, 179: 1, 181: 1, 182: 1, 190: 1, 194: 1, 453: 1, 458: 1, 207: 1, 466: 1, 467: 1, 474: 1, 219: 1, 478: 1, 227: 1, 235: 1, 210: 1, 239: 2, 241: 1, 242: 1, 211: 1, 245: 1, 425: 1, 504: 1, 218: 1, 254: 1}

p = Pairs()
#startTime = time.time()
unpaired, paired = p.findPairsGreedy2(reads, unpairedLens, pairedLens)
#endTime = time.time()
#print('%0.3f s' % (endTime-startTime))

print('')
print(unpaired)
print(paired)