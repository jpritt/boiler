import time
import bisect
import random
import copy

class Reads:
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
        return reads

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

        '''
        if boundaries == None:
            boundaries = [0, len(coverage)]
        else:
            boundaries = list(boundaries)
        if not len(boundaries) == 2:
            print('%d boundaries' % len(boundaries))
        boundBottom = 0
        boundTop = len(boundaries)-1
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

            print(coverage)
            print(readLens)

            closestEndpoint = None
            for length in range(1, maxLen+1):
                if (readStart+length == end) or (readStart+length < end and coverage[readStart + length] < coverage[readStart + length - 1]):
                    if length in readLens:
                        readEnd = readStart + length
                        reads.append([readStart, readEnd])

                        print('Adding read (%d, %d)' % (readStart, readEnd))

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
                print('No perfect length on step down')
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

            print('')
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

                print(coverage)
                print(readLens)

                closestEndpoint = None
                for length in range(1, maxLen+1):
                #for length in range(1, maxLen+1):
                    if (readEnd-length == start) or (readEnd-length > start and coverage[readEnd - length] > coverage[readEnd - length - 1]):
                        if length in readLens:
                            readStart = readEnd - length
                            reads.append([readStart, readEnd])

                            print('Adding read (%d, %d)' % (readStart, readEnd))

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
                    print('No perfect read on step up')
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
            print('')
        return reads

unpairedLens = {76: 3}
pairedLens = {74: 1, 62: 1}
lensLeft = {74: 1, 62: 1}
lensRight = {74: 1, 62: 1}
cov = [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3]

r = Reads()
reads = r.findReads(unpairedLens, pairedLens, lensLeft, lensRight, cov)
print(reads)