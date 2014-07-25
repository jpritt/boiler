#! /usr/bin/env python

import sys
import random


def findReadsInCoverage_v1(coverage, readLens, boundaries=None):
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

def findReadsInCoverage_v2(coverage, readLens, boundaries=None):
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
    while coverage[start] <= 0:
        start += 1
    end = len(coverage)
    while coverage[end-1] <= 0:
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


def findReadsInCoverage_BruteForce(coverage, start, end, readLens, level=0):
    #self.RLE(coverage, sys.stdout)

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


def findReadsInCoverage_v3(coverage, readLens, boundaries=None):
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

    reads = findReadsInCoverage_BruteForce(coverage, start, end, readLens)
    return reads


def findReads(readLens, lensLeft, lensRight, coverage, boundaries=None):
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


    #reads = self.findReadsInCoverage_v1(coverage, readLens, boundaries)
    reads = findReadsInCoverage_v2(coverage, readLens, boundaries)
    #reads = self.findReadsInCoverage_v3(coverage, readLens, boundaries)

    if reads == None:
        exit()

    # Pair left and right read lengths with long fragments
    # First sort reads and lengths by decreasing length
    reads.sort(key=lambda r: r[1]-r[0], reverse=True)
    lensLeftSorted = sorted(lensLeft.keys(), reverse=True)
    lensRightSorted = sorted(lensRight.keys(), reverse=True)

    pairedReads = []

    readId = 0

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







covLength = 500
coverage = [0]*covLength

readLen = 76
numReads = 25

numIters = 1000

for it in xrange(numIters):
    reads = []
    for i in xrange(numReads):
        start = random.randint(0, covLength-readLen-1)
        for j in xrange(start, start+readLen):
            coverage[j] += 1

    readLens = dict()
    readLens[readLen] = numReads

    foundReads, foundPaired = findReads(readLens, dict(), dict(), coverage)

    if len(foundPaired) > 0:
        print 'Error! %d paired reads found!' % len(foundPaired)
        print reads
        exit()

    if not len(foundReads) == numReads:
        print 'Error! %d reads found!' % len(foundReads)
        print reads
        exit()

    if sorted(reads) == sorted(foundReads):
        print 'Error!'
        print 'Originial reads: ' + str(sorted(reads))
        print 'Found reads: ' + str(sorted(foundReads))
        exit()

    print 'Correct! %d / %d done' % (it+1, numIters)