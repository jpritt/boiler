#! /usr/bin/env python

import sys
import random
import alignments

chroms = {'2L': 10000}
aligned = alignments.Alignments(chroms)



covLength = 500
coverage = [0]*covLength

readLen = 76
numReads = 25

numIters = 1

for it in range(numIters):
    reads = []
    for i in range(numReads):
        start = random.randint(0, covLength-readLen-1)
        l = random.randint(-3, 3)
        for j in range(start, start+readLen+l):
            coverage[j] += 1
        reads.append([start, start+readLen+l])

    readLens = dict()
    readLens[readLen] = numReads



    foundReads = aligned.findReadsInCoverage_v3(coverage, readLens)

    print('Originial reads: ' + str(sorted(reads)))
    print('Found reads:     ' + str(sorted(foundReads)))
    print(sorted(reads) == sorted(foundReads))