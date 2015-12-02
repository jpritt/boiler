#! /usr/bin/env python

class CrossBundleBucket:
    ''' A bucket spanning two bundles '''

    def __init__(self, bundleA, exonIdsA, bundleB, exonIdsB):
        self.bundleA = bundleA
        self.exonIdsA = exonIdsA
        self.bundleB = bundleB
        self.exonIdsB = exonIdsB

        # Distribution of lengths of all paired reads in this bucket
        self.pairedLens = dict()

        # Distribution of left and right lengths for all paired-end reads in this bucket
        self.lensLeft = dict()
        self.lensRight = dict()

    def set_length(self, length):
        self.length = length
        self.coverage = [[0, self.length]]

    def add_pair(self, readA, readB):
        length = self.length - readA.startOffset - readB.endOffset

        # update bucket coverage vector in dictionary
        self.updateCov(readA.startOffset, readA.length)
        self.updateCov(self.length-readB.endOffset-readB.length, readB.length)

        #print(self.length)
        #print('%d, %d' % (readA.startOffset, readA.endOffset))
        #print('%d, %d' % (readB.startOffset, readB.endOffset))
        #print('')

        # update lengths
        if readA.length in self.lensLeft:
            self.lensLeft[readA.length] += 1
        else:
            self.lensLeft[readA.length] = 1

        if readB.length in self.lensRight:
            self.lensRight[readB.length] += 1
        else:
            self.lensRight[readB.length] = 1

        if length in self.pairedLens:
            self.pairedLens[length] += 1
        else:
            self.pairedLens[length] = 1

    def updateCov(self, start, length):
        i = 0

        while start >= self.coverage[i][1]:
            start -= self.coverage[i][1]
            i += 1

        if start > 0:
            self.coverage = self.coverage[:i] + [ [self.coverage[i][0], start], [self.coverage[i][0], self.coverage[i][1]-start] ] + self.coverage[i+1:]
            i += 1

        rle_len = len(self.coverage)

        while i < rle_len and length >= self.coverage[i][1]:
            self.coverage[i][0] += 1
            if self.coverage[i][0] < 0:
                self.coverage[i][0] = 0

            length -= self.coverage[i][1]
            i += 1

        if i < rle_len and length > 0:
            self.coverage = self.coverage[:i] + [ [max(self.coverage[i][0]+1,0), length], [self.coverage[i][0], self.coverage[i][1]-length] ] + self.coverage[i+1:]
