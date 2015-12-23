#! /usr/bin/env python

class Bucket:
    ''' A bucket spanned by at least 1 aligned read '''

    def __init__(self, exons, length, boundaries=None):
        self.exons = exons
        self.label = ' '.join([str(e) for e in exons])
        self.length = length

        # Offsets of boundaries between subexons
        self.boundaries = boundaries

        #self.pairs = []

        #self.paired = []
        #self.unpaired = []

        # Distribution of lengths of all paired and unpaired reads crossing this junction
        self.pairedLens = dict()
        self.unpairedLens = dict()

        # Distribution of left and right lengths for all paired-end reads crossing this junction
        self.lensLeft = dict()
        self.lensRight = dict()

        #self.coverage = [[0, length]]
        self.coverage = [0] * length

    def add_paired(self, p):
        p.length = self.length - p.startOffset - p.endOffset

        # update bucket coverage vector in dictionary
        self.updateCov(p.startOffset, p.lenLeft)
        self.updateCov(self.length-p.endOffset-p.lenRight, p.lenRight)

        # update lengths
        if p.lenLeft in self.lensLeft:
            self.lensLeft[p.lenLeft] += 1
        else:
            self.lensLeft[p.lenLeft] = 1

        if p.lenRight in self.lensRight:
            self.lensRight[p.lenRight] += 1
        else:
            self.lensRight[p.lenRight] = 1

        if p.length in self.pairedLens:
            self.pairedLens[p.length] += 1
        else:
            self.pairedLens[p.length] = 1

    def add_unpaired(self, r):
        # update bucket coverage vector in dictionary
        self.updateCov(r.startOffset, r.length)

        # update unpairedLens
        if r.length in self.unpairedLens:
            self.unpairedLens[r.length] += 1
        else:
            self.unpairedLens[r.length] = 1

    '''
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
    '''

    def updateCov(self, start, length):
        for i in range(start, start+length):
            self.coverage[i] += 1

