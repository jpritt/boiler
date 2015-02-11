#! /usr/bin/env python

class Junction:
    ''' An exon junction spanned by at least 1 aligned read '''

    def __init__(self, exons, length):
        self.exons = exons
        self.label = ' '.join([str(e) for e in exons])
        self.length = length

        # Distribution of lengths of all reads crossing this junction
        self.readLens = dict()

        # Distribution of gap lengths for all paired-end reads crossing this junction
        #self.gaps = dict()

        # Distribution of left and right lengths for all paired-end reads crossing this junction
        self.lensLeft = dict()
        self.lensRight = dict()

        print length
        self.coverage = [0] * length
