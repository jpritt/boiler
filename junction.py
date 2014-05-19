#! /usr/bin/env python

class Junction:
    ''' An exon junction spanned by at least 1 aligned read '''

    def __init__(self, exons, length):
        self.exons = exons
        self.label = ' '.join([str(e) for e in exons])
        self.length = length

        self.readLens = dict()
        self.coverage = [0] * length
