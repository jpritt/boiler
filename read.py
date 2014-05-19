#! /usr/bin/env python

class Read:
    def __init__(self, exons, xs=None):
        ''' exon is a list of (start,end) tuples marking each exonic region of this read.
            xs indicates the strand on which this gene lies, as described in the SAM file format.
        '''
        self.exons = exons
        self.xs = xs