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

        self.coverage = [0] * length

    def differenceEncode(self):
        ''' Difference run-length encode this junctions coverage vector
        '''

        self.covDRLE = []

        oldVal = 0
        val = self.coverage[0]
        length = 0
        for v in self.coverage:
            if v == val:
                length += 1
            else:
                self.covDRLE.append([val - oldVal, length])
                oldVal = val
                val = v
                length = 1
        self.covDRLE.append([val - oldVal, length])
