#! /usr/bin/env python
import re

class ReadSAM:
    '''
        Support queries for SAM file
    '''

    filename = None
    chromosomes = None

    def __init__(self, filename, chromosomes):
        self.filename = filename
        self.chromosomes = chromosomes

    def getCoverage(self, chrom, start=None, end=None):
        if start == None or end == None:
            start = 0
            end = self.chromosomes[chrom]

        with open(self.filename, 'r') as f:
            coverage = [0] * (end-start)

            for line in f:
                row = line.rstrip().split('\t')

                if len(row) < 10:
                    continue

                if row[2] == chrom:
                    readStart = int(row[3])
                    cigar = row[5]

                    if readStart < end:
                        exons = self.parseCigar(cigar, readStart)

                        NH = 1.0
                        for s in row:
                            if s[:5] == 'NH:i:':
                                NH = 1.0 / float(s[5:])

                        for e in exons:
                            if e[0] < end and e[1] >= start:
                                for i in xrange(max(start,e[0]), min(end,e[1])):
                                    coverage[i-start] += NH

            return coverage


    def parseCigar(self, cigar, offset):
        ''' Parse the cigar string starting at the given index of the genome
            Returns a list of offsets for each exonic region of the read [(start1, end1), (start2, end2), ...]
        '''
        exons = []
        newExon = True

        # Parse cigar string
        match = re.search("\D", cigar)
        while match:
            index = match.start()
            length = int(''.join(cigar[:index]))

            if cigar[index] == 'N':
                # Separates contiguous exons, so set boolean to start a new one
                newExon = True
            elif cigar[index] == 'M':
                # If in the middle of a contiguous exon, append the length to it, otherwise start a new exon
                if newExon:
                    exons.append([offset, offset+length])
                    newExon = False
                else:
                    exons[-1][1] += length
            elif cigar[index] == 'D':
                # If in the middle of a contiguous exon, append the deleted length to it
                if not newExon:
                    exons[-1][1] += length

            offset += length
            cigar = cigar[index+1:]
            match = re.search("\D", cigar)

        return exons