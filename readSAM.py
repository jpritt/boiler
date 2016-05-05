#! /usr/bin/env python3
import re

class ReadSAM:
    '''
        Support queries for SAM file
    '''

    filename = None
    chromosomes = None

    def __init__(self, filename, chromosomes):
        '''
            Assumes SAM file is sorted
        '''
        self.filename = filename
        self.chromosomes = chromosomes

    def getCoverage(self, chrom, start=None, end=None):
        if start == None:
            start = 0
        else:
            start = max(start, 0)
        if end == None:
            end = self.chromosomes[chrom]
        else:
            end = min(end, self.chromosomes[chrom])

        with open(self.filename, 'r') as f:
            coverage = [0.0] * (end-start)

            for line in f:
                row = line.rstrip().split('\t')

                if len(row) < 10:
                    continue

                if row[2] == chrom:
                    readStart = int(row[3])
                    cigar = row[5]

                    if readStart < end:
                        exons = self.parseCigar(cigar, readStart)

                        if exons[-1][1] > start:
                            NH = 1.0
                            for s in row:
                                if s[:5] == 'NH:i:':
                                    NH = 1.0 / float(s[5:])

                            for e in exons:
                                if e[0] < end and e[1] >= start:
                                    for i in range(max(start,e[0]), min(end,e[1])):
                                        coverage[i-start] += NH
                    else:
                        return coverage

            return coverage

    def getReads(self, chrom, start=None, end=None):
        '''
        Assumes the SAM file is sorted
        :param chrom:
        :param start:
        :param end:
        :return:
        '''
        if start == None or end == None:
            start = 0
            end = self.chromosomes[chrom]

        offset = 0
        for k in sorted(self.chromosomes.keys()):
                if not k == chrom:
                    offset += self.chromosomes[k]
                    start += self.chromosomes[k]
                    end += self.chromosomes[k]
                else:
                    break

        unmatched = dict()

        unpaired = []
        paired = []
        sawChrom = False # Set this flag to true after we've seen the first read from chrom
        with open(self.filename, 'r') as f:
            for line in f:
                row = line.rstrip().split('\t')

                if len(row) < 10:
                    continue

                if row[2] == chrom:
                    sawChrom = True
                    readStart = int(row[3])
                    cigar = row[5]

                    if row[6] == '=':
                        pairStart = int(row[7])

                        fragmentStart = min(readStart, pairStart) + offset
                        fragmentEnd = fragmentStart + abs(int(row[8])) + offset
                        if fragmentStart >= end or fragmentEnd <= start:
                            continue

                        #if readStart >= end and pairStart >= end:
                        #    continue

                        NH = 1
                        for r in row[11 : len(row)]:
                            if r[0:5] == 'XS:A:' or r[0:5] == 'XS:a:':
                                xs = r[5]
                            elif r[0:3] == 'NH:':
                                NH = int(r[5:])


                        if row[0] in unmatched:
                            found = False
                            m = unmatched[row[0]]
                            for i in range(len(m)):
                                if m[i][0] == pairStart and m[i][2] == readStart and m[i][3] == NH:
                                    r1 = self.parseCigar(m[i][1], m[i][0]+offset)
                                    r2 = self.parseCigar(cigar, readStart+offset)

                                    if self.readOverlapsRegion(r1, start, end) or self.readOverlapsRegion(r2, start, end):
                                        paired.append([r1, r2])

                                    del unmatched[row[0]][i]
                                    break

                            if not found:
                                unmatched[row[0]].append((readStart, cigar, pairStart, NH))
                        else:
                            unmatched[row[0]] = [(readStart, cigar, pairStart, NH)]
                    else:
                        exons = self.parseCigar(cigar, readStart+offset)

                        if readStart < end and exons[-1][1] > start:
                            if self.readOverlapsRegion(exons, start, end):
                                NH = 1
                                for r in row[11 : len(row)]:
                                    if r[0:5] == 'XS:A:' or r[0:5] == 'XS:a:':
                                        xs = r[5]
                                    elif r[0:3] == 'NH:':
                                        NH = int(r[5:])
                                unpaired.append(exons)
                else:
                    if sawChrom:
                        return unpaired, paired

        return unpaired, paired

    def readOverlapsRegion(self, exons, start, end):
        '''

        :param exons: List of exon bounds
        :param start: Start of region of interest
        :param end: End of region of interest
        :return: True if any of the exons overlap the region of interest
        '''

        for e in exons:
            if e[0] >= end:
                # Passed region
                return False
            elif e[1] > start:
                return True
        return False

    def getGenes(self, chrom, start=None, end=None, overlapRadius=50):
        if start == None or end == None:
            start = 0
            end = self.chromosomes[chrom]

        genes = []
        with open(self.filename, 'r') as f:
            for line in f:
                row = line.rstrip().split('\t')

                if len(row) < 10:
                    continue

                if row[2] == chrom:
                    readStart = int(row[3])

                    if readStart >= end:
                        break

                    if row[6] == '=' and int(row[7]) > readStart:
                        genes.append((readStart, int(row[7])))
                    else:
                        cigar = row[5]
                        if readStart >= start and readStart < end:
                            exons = self.parseCigar(cigar, readStart)

                            if exons[-1][1] <= end:
                                genes.append((exons[0][0], exons[-1][1]))


            genes.sort()

            i = 0
            while i < (len(genes)-1):
                while i < (len(genes)-1) and genes[i+1][0] - genes[i][1] <= overlapRadius:
                    genes[i] = (genes[i][0], max(genes[i][1], genes[i+1][1]))
                    del genes[i+1]
                i += 1

        return genes


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

            if not cigar[index] == 'S':
                offset += length
            cigar = cigar[index+1:]
            match = re.search("\D", cigar)

        return exons
