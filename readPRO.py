#! /usr/bin/env python3
import sys

class ReadPRO:
    '''
        Parse the .pro file output by flux
    '''

    filename = None

    def __init__(self, filename):
        self.filename = filename

    def getGenes(self, chrom, start=None, end=None):
        if start == None:
            start = 0
        if end == None:
            end = sys.maxsize

        genes = set()
        with open(self.filename, 'r') as f:
            for line in f:
                row = line.rstrip().split('\t')
            
                if len(row) < 5:
                    continue

                gene = row[0].split(':')
                if gene[0] == chrom:
                    bounds = gene[1][:-1].split('-')

                    geneType = row[2]
                    if not geneType == 'CDS':
                        print(geneType)

                    genes.add((int(bounds[0]), int(bounds[1])))
        return sorted(list(genes))
