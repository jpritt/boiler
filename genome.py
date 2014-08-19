#! /usr/bin/env python

import re
#import chromosome
import read
import alignments

class Genome:
    ''' A set of reads aligned to a genome '''

    def __init__(self, header):
        ''' Initialize GenomeAlignments object

            header: SAM file header containing (among other things) chromosome information
        ''' 

        '''        
        chrom_info = []
        for line in header.split('\n'):
            if line[0:3] == '@SQ':
                chrom_info.append(line.strip().split('\t')[1:])
        '''

        self.chromosomes = dict()
        for line in header.split('\n'):
            if line[0:3] == '@SQ':
                row = line.strip().split('\t')
                self.chromosomes[row[1][3:]] = int(row[2][3:])

        self.alignments = alignments.Alignments(self.chromosomes)
        #self.chromosomes = self.init_chromosomes(chrom_info)

    def compress(self, file_prefix):
        ''' Compresses the alignments to multiple files

            file_prefix: Prefix for all output file names
        '''
        self.alignments.compress(file_prefix)

    def expand(self, file_prefix):
        ''' Expand the alignments compressed in multiple files
        '''

        self.alignments.expand(file_prefix)

    # def init_chromosomes(self, chrom_info):
    #     '''
    #         Initialize genome chromosomes

    #         chrom_info: File handle containing a list of chromosome names and lengths
    #     '''
    #     chromosomes = dict()
    #     self.chr_names = []
    #     for line in chrom_info:
    #         chrom = line[0][3:]
    #         length = int(line[1][3:])
    #         chromosomes[chrom] = chromosome.Chromosome(chrom, length)
    #     return chromosomes

    def parseAlignments(self, alignments):
        ''' Parse a file in SAM format
            Add the exonic region indices for each read to the correct ChromosomeAlignments object
            
            alignments: SAM filehandle containing aligned reads
        '''

        # paired reads that have not yet found their partner
        unpaired = dict()

        for line in alignments:
            row = line.strip().split('\t')
            if len(row) < 5:
                continue

            chromosome = str(row[2])
            
            if not row[2] in self.chromosomes.keys():
                print 'Chromosome ' + str(row[2]) + ' not found!'
                continue

            exons = self.parseCigar(row[5], int(row[3]))

            # find XS value:
            xs = None
            NH = 1
            for r in row[11 : len(row)]:
                if r[0:5] == 'XS:A:' or r[0:5] == 'XS:a:':
                    xs = r[5]
                elif r[0:3] == 'NH:':
                    NH = int(r[5:])

            if not row[6] == '*':
                if row[6] == '=':
                    pair_chrom = chromosome
                else:
                    pair_chrom = row[6]
                pair_index = int(row[7])
                self.alignments.processRead(read.Read(chromosome, exons, xs, NH), pair_chrom, pair_index)
            else:
                self.alignments.processRead(read.Read(chromosome, exons, xs, NH))

        self.alignments.finalizeExons()
        self.alignments.finalizeReads()

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


    def writeSAM(self, filename):
        ''' Write all stored alignments to a new SAM file
        '''

        with open(filename, 'w') as f:
            # write header
            f.write('@HD\tVN:1.0\tSO:unsorted\n')
            for k,v in self.chromosomes.items():
                f.write('@SQ\tSN:' + str(k) + '\tLN:' + str(v) + '\n')

            # write alignments
            self.alignments.writeSAM(f)
    