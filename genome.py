#! /usr/bin/env python

import re
import chromosome
import read

class Genome:
    ''' A set of reads aligned to a genome '''

    def __init__(self, header):
        ''' Initialize GenomeAlignments object

            header: SAM file header containing (among other things) chromosome information
        ''' 
        chrom_info = []
        for line in header.split('\n'):
            if line[0:3] == '@SQ':
                chrom_info.append(line.strip().split('\t')[1:])
        self.chromosomes = self.init_chromosomes(chrom_info)

    def compress(self, file_prefix):
        ''' Compresses the alignments to multiple files

            file_prefix: Prefix for all output file names
        '''
        for k,v in self.chromosomes.items():
            v.compress(file_prefix + '.' + k)

    def expand(self, file_prefix):
        ''' Expand the alignments compressed in multiple files
        '''

        for k,v in self.chromosomes.items():
            v.expand(file_prefix + '.' + k)

    def init_chromosomes(self, chrom_info):
        '''
            Initialize genome chromosomes

            chrom_info: File handle containing a list of chromosome names and lengths
        '''
        chromosomes = dict()
        self.chr_names = []
        for line in chrom_info:
            chrom = line[0][3:]
            length = int(line[1][3:])
            chromosomes[chrom] = chromosome.Chromosome(chrom, length)
        return chromosomes

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
            if not row[2] in self.chromosomes.keys():
                print 'Chromosome ' + str(row[2]) + ' not found!'
                continue

            exons = self.parseCigar(row[5], int(row[3]))

            # find XS value:
            xs = None
            if len(exons) > 1:
                for r in row[11 : len(row)]:
                    if r[0:5] == 'XS:A:' or r[0:5] == 'XS:a:':
                        xs = r[5]

            pair = int(row[7])
            pair_offset = int(row[8])
            self.chromosomes[row[2]].processRead(read.Read(exons, xs), pair, pair_offset)

            '''
            if pair == 0 and pair_offset == 0:
                # unpairs read
                if len(exons) == 1:
                    self.chromosomes[row[2]].addUnspliced(read.Read(exons))
                else:
                    # find XS value:
                    for r in row[11 : len(row)]:
                        if r[0:5] == 'XS:A:' or r[0:5] == 'XS:a:':
                            xs = r[5]

                    self.chromosomes[row[2]].addSpliced(read.Read(exons, xs))

            else:
                if pair in unpaired:
                    # create new read

                    # Tophat alignments should be in increasing order of index
                    if unpaired[pair].exons[0][0] > exons[0][0]:
                        print 'Error! Tophat file not sorted!'

                    paired_exons = 
                    paired_read = read.Read()
            '''

        # finalize exon list in each chromosome
        for c in self.chromosomes.values():
            c.finalizeExons()
            c.finalizeReads()

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

    '''
    def testExpand(self):
        print 'expanding'
        cov = [1,2,2,2,1,1,2,3,4,4,4,2]
        reads = dict()
        reads[3] = 1
        reads[4] = 1
        reads[5] = 3
        reads[6] = 1

        print self.chromosomes[self.chromosomes.keys()[0]].findReads(reads, cov)
        exit()
    '''

    def writeSAM(self, filename):
        ''' Write all stored alignments to a new SAM file
        '''

        with open(filename, 'w') as f:
            # write header
            f.write('@HD\tVN:1.0\tSO:unsorted\n')
            for c in self.chromosomes.values():
                f.write('@SQ\tSN:' + c.name + '\tLN:' + str(c.length) + '\n')

            # write alignments
            for c in self.chromosomes.values():
                c.writeSAM(f)
    