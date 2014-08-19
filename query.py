#! /usr/bin/env python
import sys
import argparse
import genome
import compressedAlignment
import time

def getTrueCoverage(sam_file, chrom, start=None, end=None):
    with open(sam_file, 'r') as f:
        alignments = f.read().split('\n')
        header = ''
        for line in alignments:
            if line[0] == '@':
                header += line + '\n'
            else:
                break

        chromosomes = dict()
        for line in header.split('\n'):
            if line[0:3] == '@SQ':
                row = line.strip().split('\t')
                chromosomes[row[1][3:]] = int(row[2][3:])

        if start == None:
            start = 0
        if end == None:
            end = len(chromosomes[chrom])

        coverage = [0] * (end-start)

        for row in alignments:
            read = row.rstrip().split()
            if read[2] == chrom:
                

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--alignments', type=str, required=True, 
        help='Full path of SAM file containing aligned reads')
    
    args = parser.parse_args(sys.argv[1:])

    alignments = args.alignments
    form = alignments[-3:len(alignments)]
    if (form != 'sam'):
        print 'Only .sam files are supported'

    sys.stdout.write('Parsing alignments... ')
    start = time.time()
    with open(args.alignments) as f:
        alignments = f.read().split('\n')
        header = ''
        for line in alignments:
            if line[0] == '@':
                header += line + '\n'
            else:
                break

        origGenome = genome.Genome(header)
        origGenome.parseAlignments(alignments)
    end = time.time()
    parseTime = end-start
    print '%fs' % parseTime


    sys.stdout.write('Compressing... ')
    start = time.time()
    file_prefix = 'compressed/compressed'
    origGenome.compress(file_prefix)
    end = time.time()
    compressTime = end-start
    print '%fs' % compressTime

    sys.stdout.write('Querying...')
    compressed = compressedAlignment.CompressedAlignment(header, 'compressed/compressed.spliced.txt', 'compressed/compressed.unspliced.txt')

    