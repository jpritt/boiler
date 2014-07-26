#! /usr/bin/env python
import sys
import argparse
import genome
import time

def read_chromosome_info(chrom_info):
    ''' Read chromosome names and lengths from the given file
    '''
    chromosome_info = []
    for line in chrom_info:
        row = line.rstrip().split('\t')
        chromosome_info.append([row[0], int(row[1])])
    return chromosome_info

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


    #compress.compressReads(alignments, chr_info, 'compressed')

    #origGenome.writeSAM('minimal.sam')

    sys.stdout.write('Expanding... ')

    genomeExpanded = genome.Genome(header)
    start = time.time()
    genomeExpanded.expand(file_prefix)
    end = time.time()
    expandTime = end-start
    print '%fs' % expandTime
    
    genomeExpanded.writeSAM('expanded.sam')

