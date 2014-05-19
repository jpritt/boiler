#! /usr/bin/env python
import sys
import argparse
import genome

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

    print 'Parsing alignments'
    with open(args.alignments) as alignments:
        header = ''
        for line in alignments:
            if line[0] == '@':
                header += line
            else:
                break

        origGenome = genome.Genome(header)
        origGenome.parseAlignments(alignments)

    print 'Compressing'
    file_prefix = 'compressed/compressed'
    origGenome.compress(file_prefix)

    #compress.compressReads(alignments, chr_info, 'compressed')

    origGenome.writeSAM('minimal.sam')

    print 'Expanding'
    genomeExpanded = genome.Genome(header)
    genomeExpanded.expand(file_prefix)
    genomeExpanded.writeSAM('expanded.sam')

