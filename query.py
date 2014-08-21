#! /usr/bin/env python
import sys
import argparse
import genome
import compressedAlignment
import time
import re

def parseCigar(cigar, offset):
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

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--alignments', type=str, required=True, 
        help='Full path of SAM file containing aligned reads')
    parser.add_argument('--compressed', type=str, required=False, 
        help='Full path of directory containing compressed reads')
    
    args = parser.parse_args(sys.argv[1:])

    with open(args.alignments) as f:
        form = args.alignments[-3:len(args.alignments)]
        if (form != 'sam'):
            print 'Only .sam files are supported'

        sys.stdout.write('Parsing alignments... ')
        start = time.time()

        alignments = f.read().split('\n')
        header = ''
        for line in alignments:
            if line[0] == '@':
                header += line + '\n'
            else:
                break

        if args.compressed == None:
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
        else:
            file_prefix = args.compressed


        print 'Calculating coverage...'

        # Parse chromosome lengths
        chromosomes = dict()
        for line in header.split('\n'):
            if line[0:3] == '@SQ':
                row = line.strip().split('\t')
                chromosomes[row[1][3:]] = int(row[2][3:])

        # Offset of each chromosome from the start of the genome
        chromOffsets = dict()
        nextOffset = 0
        for c in chromosomes.keys():
            chromOffsets[c] = nextOffset
            nextOffset += chromosomes[c]

        genomeLen = nextOffset
        coverage = [0] * genomeLen

        with open(args.alignments) as f:
            for line in f:
                row = line.rstrip().split('\t')
                if len(row) < 11:
                    continue

                chrom = row[2]
                start = int(row[3])
                cigar = row[5]
                exons = parseCigar(cigar, start + chromOffsets[chrom])

                for e in exons:
                    for i in xrange(e[0], e[1]):
                        coverage[i] += 1

        sys.stdout.write('Querying...')
        compressed = compressedAlignment.CompressedAlignment(header, file_prefix + '.spliced.txt', file_prefix + '.unspliced.txt')

        intervals = [['2R', 1000, 1050]]

        for i in intervals:
            chrom = i[0]
            start = i[1]
            end = i[2]

            trueCov = coverage[start+chromOffsets[chrom] : end+chromOffsets[chrom]]
            predCov = compressed.getCoverage(chrom, start, end)

            for x in xrange(len(trueCov)):
                print trueCov[x] + '\t' + predCov[x]

        