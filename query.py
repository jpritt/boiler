#! /usr/bin/env python
import sys
import argparse

import readSAM
import readCompressed
import readPRO
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

def go(args):
    alignmentsFile = args['alignments']
    with open(alignmentsFile) as f:
        form = args['alignments'][-3:len(alignmentsFile)]
        if (form != 'sam'):
            print 'Only .sam files are supported'

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


        sam = readSAM.ReadSAM(alignmentsFile, chromosomes)
        compressed = readCompressed.ReadCompressed(args['compressed'] + '.unspliced.txt', args['compressed'] + '.spliced.txt', chromosomes)
        pro = readPRO.ReadPRO(args['pro'])

        #intervals = [['2R', None, None],['2L', None, None],['3R', None, None],['3L', None, None],['4', None, None],['M', None, None],['X', None, None]]
        intervals = [['2L', None, None]]
        
        
        print 'Querying coverage...'
        for i in intervals:
            chrom = i[0]
            start = i[1]
            end = i[2]

            if start == None or end == None:
                print '%s: (%d, %d)' % (chrom, 0, chromosomes[chrom])
            else:
                print '%s: (%d, %d)' % (chrom, start, end)
            if start == None or end == None:
                length = chromosomes[chrom]
            else:
                length = end-start
            startTime = time.time()
            trueCov = sam.getCoverage(chrom, start, end)
            endTime = time.time()
            trueTime = float(endTime - startTime)
            print '%fs, %d bases (%f bases/s)' % (trueTime, length, float(length)/trueTime)

            startTime = time.time()
            predCov = compressed.getCoverage(chrom, start, end)
            endTime = time.time()

            correct = 0
            for x in xrange(len(trueCov)):
                if abs(trueCov[x] - predCov[x]) < 0.0001:
                    correct += 1
                #else:
                #    for y in xrange(x-10, x+40):
                #        print '%0.2f\t%0.2f' % (trueCov[y], predCov[y])
                #    exit()
            predTime = float(endTime - startTime)
            print '%fs, %d bases (%f bases/s)' % (predTime, length, float(length)/predTime)
            print '%0.3f correct' % (float(correct) / float(len(trueCov)))
            print ''
        
        '''
        print 'Querying genes...'
        for i in intervals:
            chrom = i[0]
            start = i[1]
            end = i[2]

            if start == None or end == None:
                print '%s: (%d, %d)' % (chrom, 0, chromosomes[chrom])
            else:
                print '%s: (%d, %d)' % (chrom, start, end)

            trueGenes = pro.getGenes(chrom, start, end)

            startTime = time.time()
            predGenes = compressed.getGenes(chrom, start, end)
            endTime = time.time()

            correct = 0
            for i in trueGenes:
                if i in predGenes:
                    correct += 1

            predTime = float(endTime - startTime)
            print '%0.3fs' % predTime
            print '%d / %d of true genes were predicted' % (correct, len(trueGenes))
            print '%d / %d of predicted genes were true' % (correct, len(predGenes))

            print trueGenes[:10]
            print predGenes[:10]
            print ''
        '''
def go_profile(args):
   pr = None
   if args['profile']:
       import cProfile
       import pstats
       import StringIO
       pr = cProfile.Profile()
       pr.enable()
   go(args)
   if args['profile']:
       pr.disable()
       s = StringIO.StringIO()
       sortby = 'tottime'
       ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
       ps.print_stats(30)
       print s.getvalue()

if __name__ == '__main__':
    # Print file's docstring if -h is invokedc
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--alignments', type=str, required=True, 
        help='Full path of SAM file containing aligned reads')
    parser.add_argument('--compressed', type=str, required=True, 
        help='Full path of directory containing compressed reads')
    parser.add_argument('--pro', type=str, required=True, 
        help='Full path of flux .pro output file')
    parser.add_argument("--profile", help="Run speed profiling",
        action="store_true")
    
    args = parser.parse_args(sys.argv[1:])
    go_profile(vars(args))
