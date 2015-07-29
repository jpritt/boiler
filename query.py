#! /usr/bin/env python3
import sys
import argparse

import readSAM
import expand
import readPRO
import time
import re
import random

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
            print('Only .sam files are supported')

        header = ''
        while True:
            line = f.readline()
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
    expander = expand.Expander()
    pro = readPRO.ReadPRO(args['pro'])

    #intervals = [['2R', None, None],['2L', None, None],['3R', None, None],['3L', None, None],['4', None, None],['M', None, None],['X', None, None]]
    #intervals = [['2R', 100000, 200000]]

    lens = [1000, 10000, 100000, 1000000, 10000000, 20000000]
    #lens = [10000,100000]
    chroms = ['2R', '2L', '3R', '3L', 'X']
    chromLens = [21146708, 23011544, 27905053, 24543557, 22422827]

    '''
    print('Querying coverage')
    for l in lens:
        timeTrue = 0.0
        timePred = 0.0

        numIters = 3
        for c in range(len(chroms)):
            chrom = chroms[c]

            for _ in range(numIters):
                start = random.randint(0, chromLens[c]-l)

                startTime = time.time()
                trueCov = sam.getCoverage(chrom, start, start+l)
                endTime = time.time()
                timeTrue += float(endTime - startTime)

                startTime = time.time()
                predCov = expander.getCoverage(args['compressed'], chrom, start, start+l)
                endTime = time.time()
                timePred += float(endTime - startTime)

                if not len(trueCov) == len(predCov):
                    print('Error! Coverages not the same length!')
                    exit()

                for x in range(len(trueCov)):
                    if abs(trueCov[x] - predCov[x]) > 0.0001:
                        print('Error!')
                        print('%s (%d, %d)' % (chrom, start, start+l))
                        print(x)
                        for n in range(x-5,x+5):
                            print(str(trueCov[n]) + '\t' + str(predCov[n]))
                        exit()

        timeTrueAvg = timeTrue / float(numIters*len(chroms))
        timePredAvg = timePred / float(numIters*len(chroms))
        print('Length %d:' % l)
        print('  SAM:        %0.3f s' % timeTrueAvg)
        print('  Compressed: %0.3f s' % timePredAvg)
    '''

    print('Querying reads')
    print(sam.getReads('2R', 9779903, 9780903))
    print('')
    u,p = expander.getReads(args['compressed'], '2R', 9779903, 9780903)
    exit()
    for l in lens:
        timeTrue = 0.0
        timePred = 0.0

        numIters = 3
        for c in range(len(chroms)):
            chrom = chroms[c]

            for _ in range(numIters):
                start = random.randint(0, chromLens[c]-l)

                startTime = time.time()
                print('%s: %d - %d' % (chrom, start, start+l))
                trueReads = sam.getReads(chrom, start, start+l)
                endTime = time.time()
                timeTrue += float(endTime - startTime)

                startTime = time.time()
                predUnpaired, predPaired = expander.getReads(args['compressed'], chrom, start, start+l)
                endTime = time.time()
                timePred += float(endTime - startTime)

                print('%d true reads: %s' % (len(trueReads), str(trueReads)))
                #print('%d pred reads: %s' % (len(predReads), str(predReads)))
                print('')
            exit()


        timeTrueAvg = timeTrue / float(numIters*len(chroms))
        timePredAvg = timePred / float(numIters*len(chroms))
        print('Length %d:' % l)
        print('  SAM:        %0.3f s' % timeTrueAvg)
        print('  Compressed: %0.3f s' % timePredAvg)

        
def go_profile(args):
   pr = None
   if args['profile']:
       import cProfile
       import pstats
       import io
       pr = cProfile.Profile()
       pr.enable()
   go(args)
   if args['profile']:
       pr.disable()
       s = io.StringIO()
       sortby = 'tottime'
       ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
       ps.print_stats(30)
       print(s.getvalue())

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
