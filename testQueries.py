#! /usr/bin/env python3
import sys
import argparse

import readSAM
import expand
import time
import re
import math
import random
import os

def queryCoverageRanges(filename, sam, expander):
    lens = [1000, 10000, 100000, 1000000, 10000000, 20000000]
    chroms = ['2R', '2L', '3R', '3L', 'X']
    chromLens = [21146708, 23011544, 27905053, 24543557, 22422827]

    for l in lens:
        timeTrue = 0.0
        timePred = 0.0

        numIters = 5
        for c in range(len(chroms)):
            chrom = chroms[c]

            for _ in range(numIters):
                start = random.randint(0, chromLens[c]-l)

                startTime = time.time()
                trueCov = sam.getCoverage(chrom, start, start+l)
                endTime = time.time()
                timeTrue += float(endTime - startTime)

                startTime = time.time()
                predCov = expander.getCoverage(filename, chrom, start, start+l)
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

def queryCoverageInBundles(filename, expander, chrom, bamFile, bedtoolsPath, chromsFile, limit, sam=None):
    if chrom == 'all':
        print('Querying coverage for all chromsomes')
        chroms = expander.getChromosomes(filename)
    else:
        print('Querying coverage for chromsome ' + chrom)
        chroms = [chrom]

    lens = []
    timesTrue = []
    timesPred = []
    for c in chroms:
        #print('Getting bundle coverages for chromsome %s' % c)
        bundles = expander.getGeneBounds(filename, c)
        num_bundles = len(bundles)

        if limit and limit < num_bundles:
            r = num_bundles / limit
            for i in range(limit):
                start = int(i * r)
                end = int((i+1) * r) - 1
                b = bundles[random.randint(start,end)]

                lens.append(b[1]-b[0])

                startTime = time.time()
                #trueCov = sam.getCoverage(c, b[0], b[1])
                os.system("samtools view -b -h " + bamFile + " " + c + ":" + str(b[0]) + "-" + str(b[1]) + " | " + bedtoolsPath + "/bin/genomeCoverageBed -bga -split -ibam stdin -g " + chromsFile + " > coverage.txt")
                endTime = time.time()
                timesTrue.append(endTime - startTime)

                startTime = time.time()
                predCov = expander.getCoverage(filename, c, b[0], b[1])
                endTime = time.time()
                timesPred.append(endTime - startTime)

        else:
            for b in bundles:
                lens.append(b[1]-b[0])

                startTime = time.time()
                #trueCov = sam.getCoverage(c, b[0], b[1])
                os.system("samtools view -b -h " + bamFile + " " + c + ":" + str(b[0]) + "-" + str(b[1]) + " | " + bedtoolsPath + "/bin/genomeCoverageBed -bga -split -ibam stdin -g " + chromsFile + " > coverage.txt")
                endTime = time.time()
                timesTrue.append(endTime - startTime)

                startTime = time.time()
                predCov = expander.getCoverage(filename, c, b[0], b[1])
                endTime = time.time()
                timesPred.append(endTime - startTime)

    print('Average SAM query time:    %fs' % (sum(timesTrue) / len(timesTrue)))
    print('Average Boiler query time: %fs' % (sum(timesPred) / len(timesPred)))

    return lens, timesTrue, timesPred

def queryReadsInBundles(filename, expander, chrom, bamFile, bedtoolsPath, chromsFile, limit):
    if chrom == 'all':
        print('Querying reads for all chromsomes')
        chroms = expander.getChromosomes(filename)
    else:
        print('Querying reads for chromsome ' + chrom)
        chroms = [chrom]

    lens = []
    timesTrue = []
    timesPred = []
    for c in chroms:
        #print('Getting bundle coverages for chromsome %s' % c)
        bundles = expander.getGeneBounds(filename, c)
        num_bundles = len(bundles)

        if limit and limit < num_bundles:
            r = num_bundles / limit
            for i in range(limit):
                start = int(i * r)
                end = int((i+1) * r) - 1
                b = bundles[random.randint(start,end)]

                lens.append(b[1]-b[0])

                startTime = time.time()
                os.system("samtools view -h -o reads.sam " + bamFile + " " + c + ":" + str(b[0]) + "-" + str(b[1]))
                endTime = time.time()
                timesTrue.append(endTime - startTime)

                startTime = time.time()
                _, predUnpaired, predPaired = expander.getReads(filename, c, b[0], b[1])
                endTime = time.time()
                timesPred.append(endTime - startTime)
        else:
            for b in bundles:
                lens.append(b[1]-b[0])

                startTime = time.time()
                os.system("samtools view -h -o reads.sam " + bamFile + " " + c + ":" + str(b[0]) + "-" + str(b[1]))
                endTime = time.time()
                timesTrue.append(endTime - startTime)

                startTime = time.time()
                _, predUnpaired, predPaired = expander.getReads(filename, c, b[0], b[1])
                endTime = time.time()
                timesPred.append(endTime - startTime)


    print('Average SAM query time:    %fs' % (sum(timesTrue) / len(timesTrue)))
    print('Average Boiler query time: %fs' % (sum(timesPred) / len(timesPred)))

    return lens, timesTrue, timesPred

def go(args, mode):
    form = args['alignments'][-3:len(args['alignments'])]
    if form == 'sam':
        samFile = args['alignments']
        bamFile = args['alignments'][:-3] + 'bam'
    elif form == 'bam':
        bamFile = args['alignments']
        samFile = args['alignments'][:-3] + 'sam'
    else:
        print('Only .sam and .bam files are supported')
        exit()

    with open(samFile) as f:
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

    with open(args['chroms'], 'w') as f:
        for k,v in chromosomes.items():
            f.write(k + '\t' + str(v) + '\n')

    #sam = readSAM.ReadSAM(samFile, chromosomes)
    expander = expand.Expander()
    #pro = readPRO.ReadPRO(args['pro'])


    if mode == 0:
        lens, timesTrue, timesPred = queryCoverageInBundles(args['compressed'], expander, args['chrom'], bamFile, args['bedtools_path'], args['chroms'], args['limit'])
    else:
        lens, timesTrue, timesPred = queryReadsInBundles(args['compressed'], expander, args['chrom'], bamFile, args['bedtools_path'], args['chroms'], args['limit'])

    if args['output']:
        with open(args['output'], 'w') as f:
            for i in range(len(lens)):
                f.write(str(lens[i]) + '\t' + str(timesTrue[i]) + '\t' + str(timesPred[i]) + '\n')

    if args['plot']:
        import matplotlib.pyplot as plt

        if mode == 0:
            s = 'cov'
        elif mode == 0:
            s = 'reads'

        plt.scatter(lens, timesTrue)
        plt.xlabel('Bundle Length')
        plt.ylabel('SAM Query Time (s)')
        plt.savefig('sam_query_time_' + s + '.png')
        plt.clf()

        plt.scatter(lens, timesPred)
        plt.xlabel('Bundle Length')
        plt.ylabel('Boiler Query Time (s)')
        plt.savefig('boiler_query_time_' + s + '.png')
        plt.clf()

        plt.scatter(timesTrue, timesPred)
        plt.xlabel('SAM Query Time (s)')
        plt.ylabel('Boiler Query Time (s)')
        plt.savefig('boiler_sam_query_time_' + s + '.png')
        plt.clf()


    os.remove(args['chroms'])
        
def go_profile(args, mode):
   pr = None
   if args['profile']:
       import cProfile
       import pstats
       import io
       pr = cProfile.Profile()
       pr.enable()
   go(args, mode)
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
        help='Full path of SAM or BAM file containing aligned reads')
    parser.add_argument('--compressed', type=str, required=True, 
        help='Full path of directory containing compressed reads')
    #parser.add_argument('--pro', type=str, required=True, help='Full path of flux .pro output file')
    parser.add_argument('--chrom', type=str, required=True, help="Chromsome to query. Write 'all' to query bundles on all chromosomes")
    parser.add_argument('--chroms', type=str, help="Temporary file to write chromosomes to. Default = chroms.genome")
    parser.add_argument('--bedtools-path', type=str, required=True, help="Path to bedtools main directory")
    parser.add_argument("--profile", help="Run speed profiling",
        action="store_true")
    parser.add_argument("--output", type=str, help="File to store timing results")
    parser.add_argument("--plot", help="Plot timing comparison graphs", action="store_true")
    parser.add_argument("--timings", help="If present, read timings file and simply plot results", action="store_true")
    parser.add_argument("--mode", type=str, required=True, help="Either 'cov' or 'reads'")
    parser.add_argument("--limit", type=int, help="Maximum number of bundles from each chromosome to query")
    
    args = parser.parse_args(sys.argv[1:])

    if args.mode == 'cov':
        mode = 0
    elif args.mode == 'reads':
        mode = 1
    else:
        print('Please choose either cov or reads for mode')
        exit()

    if not args.chroms:
        args.chroms = 'chroms.' + args.chrom

    if args.profile:
        go_profile(vars(args), mode)
    else:
        go(vars(args), mode)
