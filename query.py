#! /usr/bin/env python3
import sys
import argparse

import readSAM
import expand
import readPRO
import time
import re
import math
import random
#import matplotlib.pyplot as plt

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

    expander = expand.Expander()

    sam = readSAM.ReadSAM(alignmentsFile, chromosomes)
    expander = expand.Expander()
    pro = readPRO.ReadPRO(args['pro'])

    '''
    trueCov = sam.getCoverage('3R', 24150663, 24151663)
    predCov = expander.getCoverage(args['compressed'], '3R', 24150663, 24151663)
    print(trueCov[548])
    print(predCov[548])

    for x in range(len(trueCov)):
        if abs(trueCov[x] - predCov[x]) > 0.0001:
            print('Error!')
            print('%s (%d, %d)' % (chrom, start, start+l))
            print(x)
            for n in range(x-5,x+5):
                print(str(trueCov[n]) + '\t' + str(predCov[n]))
            exit()
    exit()
    '''

    #s = time.time()
    #trueCov = sam.getCoverage('2L', 15609942, 16609942)
    #predCov = expander.getCoverage(args['compressed'], '2L', 15609942, 16609942)

    #for x in range(len(trueCov)):
    #    if abs(trueCov[x] - predCov[x]) > 0.0001:
    #        print('Error!')
    #        print(x)
    #        for n in range(x-5,x+5):
    #            print(str(trueCov[n]) + '\t' + str(predCov[n]))
    #        exit()
    #e = time.time()
    #print('%f s' % (e-s))
    #exit()

    lens = [1000, 10000, 100000, 1000000, 10000000, 20000000]
    #lens = [10000000]
    chroms = ['2R', '2L', '3R', '3L', 'X']
    chromLens = [21146708, 23011544, 27905053, 24543557, 22422827]


    print('Querying coverage')

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


    print('Querying reads')


    # Index i in these lists corresponds to the range of lengths [10^i, 10^(i+1))
    true_times = []
    pred_times = []
    counts = []

    pro = readPRO.ReadPRO(args['pro'])

    num_genes = 0
    for c in range(len(chroms)):
        chrom = chroms[c]

        genes = pro.getGenes(chrom)
        num_genes += len(genes)
        print('%d genes in chromosome %s' % (len(genes), chrom))

        offset = 0
        for k in sorted(chromosomes.keys()):
            if not k == chrom:
                offset += chromosomes[k]
            else:
                break

        for g in genes:
            start = g[0]
            end = g[1]

            bin = int(math.log(end - start,10))
            if bin > len(true_times)-1:
                add = bin - len(true_times) + 1
                true_times += [0] * add
                pred_times += [0] * add
                counts += [0] * add

            startTime = time.time()
            trueUnpaired, truePaired = sam.getReads(chrom, start, end)
            endTime = time.time()

            true_times[bin] += float(endTime - startTime)

            startTime = time.time()
            predUnpaired, predPaired = expander.getReads(args['compressed'], chrom, start, end)
            endTime = time.time()


            pred_times[bin] += float(endTime - startTime)
            counts[bin] += 1

            pred = predUnpaired[:]
            correctUnpaired = 0
            for r in trueUnpaired:
                for i in range(len(pred)):
                    if r == pred[i]:
                        correctUnpaired += 1
                        del pred[i]
                        break

            pred = predPaired[:]
            correctPaired = 0
            for r in truePaired:
                for i in range(len(pred)):
                    if r == pred[i]:
                        correctPaired += 1
                        del pred[i]
                        break

            #print('Unpaired:')
            #print('%d true reads, %d predicted, %d correct' % (len(trueUnpaired), len(predUnpaired), correctUnpaired))
            #print('Paired:')
            #print('%d true reads, %d predicted, %d correct' % (len(truePaired), len(predPaired), correctPaired))

            '''
            if (len(trueUnpaired) + len(truePaired) + len(predUnpaired) + len(predPaired)) > 0:
                cov_true_total = 0
                for r in trueUnpaired:
                    for e in r:
                        cov_true_total += e[1] - e[0]
                for p in truePaired:
                    for r in p:
                        for e in r:
                            cov_true_total += e[1] - e[0]

                cov_pred_total = 0
                for r in predUnpaired:
                    for e in r:
                        cov_pred_total += e[1] - e[0]
                for p in predPaired:
                    for r in p:
                        for e in r:
                            cov_pred_total += e[1] - e[0]

                print('Total Coverage:  %d, %d' % (cov_true_total, cov_pred_total))

                start += offset
                end += offset
                cov_true_region = 0
                for r in trueUnpaired:
                    for e in r:
                        if e[1] > start and e[0] < end:
                            cov_true_region += min(end, e[1]) - max(start, e[0])
                for p in truePaired:
                    for r in p:
                        for e in r:
                            if e[1] > start and e[0] < end:
                                cov_true_region += min(end, e[1]) - max(start, e[0])

                cov_pred_region = 0
                for r in predUnpaired:
                    for e in r:
                        if e[1] > start and e[0] < end:
                            cov_pred_region += min(end, e[1]) - max(start, e[0])
                for p in predPaired:
                    for r in p:
                        for e in r:
                            if e[1] > start and e[0] < end:
                                cov_pred_region += min(end, e[1]) - max(start, e[0])

                print('Region Coverage: %d, %d' % (cov_true_region, cov_pred_region))

                true_avg_len = 0
                for r in trueUnpaired:
                    true_avg_len += r[-1][1] - r[0][0]
                for p in truePaired:
                    true_avg_len += p[1][-1][1] - p[0][0][0]
                true_avg_len /= float(len(trueUnpaired) + len(truePaired))

                pred_avg_len = 0
                for r in predUnpaired:
                    pred_avg_len += r[-1][1] - r[0][0]
                for p in predPaired:
                    pred_avg_len += p[1][-1][1] - p[0][0][0]
                pred_avg_len /= float(len(predUnpaired) + len(predPaired))

                print('Average Length:  %d, %d' % (true_avg_len, pred_avg_len))

            print('')
            '''

    for i in range(len(true_times)):
        if counts[i] > 0:
            true_times[i] /= float(counts[i])
            pred_times[i] /= float(counts[i])

    '''
    w = 0.33
    xs = [0] * len(true_times)
    for i in range(len(true_times)):
        xs[i] = i - w

    # plot times
    #a, = plt.plot(xs, true_times)
    #b, = plt.plot(xs, pred_times)

    a = plt.bar(xs, true_times, width=w, color='b')
    for i in range(len(xs)):
        xs[i] = i
    b = plt.bar(xs, pred_times, width=w, color='r')
    plt.xlabel('Region Length')
    plt.ylabel('Avg Time (s)')
    plt.xlim([-w, xs[-1]+w])
    plt.legend([a,b], ['True', 'Compressed'], loc=2)
    plt.title('Average Read Query Time')
    plt.savefig('read_query_time.png')
    plt.clf()
    '''

    with open('times.txt', 'w') as f:
        f.write('\t'.join([str(t) for t in true_times]))
        f.write('\n')
        f.write('\t'.join([str(t) for t in pred_times]))
        f.write('\n')
        f.write('\t'.join([str(t) for t in counts]))
        f.write('\n')
        
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
