#! /usr/bin/env python
import sys
import re
import argparse

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

def genCigar(exons):
    cigar = [str(exons[0][1] - exons[0][0]) + 'M']
    for i in range(1, len(exons)):
        if exons[i][0] - exons[i-1][1] == 0:
            prevLen = int(cigar[-1][:-1])
            cigar[-1] = str(prevLen + exons[i][1] - exons[i][0]) + 'M'
        else:
            cigar += [str(exons[i][0] - exons[i-1][1]) + 'N']
            cigar += [str(exons[i][1] - exons[i][0]) + 'M']
    return ''.join(cigar)

def bin(data, width):
    new_data = [0] * int(len(data) / width)
    for i in range(len(new_data)):
        new_data[i] = sum(data[i*width:(i+1)*width])
    return new_data

def getChromReads(f, chr, fragment_lengths, countUnpaired, countPaired):
    reads = []
    for line in f:
        row = line.rstrip().split('\t')
        if len(row) < 6:
            continue

        if row[2] == chr:
            if row[6] == '*':
                countUnpaired += 1
            else:
                countPaired += 1

                if row[6] == '=':
                    l = int(row[8])

                    if l > 0:
                        if l >= len(fragment_lengths):
                            fragment_lengths += [0] * (l + 1 - len(fragment_lengths))
                        fragment_lengths[l] += 1


            reads.append((row[2], int(row[3]), genCigar(parseCigar(row[5],0)), row[6], int(row[7])))

    return reads, countUnpaired, countPaired

def compareSAMs(file1, file2, plot):
    fragment_lengths1 = [0] * 100000
    fragment_lengths2 = [0] * 100000
    countUnpaired1 = 0
    countPaired1 = 0
    countUnpaired2 = 0
    countPaired2 = 0

    chroms = []

    with open(file1, 'r') as f:
        for line in f:
            row = line.rstrip().split('\t')
            if line[0] == '@':
                if row[0] == '@SQ':
                    chroms.append(row[1].split(':')[1])
            else:
                break

    TP_no_pairs = 0
    FN_no_pairs = 0
    FP_no_pairs = 0
    TP_pairs = 0
    FN_pairs = 0
    FP_pairs = 0
    for c in chroms:
        with open(file1, 'r') as f:
            reads1, countUnpaired1, countPaired1 = getChromReads(f, c, fragment_lengths1, countUnpaired1, countPaired1)
        with open(file2, 'r') as f:
            reads2, countUnpaired2, countPaired2 = getChromReads(f, c, fragment_lengths2, countUnpaired2, countPaired2)

        reads1.sort()
        reads2.sort()

        # No pairing
        tp = 0
        id = 0
        for read in reads1:
            r = read[:3]
            while id < len(reads2) and reads2[id][:3] < r:
                id += 1

            if id < len(reads2) and reads2[id][:3] == r:
                tp += 1
                id += 1
            #elif id > 0 and id < len(reads2):
            #    print(r)
            #    print(reads2[id-1])
            #    print(reads2[id])

        TP_no_pairs += tp
        FN_no_pairs += len(reads1) - tp
        FP_no_pairs += len(reads2) - tp

        # Pairing
        tp = 0
        id = 0
        for r in reads1:
            while id < len(reads2) and reads2[id] < r:
                id += 1

            if id < len(reads2) and reads2[id] == r:
                    tp += 1
                    id += 1
        TP_pairs += tp
        FN_pairs += len(reads1) - tp
        FP_pairs += len(reads2) - tp

    print('No pairing:')
    print('TP: %d' % TP_no_pairs)
    print('FN: %d' % FN_no_pairs)
    print('FP: %d' % FP_no_pairs)
    print('Recall    = TP / TP+FN = %f' % (float(TP_no_pairs) / (TP_no_pairs+FN_no_pairs)))
    print('Precision = TP / TP+FP = %f' % (float(TP_no_pairs) / (TP_no_pairs+FP_no_pairs)))
    print('')

    print('Pairing:')
    print('TP: %d' % TP_pairs)
    print('FN: %d' % FN_pairs)
    print('FP: %d' % FP_pairs)
    print('Recall    = TP / TP+FN = %f' % (float(TP_pairs) / (TP_pairs+FN_pairs)))
    print('Precision = TP / TP+FP = %f' % (float(TP_pairs) / (TP_pairs+FP_pairs)))
    print('')


    print('File 1: %d unpaired, %d paired' % (countUnpaired1, countPaired1))
    print('File 2: %d unpaired, %d paired' % (countUnpaired2, countPaired2))
    print('')

    fragment_lengths1 = bin(fragment_lengths1, 1000)
    fragment_lengths2 = bin(fragment_lengths2, 1000)
    return fragment_lengths1, fragment_lengths2

def plot_frag_lens(fragment_lengths1, fragment_lengths2):
    from matplotlib import pyplot as plt

    x_range = 100
    a, = plt.plot(range(x_range), fragment_lengths1[:x_range])
    b, = plt.plot(range(x_range), fragment_lengths2[:x_range])
    plt.xlabel('Fragment Length (kb)')
    plt.ylabel('Frequency')
    plt.yscale('log')
    plt.legend([a,b], ['Original', 'Compressed'])
    plt.savefig('frag_len_dist.png')
    plt.clf()

    diffs = [0] * x_range
    for i in range(x_range):
        diffs[i] = fragment_lengths1[i] - fragment_lengths2[i]

    plt.plot(range(x_range), diffs)
    plt.xlabel('Fragment Length (kb)')
    plt.ylabel('Original - Compressed Frequency')
    plt.savefig('frag_len_diff.png')
    plt.clf()

def write_frag_lens(filename, fragment_lengths1, fragment_lengths2):
    with open(filename, 'w') as f:
        f.write(','.join([str(l) for l in fragment_lengths1]) + '\n')
        f.write(','.join([str(l) for l in fragment_lengths2]) + '\n')

def read_frag_lens(filename):
    with open(filename, 'r') as f:
        fragment_lengths1 = [int(l) for l in f.readline().rstrip().split(',')]
        fragment_lengths2 = [int(l) for l in f.readline().rstrip().split(',')]
    return fragment_lengths1, fragment_lengths2

if __name__ == '__main__':

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sam1', type=str, help='Full path of first SAM file to compare')
    parser.add_argument('--sam2', type=str, help='Full path of first SAM file to compare')
    parser.add_argument('--out-frags', type=str, help="File to write fragment length distributions to")
    parser.add_argument('--plot', help="Plot fragment length distribution graphs", action="store_true")

    args = parser.parse_args(sys.argv[1:])

    if args.sam1 and args.sam2:
        frags1, frags2 = compareSAMs(args.sam1, args.sam2, args.plot)

        if args.out_frags:
            write_frag_lens(args.out_frags, frags1, frags2)
    elif args.out_frags:
        frags1, frags2 = read_frag_lens(args.out_frags)
    else:
        print('Include either 2 SAM files to compare, or the path to a fragment length distribution')
        exit()

    if args.plot:
        plot_frag_lens(frags1, frags2)
