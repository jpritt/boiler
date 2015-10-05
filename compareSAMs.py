#! /usr/bin/env python
import sys
import re
from matplotlib import pyplot as plt

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


            reads.append((row[2], int(row[3]), row[5], row[6], int(row[7])))

    return reads, countUnpaired, countPaired

def compareSAMs(file1, file2):
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
        for r in reads1:
            while id < len(reads2) and reads2[id][0] == r[0] and reads2[id][1] < r[1]:
                id += 1

            if reads2[id][0] == r[0] and reads2[id][1] == r[1] and reads2[id][2] == r[2]:
                tp += 1
                id += 1

        TP_no_pairs += tp
        FN_no_pairs += len(reads1) - tp
        FP_no_pairs += len(reads2) - tp

        tp = 0
        id = 0
        for r in reads1:
            while id < len(reads2) and reads2[id] < r:
                id += 1

            if reads2[id] == r:
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

    '''
    for i in range(min(len(fragment_lengths1), len(fragment_lengths2))):
        print('%d\t%d' % (fragment_lengths1[i], fragment_lengths2[i]))
    if len(fragment_lengths1) > len(fragment_lengths2):
        for i in range(len(fragment_lengths2), len(fragment_lengths1)):
            print('%d\t%d' % (fragment_lengths1[i], 0))
    elif len(fragment_lengths2) > len(fragment_lengths1):
        for i in range(len(fragment_lengths1), len(fragment_lengths2)):
            print('%d\t%d' % (0, fragment_lengths2[i]))
    print('')
    '''

    #for i in range(30):
    #    print('%d\t%d' % (fragment_lengths1[i], fragment_lengths2[i]))

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

compareSAMs(sys.argv[1], sys.argv[2])
