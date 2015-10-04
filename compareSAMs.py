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

def compareSAMs(file1, file2):
    fragment_lengths1 = [0] * 100000
    reads1 = []
    countUnpaired = 0
    countPaired = 0

    header = True
    chromOffsets = dict()
    covLen = 0

    with open(file1, 'r') as f:
        for line in f:
            row = line.rstrip().split('\t')
            if header and line[0] == '@':
                if row[0] == '@SQ':
                    chrom = row[1].split(':')[1]
                    length = int(row[2].split(':')[1])

                    chromOffsets[chrom] = length
                    covLen += length

            else:
                if header:
                    # End of header
                    header = False
                    cov1 = [0] * covLen
                    cov2 = [0] * covLen

                if len(row) < 6:
                    continue

                if row[6] == '*':
                    countUnpaired += 1
                else:
                    countPaired += 1

                    if row[6] == '=':
                        l = int(row[8])

                        if l > 0:
                            if l >= len(fragment_lengths1):
                                fragment_lengths1 += [0] * (l + 1 - len(fragment_lengths1))
                            fragment_lengths1[l] += 1

                exons = parseCigar(row[5], chromOffsets[row[2]]+int(row[3]))
                for e in exons:
                    for c in range(e[0], e[1]):
                        cov1[c] += 1

                read = (row[2], int(row[3]), row[5], row[6], int(row[7]))
                reads1.append(read)
    print('%d unpaired, %d paired' % (countUnpaired, countPaired))

    fragment_lengths2 = [0] * len(fragment_lengths1)
    reads2 = []
    countUnpaired = 0
    countPaired = 0
    with open(file2, 'r') as f:
        for line in f:
            row = line.rstrip().split('\t')

            if len(row) < 6:
                continue

            if row[6] == '*':
                countUnpaired += 1
            else:
                countPaired += 1

                if row[6] == '=':
                    l = int(row[8])

                    if l > 0:
                        if l >= len(fragment_lengths2):
                            fragment_lengths2 += [0] * (l + 1 - len(fragment_lengths2))
                            fragment_lengths1 += [0] * (l + 1 - len(fragment_lengths1))
                        fragment_lengths2[l] += 1

            exons = parseCigar(row[5], chromOffsets[row[2]]+int(row[3]))
            for e in exons:
                for c in range(e[0], e[1]):
                    cov2[c] += 1

            read = (row[2], int(row[3]), row[5], row[6], int(row[7]))
            reads2.append(read)
    print('%d unpaired, %d paired' % (countUnpaired, countPaired))
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

    reads1.sort()
    reads2.sort()

    print('No pairing:')
    tp = 0
    id = 0
    for r in reads1:
        while id < len(reads2) and reads2[id][0] == r[0] and reads2[id][1] < r[1]:
            id += 1

        if reads2[id][0] == r[0] and reads2[id][1] == r[1] and reads2[id][2] == r[2]:
            tp += 1
            id += 1

    print('TP: %d' % tp)
    print('FN: %d' % (len(reads1)-tp))
    print('FP: %d' % (len(reads2)-tp))
    print('Recall    = TP / TP+FN = %f' % (float(tp) / len(reads1)))
    print('Precision = TP / TP+FP = %f' % (float(tp) / len(reads2)))
    print('')

    print('Pairing:')
    tp = 0
    id = 0
    for r in reads1:
        while id < len(reads2) and reads2[id] < r:
            id += 1

        if reads2[id] == r:
                tp += 1
                id += 1

    print('TP: %d' % tp)
    print('FN: %d' % (len(reads1)-tp))
    print('FP: %d' % (len(reads2)-tp))
    print('Recall    = TP / TP+FN = %f' % (float(tp) / len(reads1)))
    print('Precision = TP / TP+FP = %f' % (float(tp) / len(reads2)))
    print('')

compareSAMs(sys.argv[1], sys.argv[2])
