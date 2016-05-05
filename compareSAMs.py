#! /usr/bin/env python
import sys
import re
import argparse

def conflicts(exonsA, exonsB):
    '''
        Returns true if any of the exons from A or B overlaps one of the introns from the other set of exons
    '''
    
    for e in exonsB:
        if e[0] > exonsA[-1][0]:
            break

        for i in range(len(exonsA)-1):
            if e[0] >= exonsA[-i-1][0]:
                break
            elif e[1] > exonsA[-i-2][1]:
                # Exon in B overlaps an intron in A
                return 1

    countA = len(exonsA)
    for i in range(countA):
        e = exonsA[countA-i-1]
        if e[1] < exonsB[0][1]:
            break

        for i in range(len(exonsB)-1):
            if e[1] <= exonsB[i][1]:
                break
            elif e[1] > exonsB[i][1] and e[0] < exonsB[i+1][0]:
                # Exon in A overlaps an intron in B
                return 2

    return 0

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

        # Skip soft clipping
        if not cigar[index] == 'S':
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

def getChromReads(f, chr, fragment_lengths, countUnpaired, countPaired, countDiscordant, bin=1000):
    reads = []
    unmatched = dict()
    for line in f:
        row = line.rstrip().split('\t')
        if len(row) < 6:
            continue

        flags = int(row[1])
        if (flags & 4):
            continue

        if row[2] == chr:
            exons = parseCigar(row[5],int(row[3]))
            if row[6] == '*' or (flags & 8):
                countUnpaired += 1
            else:
                countPaired += 1

                if row[6] == '=':
                    name = row[0]

                    if name in unmatched:
                        foundMatch = False
                        for i in range(len(unmatched[name])):
                            match = unmatched[name][i]
                            if int(row[3]) == match[2] and int(row[7]) == match[0] and not conflicts(exons, match[1]):
                                #l = abs(int(row[8]))
                                #print(chr)
                                #print(name)
                                #print(exons)
                                #print(match[1])
                                l = max(exons[-1][1], match[1][-1][1]) - min(exons[0][0], match[1][0][0])
                                #print('%d, %d, %d' % (int(row[8]), match[3], l))
                                #exit()

                                if l > 0:
                                    lbin = int(l / bin)
                                    if lbin >= len(fragment_lengths):
                                        fragment_lengths += [0] * (lbin + 1 - len(fragment_lengths))
                                    fragment_lengths[lbin] += 1

                                del unmatched[name][i]
                                foundMatch = True
                                break

                        if not foundMatch:
                            unmatched[name].append((int(row[3]), exons, int(row[7])))
                    else:
                        unmatched[name] = [(int(row[3]), exons, int(row[7]))]

            reads.append((row[2], int(row[3]), genCigar(exons), row[6], int(row[7])))

    for k,v in unmatched.items():
        countDiscordant += len(v)
    return reads, countUnpaired, countPaired, countDiscordant

def compareSAMs(file1, file2):
    fragment_lengths1 = [0] * 100000
    fragment_lengths2 = [0] * 100000
    countUnpaired1 = 0
    countPaired1 = 0
    countUnpaired2 = 0
    countPaired2 = 0
    countDiscordant1 = 0
    countDiscordant2 = 0

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
            reads1, countUnpaired1, countPaired1, countDiscordant1 = getChromReads(f, c, fragment_lengths1, countUnpaired1, countPaired1, countDiscordant1)
        with open(file2, 'r') as f:
            reads2, countUnpaired2, countPaired2, countDiscordant2 = getChromReads(f, c, fragment_lengths2, countUnpaired2, countPaired2, countDiscordant2)

        reads1.sort()
        reads2.sort()

        #print(reads1[:10])
        #print(reads2[:10])
        #exit()

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

        # Set index to 1 when we've matched to a read
        counted = [0] * len(reads2)

        for r in reads1:
            while id < len(reads2) and reads2[id] < r:
                id += 1

            if id < len(reads2) and reads2[id] == r:
                    tp += 1
                    id += 1
            #else:
            #    print(r)

        TP_pairs += tp
        FN_pairs += len(reads1) - tp
        FP_pairs += len(reads2) - tp

    #print(fragment_lengths1)
    #print(fragment_lengths2)

    print('No pairing:')
    #print('TP: %d' % TP_no_pairs)
    #print('FN: %d' % FN_no_pairs)
    #print('FP: %d' % FP_no_pairs)
    print('Recall    = TP / TP+FN = %f' % (float(TP_no_pairs) / (TP_no_pairs+FN_no_pairs)))
    print('Precision = TP / TP+FP = %f' % (float(TP_no_pairs) / (TP_no_pairs+FP_no_pairs)))
    print('')

    print('Pairing:')
    #print('TP: %d' % TP_pairs)
    #print('FN: %d' % FN_pairs)
    #print('FP: %d' % FP_pairs)
    print('Recall    = TP / TP+FN = %f' % (float(TP_pairs) / (TP_pairs+FN_pairs)))
    print('Precision = TP / TP+FP = %f' % (float(TP_pairs) / (TP_pairs+FP_pairs)))
    print('')


    print('File 1: %d unpaired, %d paired, %d discordant' % (countUnpaired1, countPaired1, countDiscordant1))
    print('File 2: %d unpaired, %d paired, %d discordant' % (countUnpaired2, countPaired2, countDiscordant2))
    print('')


    #fragment_lengths1 = bin(fragment_lengths1, 1000)
    #fragment_lengths2 = bin(fragment_lengths2, 1000)
    #print(fragment_lengths1)
    #print(fragment_lengths2)

    return fragment_lengths1, fragment_lengths2

def plot_frag_lens(fragment_lengths1, fragment_lengths2):
    from matplotlib import pyplot as plt

    font = {'size': 18}
    plt.rc('font', **font)

    x_start = 0
    x_end = 100
    f, axs = plt.subplots(1,2,figsize=(20,8))
    a, = axs[1].plot(range(x_start,x_end), fragment_lengths1[x_start:x_end])
    b, = axs[1].plot(range(x_start,x_end), fragment_lengths2[x_start:x_end])
    axs[1].set_xlabel('Genomic Outer Distance (kb)')
    axs[1].set_ylabel('Frequency')
    axs[1].set_yscale('log')
    axs[1].legend([a,b], ['Original', 'Compressed'], prop={'size':14})
    #plt.savefig('frag_len_dist_geuv.pdf', format='pdf')
    #plt.clf()

    x_range = 100
    ratio = [0] * x_range
    for i in range(x_range):
        if fragment_lengths1[i] > 0:
            ratio[i] = float(fragment_lengths2[i]) / float(fragment_lengths1[i])
    axs[0].plot(range(x_range), ratio)
    axs[0].set_xlabel('Genomic Outer Distance (kb)')
    axs[0].set_ylabel('Compressed / Original Frequency')
    plt.savefig('frag_len.pdf', format='pdf', bbox_inches='tight')
    plt.clf()

def write_frag_lens(filename, fragment_lengths1, fragment_lengths2):
    print(len(fragment_lengths1))
    print(len(fragment_lengths2))
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
        frags1, frags2 = compareSAMs(args.sam1, args.sam2)

        if args.out_frags:
            write_frag_lens(args.out_frags, frags1, frags2)
    elif args.out_frags:
        frags1, frags2 = read_frag_lens(args.out_frags)
    else:
        print('Include either 2 SAM files to compare, or the path to a fragment length distribution')
        exit()

    if args.plot:
        plot_frag_lens(frags1, frags2)
