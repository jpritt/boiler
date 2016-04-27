#! /usr/bin/env python

# Convert a SAM file output by HISAT to one that Boiler can use
# Input file should be sorted by read name

import argparse
import sys

def writeReads(f_out, readsA, readsB, NH1, NH2):
    if NH1 == 0 or NH2 == 0:
        print('NH Error!')

    if not (NH1 == len(readsA)) and (NH2 == len(readsB)):
        print("NH values don't match up with number of reads!")

    chromsA = [r[2] for r in readsA]
    startsA = [int(r[3]) for r in readsA]
    chromsB = [r[2] for r in readsB]
    startsB = [int(r[3]) for r in readsB]
    NH = str(NH1 * NH2)

    for r in readsA:
        for i in range(11, len(r)):
            if r[i][:3] == 'NH:':
                r[i] = 'NH:i:' + NH
                break

        for i in range(NH2):
            if chromsB[i] == r[2]:
                r[6] = '*'
            else:
                r[6] = chromsB[i]
            r[7] = startsB[i]
            f_out.write('\t'.join(r) + '\n')

    for r in readsB:
        for i in range(11, len(r)):
            if r[i][:3] == 'NH:':
                r[i] = 'NH:i:' + NH
                break

        for i in range(NH1):
            if chromsA[i] == r[2]:
                r[6] = '*'
            else:
                r[6] = chromsA[i]
            r[7] = startsA[i]
            f_out.write('\t'.join(r) + '\n')

def processHISAT(f_in, f_out):
    currName = None
    readsA = []
    readsB = []
    NH1 = 0
    NH2 = 0
    for line in f_in:
        if line[0] == '@':
            f_out.write(line)
            continue

        row = line.rstrip().split('\t')
        if not currName == row[0]:
            if readsA:
                writeReads(f_out, readsA, readsB, NH1, NH2)
            readsA = []
            readsB = []
            NH1 = 0
            NH2 = 0

        flags = int(row[1])
        NH = 0
        for r in row[11:]:
            if r[:3] == 'NH:':
                NH = int(r[5:])
                break

        if NH == 0:
            print('NH not found error!')
        if not NH1 and (flags & 64):
            NH1 = NH
            readsA.append(row)
        elif not NH2 and (flags & 128):
            NH2 = NH
            readsB.append(row)
        else:
            print('Read does not have a template flag set!')

    if readsA:
        writeReads(f_out, readsA, readsB, NH1, NH2)


if __name__ == '__main__':

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', type=str, required=True, help='Full path of SAM file to process')
    parser.add_argument('--output', type=str, required=True, help='Path and filename of SAM file to write')
    args = parser.parse_args(sys.argv[1:])

    with open(args.input, 'r') as f_in:
        with open(args.output, 'w') as f_out:
            processHISAT(f_in, f_out)