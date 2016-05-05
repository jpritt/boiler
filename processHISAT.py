#! /usr/bin/env python

# Convert a SAM file output by HISAT to one that Boiler can use
# Input file should be sorted by read name

import argparse
import sys

def oneToOne(readsA, readsB):
    if not len(readsA) == len(readsB):
        return False

    rA = [(int(r[3]), int(r[7])) for r in readsA]
    rB = [(int(r[7]), int(r[3])) for r in readsB]

    for r in rA:
        if r in rB:
            i = rB.index(r)
            del rB[i]
        else:
            return False
    return True

def writeReads(f_out, readsA, readsB):
    numA = len(readsA)
    numB = len(readsB)
    if numA == 0:
        print('Error! No left mates found for %s' % readsB[0][0])
        exit()
    if numB == 0:
        print('Error! No right mates found for %s' % readsA[0][0])
        exit()

    #if not (numA == numB) and numA > 1 and numB > 1:
    '''
    if numA == 2 and (numA == numB):
        for r in readsA:
            print(r[0] + '\t' + r[2] + '\t' + r[3] + '\t\t' + r[6] + '\t' + r[7])
        print('')
        for r in readsB:
            print(r[0] + '\t' + r[2] + '\t' + r[3] + '\t\t' + r[6] + '\t' + r[7])
        print('')
        print('')        

        exit()
    '''

    chromsA = [r[2] for r in readsA]
    startsA = [r[3] for r in readsA]
    chromsB = [r[2] for r in readsB]
    startsB = [r[3] for r in readsB]
    NH = str(numA * numB)

    for r in readsA:
        # Set NH value
        NH_set = False
        for i in range(11, len(r)):
            if r[i][:3] == 'NH:':
                r[i] = 'NH:i:' + NH
                NH_set = True
                break
        if not NH_set:
            r.append('NH:i:' + NH)

        # Set cigar string if necessary
        if r[5] == '*':
            for i in range(numB):
                if readsB[i][2] == r[2] and readsB[i][3] == r[3] and readsB[i][7] == r[7] and readsB[i][3] == r[7]:
                    r[5] = readsB[i][5]
                    break
            if r[5] == '*':
                print('Unable to set cigar string for %s' % r[0])
                exit()

        # Pair with all mates
        for i in range(numB):
            if chromsB[i] == r[2]:
                r[6] = '='
            else:
                r[6] = chromsB[i]
            r[7] = startsB[i]
            f_out.write('\t'.join(r) + '\n')

    for r in readsB:
        # Set NH value
        NH_set = False
        for i in range(11, len(r)):
            if r[i][:3] == 'NH:':
                r[i] = 'NH:i:' + NH
                NH_set = True
                break
        if not NH_set:
            r.append('NH:i:' + NH)

        # Set cigar string if necessary
        if r[5] == '*':
            for i in range(numA):
                if readsA[i][2] == r[2] and readsA[i][3] == r[3] and readsA[i][7] == r[7] and readsA[i][3] == r[7]:
                    r[5] = readsA[i][5]
                    break
            if r[5] == '*':
                print('Unable to set cigar string for %s' % r[0])
                exit()

        # Pair with all mates
        for i in range(numA):
            if chromsA[i] == r[2]:
                r[6] = '='
            else:
                r[6] = chromsA[i]
            r[7] = startsA[i]
            f_out.write('\t'.join(r) + '\n')

def processHISAT(f_in, f_out):
    currName = None
    readsA = []
    readsB = []
    last = None
    for line in f_in:
        if line[0] == '@':
            f_out.write(line)
            continue

        row = line.rstrip().split('\t')
        if not currName == row[0]:
            if readsA:
                if oneToOne(readsA, readsB):
                    NH = str(len(readsA))
                    for r in readsA:
                        for i in range(11, len(r)):
                            if r[i][:2] == 'NH':
                                r[i] = 'NH:i:' + NH
                                break
                        f_out.write('\t'.join(r) + '\n')
                    NH = str(len(readsB))
                    for r in readsB:
                        for i in range(11, len(r)):
                            if r[i][:2] == 'NH':
                                r[i] = 'NH:i:' + NH
                                break
                        f_out.write('\t'.join(r) + '\n')
                else:
                    writeReads(f_out, readsA, readsB)
            currName = row[0]
            readsA = []
            readsB = []

        if row[2] == '*' and row[3] == '0' and row[5] == '*':
            continue

        if row[6] == '*':
            #print(row[0])
            #exit()
            f_out.write(line)

        '''
        if row[5] == '*':
            if not row[2] == last[2] and row[6] == last[6]:
                print(row[0])
                exit()
            row[5] = last[5]
            for r in row[11:]:
                if r[:2] == 'XS':
                    print('Row has XS')
                    print(line)
                    exit()
            for r in last[11:]:
                if r[:2] == 'XS':
                    row.append(r)
        '''

        flags = int(row[1])

        if (flags & 64):
            readsA.append(row)
        elif (flags & 128):
            readsB.append(row)
        else:
            print('Read does not have a template flag set!')
            print(line)
            exit()

        last = row[:]

    if readsA:
        if oneToOne(readsA, readsB):
            NH = str(len(readsA))
            for r in readsA:
                for i in range(11, len(r)):
                    if r[i][:2] == 'NH':
                        r[i] = 'NH:i:' + NH
                        break
                f_out.write('\t'.join(r) + '\n')
            NH = str(len(readsB))
            for r in readsB:
                for i in range(11, len(r)):
                    if r[i][:2] == 'NH':
                        r[i] = 'NH:i:' + NH
                        break
                f_out.write('\t'.join(r) + '\n')
        else:
            writeReads(f_out, readsA, readsB)


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
