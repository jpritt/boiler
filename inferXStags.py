#! /usr/bin/env python3

import sys
import re
import random

def printByRange(line):
    chrom = '2L'
    start = 21578320
    end = 21578489

    row = line.rstrip().split('\t')
    if len(row) < 6:
        print('\t'.join(row))
    else:
        row[0] = chrom
        row[9] = '*'
        row[10] = '*'
        id1 = int(row[3])
        chrom1 = row[2]
        if chrom1 == chrom and id1 >= start and id1 <= end:
            print('\t'.join(row))
        else:
            id2 = int(row[7])
            if chrom1 == chrom and row[6] == '=' and id2 >= start and id2 <= end:
                print('\t'.join(row))
            elif row[6] == chrom and id2 >= start and id2 <= end:
                print('\t'.join(row))

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
                return True

    countA = len(exonsA)
    for i in range(countA):
        e = exonsA[countA-i-1]
        if e[1] < exonsB[0][0]:
            break

        for i in range(len(exonsB)-1):
            if e[1] <= exonsB[i][1]:
                break
            elif e[1] > exonsB[i][1] and e[0] < exonsB[i+1][0]:
                return True

    return False

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
                # B overlaps an intron of A
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
                # A overlaps an intron of B
                return 2

    return 0


def filterFlags(flags, filters):
    for f in filters:
        if flags & f:
            flags -= f
    return flags

with open(sys.argv[1], 'r') as f:

    unmatched = dict()
    for line in f:
        row = line.rstrip().split('\t')
        if len(row) < 6:
            print('\t'.join(row))
        else:
            XS = None
            for r in row[11:]:
                if r[:2] == 'XS':
                    XS = r[5]
            
            if not row[6] == '*':
                name = row[0]
                exonsB = parseCigar(row[5], int(row[3]))
                if name in unmatched:
                    foundMatch = False
                    for j in range(len(unmatched[name])):
                        match = unmatched[name][j][0]

                        exonsA = parseCigar(match[5], int(match[3]))
                        if ((row[6] == '=' and match[6] == '=' and row[2] == match[2]) or (row[6] == match[2] and row[2] == match[6])) and row[7] == match[3] and match[7] == row[3]:
                            foundMatch = True

                            #if XS and unmatched[name][j][1] and not (XS == unmatched[name][j][1]):
                            if (not XS and unmatched[name][j][1]):
                                row.append('XS:A:' + unmatched[name][j][1])
                            elif (XS and not unmatched[name][j][1]):
                                match.append('XS:A:' + XS)
                                # unpair
                                #row[6] = '*'
                                #row[7] = '0'
                                #row[8] = '0'
                                #match[6] = '*'
                                #match[7] = '0'
                                #match[8] = '0'

                                # remove XS values
                                #for n in range(len(row)):
                                #    if row[n][:2] == 'XS':
                                #        del row[n]
                                #        break
                                #for n in range(len(match)):
                                #    if match[n][:2] == 'XS':
                                #        del match[n]
                                #        break

                            print('\t'.join(row))
                            print('\t'.join(match))

                            del unmatched[name][j]
                            break
                    if not foundMatch:
                        unmatched[name].append((row, XS))
                else:
                    unmatched[name] = [(row, XS)]
            else:
                print('\t'.join(row))
    
    count = 0
    for k,v in unmatched.items():
        count += len(v)
        for r in v:
            print('\t'.join(r[0]))

