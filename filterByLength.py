#! /usr/bin/env python

import sys
import re

def filterByLength(inputFile, outputFile, filterLength):
    with open(inputFile, 'r') as fin:
        with open(outputFile, 'w') as fout:
            i = 0
            for line in fin:
                row = line.rstrip().split('\t')

                if len(row) < 6:
                    fout.write(line)
                else:
                    # Parse cigar string
                    cigar = row[5]
                    origCigar = cigar
                    totalLength = 0
                    match = re.search("\D", cigar)
                    while match:
                        index = match.start()
                        length = int(''.join(cigar[:index]))

                        if cigar[index] == 'M' or cigar[index] == 'D':
                            totalLength += length

                        cigar = cigar[index+1:]
                        match = re.search("\D", cigar)


                    #if totalLength == filterLength:
                    if True:
                        fout.write(str(i) + '\t' + '0' + '\t' + row[2] + '\t' + row[3] + '\t0\t' + row[5] + '\t' + row[6] + '\t' + row[7] + '\t' + row[8] + '\t*\t*')
                        
                        for x in row[11:]:
    	                    if x[:3] == 'XS:':
        	                    fout.write('\tXS:A:+')
                            elif x[:3] == 'NH:':
                                fout.write('\t' + x)
                        fout.write('\n')
                i += 1


''' Usage: ./filterByLength input.sam output.sam length
'''
filterByLength(sys.argv[1], sys.argv[2], int(sys.argv[3]))