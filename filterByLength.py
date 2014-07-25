#! /usr/bin/env python

import sys
import re

def filterByLength(inputFile, outputFile, filterLength):
    with open(inputFile, 'r') as fin:
        with open(outputFile, 'w') as fout:
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
                        fout.write(row[2] + '\t0\t' + '\t'.join(row[2:9]) + '\t*\t*')
                        #fout.write('\t'.join(row[:11]))

                        if str(row[7]) == 0 and str(row[8]) == 0:
                            fout.write('\tXS:A:+')
                        else:
                            for x in row[11:]:
    	                        if x[:3] == 'XS:':
        	                        fout.write('\tXS:A:+')
                                    #fout.write('\t' + x)
                        fout.write('\n')


''' Usage: ./filterByLength input.sam output.sam length
'''
filterByLength(sys.argv[1], sys.argv[2], int(sys.argv[3]))