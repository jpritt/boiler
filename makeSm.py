#! /usr/bin/env python

import sys
import re

def makeSm(inputFile, outputFile, lines):
    ''' Truncate the file to the given number of lines
    '''

    with open(inputFile, 'r') as fin:
        with open(outputFile, 'w') as fout:
            lineNum = 0
            for line in fin:
                fout.write(line)

                lineNum += 1
                if lineNum >= lines:
                    return


''' Usage: ./filterByLength input.sam output.sam length
'''
makeSm(sys.argv[1], sys.argv[2], int(sys.argv[3]))