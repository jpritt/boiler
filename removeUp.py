#! /usr/bin/env python

import sys
import re

def removeUnmapped(input, output):
    with open(input, 'r') as f:
        with open(output, 'w') as f2:
            for line in f:
                if line[0] == '@':
                    f2.write(line)
                else:
                    row = line.rstrip().split('\t')
                    if row[2] == '*' and row[3] == '0':
                        continue
                    for i in range(11, len(row)):
                        if row[i][:2] == 'YT':
                            if row[i][5:] == 'UP':
                                del row[i]
                                break
                    f2.write('\t'.join(row) + '\n')

if __name__ == '__main__':
    removeUnmapped(sys.argv[1], sys.argv[2])
