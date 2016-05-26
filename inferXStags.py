#! /usr/bin/env python3

import sys
import re
import random

'''
For pairs with one mate assigned an XS tag and the other with no XS tag, add that tag to the untagged mate as well
'''

def inferTags(input, output):
    with open(input, 'r') as f:
        with open(output, 'w') as f2:
            unmatched = dict()
            for line in f:
                row = line.rstrip().split('\t')
                if line[0] == '@':
                    f2.write(line)
                else:
                    XS = None
                    for r in row[11:]:
                        if r[:2] == 'XS':
                            XS = r[5]
            
                    if not row[6] == '*':
                        name = row[0]
                        if name in unmatched:
                            foundMatch = False
                            for j in range(len(unmatched[name])):
                                match = unmatched[name][j][0]

                                if ((row[6] == '=' and match[6] == '=' and row[2] == match[2]) or (row[6] == match[2] and row[2] == match[6])) and row[7] == match[3] and match[7] == row[3]:
                                    foundMatch = True

                                    #if XS and unmatched[name][j][1] and not (XS == unmatched[name][j][1]):
                                    if (not XS and unmatched[name][j][1]):
                                        row.append('XS:A:' + unmatched[name][j][1])
                                    elif (XS and not unmatched[name][j][1]):
                                        match.append('XS:A:' + XS)

                                    f2.write('\t'.join(row) + '\n')
                                    f2.write('\t'.join(match) + '\n')

                                    del unmatched[name][j]
                                    break
                            if not foundMatch:
                                unmatched[name].append((row, XS))
                        else:
                            unmatched[name] = [(row, XS)]
                    else:
                        f2.write(line)
    
            count = 0
            for k,v in unmatched.items():
                count += len(v)
                for r in v:
                    f2.writet('\t'.join(r[0]) + '\n')

if __name__ == '__main__':
    inferTags(sys.argv[1], sys.argv[2])
