#! /usr/bin/env python3

import sys
import re
import random

'''
For pairs with one mate assigned an XS tag and the other with no XS tag, add that tag to the untagged mate as well
'''

with open(sys.argv[1], 'r') as f:
    unmatched = dict()
    for line in f:
        row = line.rstrip().split('\t')
        if line[0] == '@':
            print('\t'.join(row))
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

