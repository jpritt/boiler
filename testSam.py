#! /usr/bin/env python

import sys

def parseFlags(flagNum):
    flags = [0] * 12
    for i in xrange(1,13):
        num = 2 ** (12-i)
        if flagNum >= num:
            flags[i-1] = 1
            flagNum -= num
    return flags

with open(sys.argv[1], 'r') as f:
    #flagSums = [0]*12
    flagsCountNone = dict()
    flagsCountPlus = dict()
    flagsCountMinus = dict()
    for line in f:
        row = line.rstrip().split('\t')
        if len(row) > 6:
            #if parseFlags(int(row[1]))[9] == 1:
            #    print line.rstrip()

            '''
            if row[6] == '*':
                pass
                #flag = int(row[1])
                #if flag in flagsCount:
                #    flagsCount[flag] += 1
                #else:
                #    flagsCount[flag] = 1
                
                #flags = parseFlags(int(row[1]))
                #for i in xrange(12):
                #    flagSums[i] += flags[i]
                #print row[1] + ':\t' + ','.join([str(c) for c in flags])
            else:
                flag = int(row[1])
                if flag in flagsCount:
                    flagsCount[flag] += 1
                else:
                    flagsCount[flag] = 1

                flagBits = parseFlags(flag)
                if flagBits[5] == 1 and flagBits[8] == 1 and int(row[3]) > int(row[7]):
                    print '\t'.join(row)
                elif flagBits[6] == 1 and flagBits[7] == 1 and int(row[3]) < int(row[7]):
                    print '\t'.join(row)
            '''
            flag = int(row[1])
            for r in row[10:]:
                if r[:3] == 'XS:':
                    if r[-1] == '+':
                        if flag in flagsCountPlus:
                            flagsCountPlus[flag] += 1
                        else:
                            flagsCountPlus[flag] = 1
                    else:
                        if flag in flagsCountMinus:
                            flagsCountMinus[flag] += 1
                        else:
                            flagsCountMinus[flag] = 1
                    break
            if flag in flagsCountNone:
                flagsCountNone[flag] += 1
            else:
                flagsCountNone[flag] = 1


    print 'XS:A:+'
    print sorted(flagsCountPlus.keys())
    #for k,v in flagsCountPlus.items():
    #    print '%d (%d):\t' % (k,v) + ','.join([str(c) for c in parseFlags(k)])

    print 'XS:A:-'
    print sorted(flagsCountMinus.keys())
    #for k,v in flagsCountMinus.items():
    #    print '%d (%d):\t' % (k,v) + ','.join([str(c) for c in parseFlags(k)])

    print 'No XS'
    print sorted(flagsCountNone.keys())
    #for k,v in flagsCountNone.items():
    #    print '%d (%d):\t' % (k,v) + ','.join([str(c) for c in parseFlags(k)])