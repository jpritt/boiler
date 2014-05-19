#! /usr/bin/env python
import sys
import csv
import math
import re
import string
import numpy as np
import pickle
import sPickle


def compareGTFs(file1, file2):
    lines1 = []
    with open(file1, 'r') as tsv:
        for line in tsv:
            row = line.strip().split('\t')
            if len(row) >= 5:
                covIndex = row[8].find('cov')
                covStart = row[8].find('"', covIndex) + 1
                covEnd = row[8].find('"', covStart)

                lines1.append(row[:3] + [int(row[3])] + [int(row[4])] + [float(row[8][covStart:covEnd])])

    lines2 = []
    with open(file2, 'r') as tsv:
        for line in tsv:
            row = line.strip().split('\t')
            if len(row) >= 5:
                lines2.append(row[:3] + [int(row[3])] + [int(row[4])])

    compareAll(lines1, lines2)

def compareAll(lines1, lines2):
    ''' Compare 
    '''
    numLines1 = len(lines1)
    numLines2 = len(lines2)
    
    testedLines1 = 0.0
    matching = 0.0

    #numTranscripts = 0

    covThreshold = 0

    print "%d lines in file 1" % numLines1
    print "%d lines in file 2" % numLines2

    index = 0
    for l1 in lines1:
        index += 1
        if index % 1000 == 0:
            print "%d / %d lines done" % (index, numLines1)

        if l1[5] <= covThreshold:
            continue

        testedLines1 += 1
        closestScore = 0.0
        
        for l2 in lines2:
            if l1[0] == l2[0] and l1[2] == l2[2]:
                if l1[2] == 'transcript':
                    threshold = 100
                else:
                    threshold = 10
                score = exonsMatch((l1[3], l1[4]), (l2[3], l2[4]), threshold)
                if score > closestScore:
                    closestScore = score
                    closestMatch = l2
        if closestScore > 0:
            matching += closestScore
            #matching += l1[5]

    #print str(numTranscripts) + ' transcripts'
    print 'Using %d of %d lines' % (testedLines1, numLines1)
    print 'TP = ' + str(matching)
    print 'T  = ' + str(testedLines1)
    #print 'P  = ' + str(numLines2)

    #print 'Precision = TP/P = ' + str(matching / numLines2)
    print 'Recall    = TP/T = ' + str(matching / testedLines1)


def exonsMatch(exon1, exon2, threshold):
    ''' Returns a value between 0 and 1, where 1 indicates that the 2 exons are a perfect match and 0 indicates no match
    '''

    dx = min(abs(exon1[0] - exon2[0]), threshold)
    dy = min(abs(exon1[1] - exon2[1]), threshold)

    return 1 - dx / (2.0*threshold) - dy / (2.0*threshold)


compareGTFs(sys.argv[1], sys.argv[2])



''' Unfinished/Old functions '''

def compareGTFsExact(file1, file2):
    ''' Find the number of lines with the exact same chromosome, start, and end indexes between the 2 files
    '''

    reads1 = []
    with open(file1, 'r') as tsv:
        for line in tsv:
            row = line.strip().split('\t')
            if len(row) >= 6:
                reads1.append( (row[0], row[2], row[3], row[4]) )
                #reads1.append( (row[0], row[2], row[3], row[4], row[5], row[6], row[7], row[8]) )
    
    reads2 = []
    with open(file2, 'r') as tsv:
        for line in tsv:
            row = line.strip().split('\t')
            if len(row) >= 6:
                reads2.append( (row[0], row[2], row[3], row[4]) )
                #reads2.append( (row[0], row[2], row[3], row[4], row[5], row[6], row[7], row[8]) )

    print str(len(reads1)) + ' lines in file 1'
    print str(len(reads2)) + ' lines in file 2'

    reads1 = np.array(reads1)
    order = reads1[:,0].argsort()
    reads1 = np.take(reads1, order, 0)
    reads1 = reads1.tolist()

    reads2 = np.array(reads2)
    order = reads2[:,0].argsort()
    reads2 = np.take(reads2, order, 0)
    reads2 = reads2.tolist()

    matches = 0

    f = open('unmatched1.txt', 'w')
    for read in reads1:
        if read in reads2:
            matches += 1
            reads2.remove(read)
        else:
            f.write(str(read[0]) + '\t' + str(read[1]) + '\t' + str(read[2]) + '\t' + str(read[3]) + '\n')
    f.close()

    f = open('unmatched2.txt', 'w')
    for read in reads2:
        f.write(str(read[0]) + '\t' + str(read[1]) + '\t' + str(read[2]) + '\t' + str(read[3]) + '\n')
    f.close()

    print str(matches) + ' matches'


def parseGTF(file):
    ''' Get list of transcripts and exons from a gtf file
    '''
    transcripts = []
    exons = []

    with open(file1, 'r') as tsv:
        for line in tsv:
            row = line.strip().split('\t')
            if len(row) >= 6:
                if row[3] == 'transcript':
                    transcripts.append([])
                else:
                    exons.append( (row[0], row[3], row[4]) )
                    transcripts[-1].append(len(exons)-1)
    return transcripts, exons


def compareExons(exons1, exons2):
    ''' Matches the exons in the 2 lists, returns the distance between all the exons and the corresponding indices of matching exons
    '''

    # Indicates the index in exons2 of the match to the exon in exons1
    matches = [(0,0)]*len(exons)

    for i in xrange(exons1):
        e1 = exons1[i]
        bestMatch = -1
        for j in xrange(exons2):
            e2 = exons2[j]
            if e1[0] == e2[0] and e1[2] > e2[1] and e1[1] < e2[2]:
                dist = exonsMatch(e1,e2)
                if dist < bestMatch:
                  bestMatch = dist
                  closestExon = j
        if bestMatch >= 0:
            matches[i] = (closestExon, bestMatch)
        else:
            matches[i] = (-1, 0)