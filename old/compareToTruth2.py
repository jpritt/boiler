#! /usr/bin/env python
import sys
import csv
import math
import re
import string
import numpy as np
import matplotlib.pyplot as plt

class Transcript:
    '''
    '''
    def __init__(self, chrom, start, end, cov, name):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.cov = cov
        self.exons = []
        self.name = name

    def scoreTranscript(self, transcript, threshold):
        ''' Compare this transcript to the given transcript and return a score representing their closeness
            0 = no match, 1 = perfect match
        '''
        if (not self.chrom == transcript.chrom) or (self.start > transcript.end) or (self.end < transcript.start):
            return 0

        #print transcript.name

        thisMatches = []
        scores = []
        for i in xrange(len(self.exons)):
            e1 = self.exons[i]
            closest = 0
            closestScore = 0

            for j in xrange(len(transcript.exons)):
                e2 = transcript.exons[j]
                currScore = self.scoreExons(e1, e2, threshold)

                if currScore > closestScore:
                    closestScore = currScore
                    closest = j

            thisMatches.append(closest)
            scores.append(closestScore)

        otherMatches = []
        for i in xrange(len(transcript.exons)):
            e1 = transcript.exons[i]
            closest = 0
            closestScore = 0

            for j in xrange(len(self.exons)):
                e2 = self.exons[j]
                currScore = self.scoreExons(e1, e2, threshold)

                if currScore > closestScore:
                    closestScore = currScore
                    closest = j

            otherMatches.append(closest)

        '''
        print 'Comparing %s and %s' % (self.name, transcript.name)
        print 'This matches:  ' + '\t'.join([str(c) for c in thisMatches])
        print 'Other matches: ' + '\t'.join([str(c) for c in otherMatches])
        print 'Scores:        ' + '\t'.join([str(c) for c in scores])
        '''

        totalScore = 0.0
        count = len(thisMatches)
        for i in xrange(len(thisMatches)):
            if otherMatches[thisMatches[i]] == i:
                totalScore += scores[i]
        for i in xrange(len(otherMatches)):
            if not thisMatches[otherMatches[i]] == i:
                count += 1
        
        '''
        print '%f / %f = %f' % (totalScore, count, totalScore/float(count))
        print ''
        '''

        return totalScore / float(count)

    def scoreExons(self, exon1, exon2, threshold):
        ''' Returns a value between 0 and 1, where 1 indicates that the 2 exons are a perfect match and 0 indicates no match
        '''

        dx = min(abs(exon1[0] - exon2[0]), threshold)
        dy = min(abs(exon1[1] - exon2[1]), threshold)

        return 1 - dx / (2.0*threshold) - dy / (2.0*threshold)


'''
Compare a Cufflinks GTF file to the .pro output by flux
'''

def compareGTFs(proFile, truthGTF, compGTF1, compGTF2):
    # file 1 is a .pro file output by flux
    transcriptsTruth = parsePro(proFile)

    with open(truthGTF, 'r') as tsv:
        for line in tsv:
            row = line.strip().split('\t')
            if len(row) < 5:
                continue

            transcriptIdIndex = row[8].find('transcript_id')
            transcriptIdStart = row[8].find('"', transcriptIdIndex) + 1
            transcriptIdEnd = row[8].find('"', transcriptIdStart)
            transcriptId = row[8][transcriptIdStart : transcriptIdEnd]

            #if row[2] == 'transcript':
            #    transcriptsTruth[transcriptId] = Transcript(row[0], int(row[3]), int(row[4]), transcriptCovs[transcriptId])

            if row[1] == 'protein_coding' and row[2] == 'exon' and transcriptId in transcriptsTruth:
                transcriptsTruth[transcriptId].exons.append( (int(row[3]), int(row[4])) )
    transcriptsTruth = transcriptsTruth.values()

    transcriptsComp1 = dict()
    with open(compGTF1, 'r') as tsv:
        for line in tsv:
            row = line.strip().split('\t')
            if len(row) < 5:
                continue

            covIndex = row[8].find('cov')
            covStart = row[8].find('"', covIndex) + 1
            covEnd = row[8].find('"', covStart)
            cov = float(row[8][covStart:covEnd])

            transcriptIdIndex = row[8].find('transcript_id')
            transcriptIdStart = row[8].find('"', transcriptIdIndex) + 1
            transcriptIdEnd = row[8].find('"', transcriptIdStart)
            transcriptId = row[8][transcriptIdStart : transcriptIdEnd]

            if row[2] == 'transcript':
                transcriptsComp1[transcriptId] = Transcript(row[0], int(row[3]), int(row[4]), cov, transcriptId)
            elif row[2] == 'exon':
                transcriptsComp1[transcriptId].exons.append( (int(row[3]), int(row[4])) )
    transcriptsComp1 = transcriptsComp1.values()

    transcriptsComp2 = dict()
    with open(compGTF2, 'r') as tsv:
        for line in tsv:
            row = line.strip().split('\t')
            if len(row) < 5:
                continue

            covIndex = row[8].find('cov')
            covStart = row[8].find('"', covIndex) + 1
            covEnd = row[8].find('"', covStart)
            cov = float(row[8][covStart:covEnd])

            transcriptIdIndex = row[8].find('transcript_id')
            transcriptIdStart = row[8].find('"', transcriptIdIndex) + 1
            transcriptIdEnd = row[8].find('"', transcriptIdStart)
            transcriptId = row[8][transcriptIdStart : transcriptIdEnd]

            if row[2] == 'transcript':
                transcriptsComp2[transcriptId] = Transcript(row[0], int(row[3]), int(row[4]), cov, transcriptId)
            elif row[2] == 'exon':
                transcriptsComp2[transcriptId].exons.append( (int(row[3]), int(row[4])) )
    transcriptsComp2 = transcriptsComp1.values()

    compareAll(transcriptsTruth, transcriptsComp1, transcriptsComp2)

def parsePro(filename):
    ''' Return a dictionary with transcript id (e.g. 0300689) pointing to coverage level
    '''
    threshold = 0.00005

    transcripts = dict()
    with open(filename, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            if len(row) < 8:
                continue

            tag = row[1]
            sep1 = row[0].find(':')
            sep2 = row[0].find('-', sep1)
            sep3 = row[0].find('W', sep2)
            chrom = row[0][:sep1]
            start = int(row[0][sep1+1:sep2])
            end = int(row[0][sep2+1:sep3])
            fraction = float(row[8])

            if fraction > threshold:
                transcripts[tag] = Transcript(chrom, start, end, fraction, tag)
    return transcripts


def compareAll(transcriptsTrue, transcriptsPredictedA, transcriptsPredictedB):
    ''' Compare 
    '''
    
    totalScore = 0.0
    threshold = 10

    # match each transcript from predicted set A to the truth
    for t in transcriptsPredictedA:
        closestScore = 0.0
        closestT = None

        for t2 in transcriptsTrue:
            if t.chrom == t2.chrom and t.start < t2.end and t.end > t2.start:
                score = t.scoreTranscript(t2, threshold)
                if score > closestScore:
                    closestScore = score
                    closestT = t2

        if not closestT == None
        


    transcriptsTrueCount = 0
    line = 0
    for t1 in transcriptsTrue:
        line += 1
        #if line % 1000 == 0:
        #    print "%d / %d transcripts done" % (line, len(transcriptsTrue))

#        if t1.cov < .000001:
#            continue


        closestScore = 0.0
        closestT = None
        
        for t2 in transcriptsPredicted:
            score = t1.scoreTranscript(t2, threshold)
            if score > closestScore:
                closestScore = score
                closestT = t2

        '''
        if closestT == None:
            print '%s\t%s\t\t%f' % (t1.name, closestT, closestScore)
        else:
            print '%s\t%s\t%f' % (t1.name, closestT.name, closestScore)
        '''

        if closestScore > 0:
            totalScore += closestScore
        transcriptsTrueCount += 1

    transcriptsPredictedCount = len(transcriptsPredicted)

    print "%d transcripts in file 1" % transcriptsTrueCount
    print "%d transcripts in file 2" % transcriptsPredictedCount

    recall = float(totalScore) / float(transcriptsTrueCount)
    print 'TP = ' + str(totalScore)
    print 'T  = ' + str(transcriptsTrueCount)
    print 'Recall    = TP/T = ' + str(recall)


    precision = float(totalScore) / float(transcriptsPredictedCount)
    print 'TP = ' + str(totalScore)
    print 'P  = ' + str(transcriptsPredictedCount)
    print 'Recall    = TP/T = ' + str(precision)


compareGTFs(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


