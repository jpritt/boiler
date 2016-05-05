#! /usr/bin/env python
import sys
from transcript import Transcript

'''
Compare a Cufflinks GTF file to the .pro output by flux
'''

def compareGTFs(truthGTF, compGTF):
    transcriptsTruth = dict()
    with open(truthGTF, 'r') as tsv:
        for line in tsv:
            row = line.strip().split('\t')
            if len(row) < 5:
                continue

            transcriptIdIndex = row[8].find('transcript_id')
            transcriptIdStart = row[8].find('"', transcriptIdIndex) + 1
            transcriptIdEnd = row[8].find('"', transcriptIdStart)
            transcriptId = row[8][transcriptIdStart : transcriptIdEnd]

            print(row[2])

            if row[2] == 'transcript':
                print('Found transcript ' + str(transcriptId))
                transcriptsTruth[transcriptId] = Transcript(row[0], int(row[3]), int(row[4]), 1, transcriptId)
            elif row[2] == 'exon' and transcriptId in transcriptsTruth:
                transcriptsTruth[transcriptId].exons.append( (int(row[3]), int(row[4])) )

    transcriptsTruth = transcriptsTruth.values()

    transcriptsComp = dict()
    with open(compGTF, 'r') as tsv:
        for line in tsv:
            row = line.strip().split('\t')
            if len(row) < 5:
                continue

            transcriptIdIndex = row[8].find('transcript_id')
            transcriptIdStart = row[8].find('"', transcriptIdIndex) + 1
            transcriptIdEnd = row[8].find('"', transcriptIdStart)
            transcriptId = row[8][transcriptIdStart : transcriptIdEnd]

            if row[2] == 'transcript':
                transcriptsComp[transcriptId] = Transcript(row[0], int(row[3]), int(row[4]), 1, transcriptId)
            elif row[2] == 'exon' and transcriptId in transcriptsComp:
                transcriptsComp[transcriptId].exons.append( (int(row[3]), int(row[4])) )

    transcriptsComp = transcriptsComp.values()

    compareAll(transcriptsTruth, transcriptsComp)

def compareAll(transcriptsTrue, transcriptsPredicted):
    ''' Compare 
    '''
    
    totalScore = 0.0
    threshold = 10

    transcriptsTrueCount = 0
    line = 0
    for t1 in transcriptsTrue:
        line += 1

        closestScore = 0
        closestT = None
        
        for t2 in transcriptsPredicted:
            score = t1.scoreTranscript(t2, threshold)
            if score > closestScore:
                closestScore = score
                closestT = t2

        if closestScore > 0:
            totalScore += closestScore
        transcriptsTrueCount += 1
    recall = float(totalScore) / float(transcriptsTrueCount)
    print('Recall    = TP/T = ' + str(recall))
    
    line = 0
    totalScore = 0
    transcriptsPredictedCount = 0
    for t1 in transcriptsPredicted:
        line += 1

        closestScore = 0
        closestT = None
        
        for t2 in transcriptsTrue:
            score = t1.scoreTranscript(t2, threshold)
            if score > closestScore:
                closestScore = score
                closestT = t2

        if closestScore > 0:
            totalScore += closestScore
        transcriptsPredictedCount += 1
    precision = float(totalScore) / float(transcriptsPredictedCount)
    print('Precision    = TP/P = ' + str(precision))

compareGTFs(sys.argv[1], sys.argv[2])

