#! /usr/bin/env python
import sys
from transcript import Transcript

'''
Compare a Cufflinks GTF file to the .pro output by flux
'''

def compareGTFs(proFile, truthGTF, compGTF):
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

            #if row[1] == 'protein_coding' and row[2] == 'exon' and transcriptId in transcriptsTruth:
            if row[2] == 'exon' and transcriptId in transcriptsTruth:
                transcriptsTruth[transcriptId].exons.append( (int(row[3]), int(row[4])) )
    transcriptsTruth = transcriptsTruth.values()

    transcriptsComp = dict()
    with open(compGTF, 'r') as tsv:
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
                transcriptsComp[transcriptId] = Transcript(row[0], int(row[3]), int(row[4]), cov, transcriptId)
            elif row[2] == 'exon':
                transcriptsComp[transcriptId].exons.append( (int(row[3]), int(row[4])) )

    transcriptsComp = transcriptsComp.values()

    compareAll(transcriptsTruth, transcriptsComp)

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
            #cov = float(row[7])
            #if cov > 0:
            #    fraction = float(row[11]) / cov

            if fraction > threshold:
                transcripts[tag] = Transcript(chrom, start, end, fraction, tag)
    return transcripts


def compareAll(transcriptsTrue, transcriptsPredicted):
    ''' Compare 
    '''
    
    threshold = 10

    totalScore1 = 0.0
    transcriptsTrueCount = 0
    line = 0
    for t1 in transcriptsTrue:
        line += 1

        closestScore = 0.0
        closestT = None
        
        for t2 in transcriptsPredicted:
            score = t1.scoreTranscript(t2, threshold)
            if score > closestScore:
                closestScore = score
                closestT = t2

        if closestScore > 0:
            totalScore1 += closestScore * t1.cov
            #totalScore1 += closestScore
        transcriptsTrueCount += t1.cov
        #transcriptsTrueCount += 1

        

    #transcriptsPredictedCount = len(transcriptsPredicted)
    totalScore2 = 0.0
    transcriptsPredictedCount = 0
    line = 0


    transcriptCovs = []

    for t1 in transcriptsPredicted:
        transcriptCovs.append((t1.cov, t1.name))
        line += 1

        closestScore = 0.0
        closestT = None
        
        for t2 in transcriptsTrue:
            score = t1.scoreTranscript(t2, threshold)
            if score > closestScore:
                closestScore = score
                closestT = t2

        if closestScore > 0:
            totalScore2 += closestScore * t1.cov
            #totalScore2 += closestScore
        transcriptsPredictedCount += t1.cov
        #transcriptsPredictedCount += 1

    #print(sorted(transcriptCovs, reverse=True)[:10])

    #print "%d transcripts in file 1" % transcriptsTrueCount
    #print "%d transcripts in file 2" % transcriptsPredictedCount

    recall = float(totalScore1) / float(transcriptsTrueCount)
    #print('TP = ' + str(totalScore1))
    #print('T  = ' + str(transcriptsTrueCount))
    print('Recall    = ' + str(recall))


    precision = float(totalScore2) / float(transcriptsPredictedCount)
    #print('TP = ' + str(totalScore2))
    #print('P  = ' + str(transcriptsPredictedCount))
    print('Precision = ' + str(precision))

    '''
    scatterX = []
    scatterY = []
    for i in xrange(len(matchingLines)):
        if matchingLines[i] >= 0:
            l1 = lines1[i]
            l2 = lines2[matchingLines[i]]
            matching += exonsMatch((l1[3], l1[4]), (l2[3], l2[4]), threshold)
            scatterX += [l1[5]]
            scatterY += [l2[5]]

    plt.figure(1)
    plt.scatter(scatterX, scatterY)
    plt.savefig('scatter.png')
    '''

compareGTFs(sys.argv[1], sys.argv[2], sys.argv[3])


