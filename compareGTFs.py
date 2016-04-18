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

            covIndex = row[8].find('cov')
            covStart = row[8].find('"', covIndex) + 1
            covEnd = row[8].find('"', covStart)
            cov = float(row[8][covStart:covEnd])

            transcriptIdIndex = row[8].find('transcript_id')
            transcriptIdStart = row[8].find('"', transcriptIdIndex) + 1
            transcriptIdEnd = row[8].find('"', transcriptIdStart)
            transcriptId = row[8][transcriptIdStart : transcriptIdEnd]

            if row[2] == 'transcript':
                transcriptsTruth[transcriptId] = Transcript(row[0], int(row[3]), int(row[4]), cov, transcriptId)
            elif row[2] == 'exon':
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
            transcripts[tag] = Transcript(chrom, start, end, fraction, tag)
    return transcripts


def compareAll(transcriptsTrue, transcriptsPredicted):
    ''' Compare 
    '''
    
    totalScore = 0.0
    threshold = 10

    transcriptsTrueCount = 0
    line = 0
    for t1 in transcriptsTrue:
        line += 1

        if t1.cov < .000001:
            continue


        closestScore = 0
        closestT = None
        
        for t2 in transcriptsPredicted:
            score = t1.scoreTranscript(t2, threshold)
            if score > closestScore:
                closestScore = score
                closestT = t2

        if weightByCov:
            if closestScore > 0:
                totalScore += closestScore * t1.cov
            transcriptsTrueCount += t1.cov
        else:
            if closestScore > 0:
                totalScore += closestScore
            transcriptsTrueCount += 1
    recall = float(totalScore) / float(transcriptsTrueCount)
    print 'Recall    = TP/T = ' + str(recall)
    
    line = 0
    totalScore = 0
    transcriptsPredictedCount = 0
    for t1 in transcriptsPredicted:
        line += 1

        if t1.cov < .000001:
            continue


        closestScore = 0
        closestT = None
        
        for t2 in transcriptsTrue:
            score = t1.scoreTranscript(t2, threshold)
            if score > closestScore:
                closestScore = score
                closestT = t2

        if weightByCov:
            if closestScore > 0:
                totalScore += closestScore * t1.cov
            transcriptsPredictedCount += t1.cov
        else:
            if closestScore > 0:
                totalScore += closestScore
            transcriptsPredictedCount += 1
    precision = float(totalScore) / float(transcriptsPredictedCount)
    print 'Precision    = TP/P = ' + str(precision)

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

weightByCov = False
if sys.argv[3] == '1':
    weightByCov = True
compareGTFs(sys.argv[1], sys.argv[2])



# def compareGTFs(file1, file2):
#     lines1 = []
#     with open(file1, 'r') as tsv:
#         for line in tsv:
#             row = line.strip().split('\t')
#             if len(row) >= 5:
#                 covIndex = row[8].find('cov')
#                 covStart = row[8].find('"', covIndex) + 1
#                 covEnd = row[8].find('"', covStart)

#                 lines1.append(row[:3] + [int(row[3])] + [int(row[4])])# + [float(row[8][covStart:covEnd])])

#     lines2 = []
#     with open(file2, 'r') as tsv:
#         for line in tsv:
#             row = line.strip().split('\t')
#             if len(row) >= 5:
#                 lines2.append(row[:3] + [int(row[3])] + [int(row[4])])

#     compareAll(lines1, lines2)

# def compareAll(lines1, lines2):
#     ''' Compare 
#     '''
#     numLines1 = len(lines1)
#     numLines2 = len(lines2)
    
#     testedLines1 = 0.0
#     matching = 0.0

#     #covThreshold = 0

#     print "%d lines in file 1" % numLines1
#     print "%d lines in file 2" % numLines2

#     index = 0
#     for l1 in lines1:
#         index += 1
#         if index % 1000 == 0:
#             print "%d / %d lines done" % (index, numLines1)

#         #if l1[5] <= covThreshold:
#         #    continue

#         testedLines1 += 1
#         closestScore = 0.0
        
#         for l2 in lines2:
#             if l1[0] == l2[0] and l1[2] == l2[2]:
#                 if l1[2] == 'transcript':
#                     threshold = 100
#                 else:
#                     threshold = 10
#                 score = exonsMatch((l1[3], l1[4]), (l2[3], l2[4]), threshold)
#                 if score > closestScore:
#                     closestScore = score
#                     closestMatch = l2
#         if closestScore > 0:
#             matching += closestScore
#             #matching += l1[5]

#     #print str(numTranscripts) + ' transcripts'
#     print 'Using %d of %d lines' % (testedLines1, numLines1)
#     print 'TP = ' + str(matching)
#     print 'T  = ' + str(testedLines1)
#     #print 'P  = ' + str(numLines2)

#     #print 'Precision = TP/P = ' + str(matching / numLines2)
#     print 'Recall    = TP/T = ' + str(matching / testedLines1)


# def exonsMatch(exon1, exon2, threshold):
#     ''' Returns a value between 0 and 1, where 1 indicates that the 2 exons are a perfect match and 0 indicates no match
#     '''

#     dx = min(abs(exon1[0] - exon2[0]), threshold)
#     dy = min(abs(exon1[1] - exon2[1]), threshold)

#     return 1 - dx / (2.0*threshold) - dy / (2.0*threshold)


# compareGTFs(sys.argv[1], sys.argv[2])
