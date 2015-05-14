#! /usr/bin/env python
import sys


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

            if row[1] == 'protein_coding' and row[2] == 'exon' and transcriptId in transcriptsTruth:
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

    print(sorted(transcriptCovs, reverse=True)[:10])

    #print "%d transcripts in file 1" % transcriptsTrueCount
    #print "%d transcripts in file 2" % transcriptsPredictedCount

    recall = float(totalScore1) / float(transcriptsTrueCount)
    print 'TP = ' + str(totalScore1)
    print 'T  = ' + str(transcriptsTrueCount)
    print 'Recall    = TP/T = ' + str(recall)


    precision = float(totalScore2) / float(transcriptsPredictedCount)
    print 'TP = ' + str(totalScore2)
    print 'P  = ' + str(transcriptsPredictedCount)
    print 'Precision = TP/P = ' + str(precision)

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


