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
        #if (not self.chrom == transcript.chrom) or (self.start > transcript.end) or (self.end < transcript.start):
        if ((self.start > transcript.end) or (self.end < transcript.start)):
            return 0

        thisMatches = []
        scores = []
        for i, e1 in enumerate(self.exons):
            closest = 0
            closestScore = 0

            for j in range(len(transcript.exons)):
                e2 = transcript.exons[j]
                currScore = self.scoreExons(e1, e2, threshold)

                if currScore > closestScore:
                    #logging.info('scoreExons returned %0.2f' % currScore)
                    closestScore = currScore
                    closest = j

            thisMatches.append(closest)
            scores.append(closestScore)

        otherMatches = []
        for i in range(len(transcript.exons)):
            e1 = transcript.exons[i]
            closest = 0
            closestScore = 0

            for j in range(len(self.exons)):
                e2 = self.exons[j]
                currScore = self.scoreExons(e1, e2, threshold)

                if currScore > closestScore:
                    #logging.info('scoreExons returned %0.2f' % currScore)
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
        for i in range(len(thisMatches)):
            if otherMatches[thisMatches[i]] == i:
                totalScore += scores[i]
        for i in range(len(otherMatches)):
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
