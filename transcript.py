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
            #if ((self.start > transcript.end) or (self.end < transcript.start)):
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

        totalScore = 0.0
        count = len(thisMatches)
        for i in range(len(thisMatches)):
            if otherMatches[thisMatches[i]] == i:
                totalScore += scores[i]
        for i in range(len(otherMatches)):
            if not thisMatches[otherMatches[i]] == i:
                count += 1

        return totalScore / float(count)

    def scoreTranscript2(self, transcript, threshold):
        if (not self.chrom == transcript.chrom) or ((self.start > transcript.end) or (self.end < transcript.start)):
            return 0

        if not len(self.exons) == len(transcript.exons):
            return 0

        score = 1
        for i in range(len(self.exons)):
            s = self.scoreExons(self.exons[i], transcript.exons[i], threshold)
            if s == 0:
                return 0
            else:
                score *= s

        return score

    
    def scoreTranscript3(self, transcript, threshold):
        #Same as scoreTranscript() above but calculate score as proportion of overlap
        if (not self.chrom == transcript.chrom) or ((self.start > transcript.end) or (self.end < transcript.start)):
            return 0

        overlap = 0
        i = 0
        for e1 in self.exons:
            while i < len(transcript.exons) and transcript.exons[i][1] <= e1[0]:
                i += 1

            if i == len(transcript.exons):
                break

            while i < len(transcript.exons) and transcript.exons[i][0] < e1[1]:
                overlap += min(transcript.exons[i][1], e1[1]) - max(transcript.exons[i][0], e1[0])
                if transcript.exons[i][0] < e1[1]:
                    i += 1

        len1 = 0
        for e1 in transcript.exons:
            len1 += e1[1] - e1[0]

        len2 = 0
        for e2 in transcript.exons:
            len2 += e2[1] - e2[0]

        if overlap > len1 or overlap > len2:
            print(self.exons)
            print(transcript.exons)
            print('Length 1 = %d' % len1)
            print('Length 2 = %d' % len2)
            print('Overlap = %d' % overlap)
            exit()

        return float(overlap * overlap) / float(len1 * len2)

    def scoreExons(self, exon1, exon2, threshold):
        ''' Returns a value between 0 and 1, where 1 indicates that the 2 exons are a perfect match and 0 indicates no match
        '''

        if threshold == 0:
            if exon1[0] == exon2[0] and exon1[1] == exon2[1]:
                return 1
            else:
                return 0
        else:
            dx = min(abs(exon1[0] - exon2[0]), threshold)
            dy = min(abs(exon1[1] - exon2[1]), threshold)

            return 1 - dx / (2.0*threshold) - dy / (2.0*threshold)
