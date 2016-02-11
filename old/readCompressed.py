import time

class ReadCompressed:
    '''
        Support queries for compressed alignments file
    '''

    unspliced = None
    spliced = None
    chromosomes = None

    def __init__(self, unsplicedFile, splicedFile, chromosomes):
        self.unspliced = unsplicedFile
        self.spliced = splicedFile

        self.chromosomes = chromosomes

    def getCoverage(self, chrom, start=None, end=None):
        if not chrom in self.chromosomes:
            print 'Error! Chromosome name not recognized!'
            print 'Chromosomes: ' + ', '.join(self.chromosomes.keys())
            exit()

        if start == None or end == None:
            start = 0
            end = self.chromosomes[chrom]

        for k in sorted(self.chromosomes.keys()):
            if not k == chrom:
                start += self.chromosomes[k]
                end += self.chromosomes[k]
            else:
                break

        segmentLength = end-start
        coverage = [0.0] * segmentLength

        exons = []
        relevantExonsStart = 0
        relevantExonsEnd = 0
        with open(self.spliced, 'r', 1) as f:
            exons = None
            junc = None

            relevant = False

            startTime = time.time()
            for line in f:
                if exons == None:
                    # First line contains exons
                    exons = [int(e) for e in line.rstrip().split('\t')]

                    relevantExonsStart = 0
                    while exons[relevantExonsStart+1] < start:
                        relevantExonsStart += 1
                    relevantExonsEnd = relevantExonsStart
                    while exons[relevantExonsEnd] < end:
                        relevantExonsEnd += 1
                else:
                    # Remaining lines grouped by junction

                    if line[0] == '>':
                        # Start of a new junction
                        key = line[1:].rstrip().split('\t')
                        junctionExons = [int(e) for e in key[:-2]]

                        relevant = False
                        passedRange = True
                        for e in junctionExons:
                            if e < relevantExonsEnd:
                                passedRange = False
                                if e >= relevantExonsStart:
                                    relevant = True
                        if passedRange:
                            break

                        NH = 1.0 / float(key[-1])

                        count = 0
                        currExon = 0
                        currExonLen = exons[junctionExons[0]+1] - exons[junctionExons[0]]
                        offsetInExon = 0

                    elif relevant:
                        if count < 3:
                            count += 1
                        else:
                            # Rest of lines contain run-length-encoded coverage vector
                            row = line.rstrip().split('\t')
                            if len(row) == 1:
                                length = 1
                            else:
                                length = int(row[1])

                            prevOffset = offsetInExon+exons[junctionExons[currExon]]
                            offsetInExon += length
                            while offsetInExon > currExonLen:
                                if junctionExons[currExon] >= relevantExonsStart and junctionExons[currExon] < relevantExonsEnd:
                                    val = float(row[0]) * NH
                                    for i in xrange(max(prevOffset, exons[junctionExons[currExon]]), exons[junctionExons[currExon]+1]):
                                        if i >= start and i < end:
                                            coverage[i-start] += val

                                currExon += 1
                                offsetInExon -= currExonLen
                                currExonLen = exons[junctionExons[currExon]+1] - exons[junctionExons[currExon]]

                            if junctionExons[currExon] >= relevantExonsStart and junctionExons[currExon] < relevantExonsEnd:
                                val = float(row[0]) * NH
                                for i in xrange(max(prevOffset, exons[junctionExons[currExon]]), exons[junctionExons[currExon]]+offsetInExon):
                                    if i >= start and i < end:
                                        coverage[i-start] += val
            endTime = time.time()
            print 'Spliced: %0.3fs' % (endTime-startTime)
                

        with open(self.unspliced, 'r', 1) as f2:
            NH = 1.0
            RLE = False

            startTime = time.time()
            for line in f2:
                #line = line.rstrip()
                #if len(line) == 1:
                #    continue

                if not RLE:
                    if line[0] == '#':
                        NH = 1.0 / float(line.rstrip()[1:])
                        RLE = True
                        offset = 0
                else:
                    if line[0] == '>':
                        RLE = False
                    else:
                        row = line.rstrip().split()
                        if len(row) == 1:
                            length = 1
                        else:
                            length = int(row[1])

                        newOffset = offset+length
                        #if offset < end and newOffset > start:
                        if newOffset > start:
                            val = float(row[0]) * NH
                            for i in xrange(max(offset-start, 0), min(newOffset-start, segmentLength)):
                                coverage[i] += val
                        offset = newOffset

                        if offset >= end:
                            RLE = False
        endTime = time.time()
        print 'Unpliced: %0.3fs' % (endTime-startTime)

        return coverage

    def getGenes(self, chrom, start=None, end=None, overlapRadius=50):
        '''
            Return all genes within the given range in the form [(x_0,y_0), (x_1,y_1), ...] with a tuple indicating the bounds of each gene
        '''

        if not chrom in self.chromosomes:
            print 'Error! Chromosome name not recognized!'
            print 'Chromosomes: ' + ', '.join(self.chromosomes.keys())
            exit()

        if start == None or end == None:
            start = 0
            end = self.chromosomes[chrom]

        chromOffset = 0
        for k in sorted(self.chromosomes.keys()):
            if not k == chrom:
                chromOffset += self.chromosomes[k]
            else:
                break
        start += chromOffset
        end += chromOffset

        # Get chunks from spliced reads
        splicedChunks = []
        with open(self.spliced, 'r') as f:
            exons = None
            relevant = False
            for line in f:
                if exons == None:
                    # First line contains exons
                    exons = [int(e) for e in line.rstrip().split('\t')]

                    relevantExonsStart = 0
                    while exons[relevantExonsStart+1] < start:
                        relevantExonsStart += 1
                    relevantExonsEnd = relevantExonsStart
                    while exons[relevantExonsEnd] < end:
                        relevantExonsEnd += 1

                    startOffset = 0
                    endOffset = 0
                else:
                    # Remaining lines grouped by junction

                    if line[0] == '>':
                        # finish processing previous junction
                        if relevant:
                            splicedChunks.append((exons[junctionExons[0]]+startOffset-chromOffset, exons[junctionExons[-1]+1]-endOffset-chromOffset))

                        # Start of a new junction
                        key = line[1:].rstrip().split('\t')
                        junctionExons = [int(e) for e in key[:-2]]

                        relevant = True
                        for e in junctionExons:
                            if e < relevantExonsStart or e >= relevantExonsEnd:
                                relevant = False

                        count = 0

                    elif relevant:
                        if count < 3:
                            count += 1
                        elif count == 3:
                            # First line of run-length-encoded coverage vector
                            row = line.rstrip().split('\t')

                            val = int(row[0])

                            if val == 0:
                                if len(row) == 1:
                                    startOffset = 1
                                else:
                                    startOffset = int(row[1])
                            else:
                                startOffset = 0

                            count += 1
                        else:
                            # All other lines are run-length-encoded coverage vector
                            row = line.rstrip().split('\t')

                            val = int(row[0])

                            if val == 0:
                                if len(row) == 1:
                                    endOffset = 1
                                else:
                                    endOffset = int(row[1])
                            else:
                                endOffset = 0


        # Get chunks from unspliced reads
        unsplicedChunks = []
        with open(self.unspliced, 'r') as f:
            RLE = False

            for line in f:
                line = line.rstrip()
                if len(line) == 0:
                    continue

                if line[0] == '#':
                    RLE = True
                    chunk = False
                    pos = 0
                elif line[0] == '>':
                    RLE = False
                elif RLE:
                    row = line.rstrip().split('\t')
                    val = int(row[0])
                    if len(row) == 1:
                        length = 1
                    else:
                        length = int(row[1])

                    if pos >= start:
                        if (not chunk) and val > 0:
                            chunk = True
                            startPos = pos
                        elif chunk and val == 0:
                            unsplicedChunks.append((startPos-chromOffset, pos-chromOffset))
                            chunk = False

                    pos += length
                    if pos >= end:
                        RLE = False

        # Combine chunks
        # TODO: Check if these are already sorted
        splicedChunks.sort()
        unsplicedChunks.sort()

        if len(splicedChunks) == 0:
            return unsplicedChunks
        elif len(unsplicedChunks) == 0:
            return splicedChunks

        i = 0
        j = 0
        genes = []
        if splicedChunks[i] < unsplicedChunks[j]:
            genes = [splicedChunks[i]]
            i += 1
        else:
            genes = [unsplicedChunks[j]]
            j += 1

        while i < len(splicedChunks) and j < len(unsplicedChunks):
            if splicedChunks[i] < unsplicedChunks[j]:
                if splicedChunks[i][0] < genes[-1][1] + overlapRadius:
                    genes[-1] = (genes[-1][0], splicedChunks[i][1])
                else:
                    genes.append(splicedChunks[i])
                i += 1
            else:
                if unsplicedChunks[j][0] < genes[-1][1] + overlapRadius:
                    genes[-1] = (genes[-1][0], unsplicedChunks[j][1])
                else:
                    genes.append(unsplicedChunks[j])
                j += 1
        while i < len(splicedChunks):
            if splicedChunks[i][0] < genes[-1][1] + overlapRadius:
                genes[-1] = (genes[-1][0], splicedChunks[i][1])
            else:
                genes.append(splicedChunks[i])
            i += 1
        while j < len(unsplicedChunks):
            if unsplicedChunks[j][0] < genes[-1][1] + overlapRadius:
                genes[-1] = (genes[-1][0], unsplicedChunks[j][1])
            else:
                genes.append(unsplicedChunks[j])
            j += 1

        return genes



    def parseCigar(self, cigar, offset):
        ''' Parse the cigar string starting at the given index of the genome
            Returns a list of offsets for each exonic region of the read [(start1, end1), (start2, end2), ...]
        '''
        exons = []
        newExon = True

        # Parse cigar string
        match = re.search("\D", cigar)
        while match:
            index = match.start()
            length = int(''.join(cigar[:index]))

            if cigar[index] == 'N':
                # Separates contiguous exons, so set boolean to start a new one
                newExon = True
            elif cigar[index] == 'M':
                # If in the middle of a contiguous exon, append the length to it, otherwise start a new exon
                if newExon:
                    exons.append([offset, offset+length])
                    newExon = False
                else:
                    exons[-1][1] += length
            elif cigar[index] == 'D':
                # If in the middle of a contiguous exon, append the deleted length to it
                if not newExon:
                    exons[-1][1] += length

            offset += length
            cigar = cigar[index+1:]
            match = re.search("\D", cigar)

        return exons