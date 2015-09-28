#! /usr/bin/env python3
import alignments
import read
import pairedread
import junction
import time
import math
import binaryIO
import bisect

class Expander:
    aligned = None
    sectionLen = 100000    

    # 0 - zlib
    # 1 - lzma
    # 2 - bz2
    compressMethod = 0

    def __init__(self):
        self.debug = False
        if self.compressMethod == 0:
            self.zlib = __import__('zlib')
        elif self.compressMethod == 1:
            self.lzma = __import__('lzma')
        elif self.compressMethod == 2:
            self.bz2 = __import__('bz2')


    def expand(self, compressedFilename, uncompressedFilename, binary=False, debug=False):
        ''' Expand both spliced and unspliced alignments
        '''

        print('Expanding')

        self.debug = debug

        self.aligned = None

        if binary:
            with open(compressedFilename, 'rb') as f:
                if self.debug:
                    start = time.time()

                chroms = binaryIO.readChroms(f)
                self.aligned = alignments.Alignments(chroms, self.debug)

                self.expandByCluster(f, uncompressedFilename)

                if self.debug:
                    end = time.time()
                    print('Decompression time: %0.3fs' % (end-start))

        else:
            print('Non-binary expanding not supported')
            exit()

    '''
    def expandSplicedBinary(self, f, exonBytes):
        readLenBytes = binaryIO.readVal(f, 1)

        # Read junctions in compressed chunks
        juncId = 0
        for chunkLen in self.junctionChunkLens:
            s = self.expandString(f.read(chunkLen))
            startPos = 0

            while startPos < len(s):
                key = self.sortedJuncs[juncId]
                juncId += 1
                exons = key[:-2]

                length = 0
                boundaries = []
                for i in range(len(exons)):
                    e = exons[i]
                    subexon_length = self.aligned.exons[e+1] - self.aligned.exons[e]

                    length += subexon_length
                    if i < len(exons)-1:
                        if i == 0:
                            boundaries.append(subexon_length)
                        else:
                            boundaries.append(boundaries[-1] + subexon_length)

                debug = False

                # Read the rest of the junction information
                junc, startPos = binaryIO.readJunction(s, junction.Junction(exons, length, boundaries), readLenBytes, startPos)
                junc.xs = key[-2]
                junc.NH = key[-1]
                junc.coverage = self.RLEtoVector(junc.coverage)

                unpaired, paired = self.aligned.findReads(junc.unpairedLens, junc.pairedLens, junc.lensLeft, junc.lensRight, junc.coverage, junc.boundaries, debug)

                juncBounds = []
                for j in junc.exons:
                    juncBounds.append([self.aligned.exons[j], self.aligned.exons[j+1]])

                # Offset of start of junction from beginning of chromosome
                juncOffset = juncBounds[0][0]

                # marks indices in junction coverage vector where exons are
                mapping = [0]
                # offsets of start of each exon in junction
                offsets = [0]
                for j in range(1,len(juncBounds)):
                    mapping.append(mapping[-1] + juncBounds[j-1][1] - juncBounds[j-1][0])
                    offsets.append(juncBounds[j][0] - juncBounds[0][0])

                for r in unpaired:
                    start = r[0]
                    i = 0
                    while i < (len(mapping)-1) and mapping[i+1] <= start:
                        i += 1
                    start += offsets[i] - mapping[i]

                    end = r[1]
                    j = 0
                    while j < (len(mapping)-1) and mapping[j+1] < end:
                        j += 1
                    end += offsets[j] - mapping[j]
                    
                    readExons = []
                    if i == j:
                        readExons.append( [start+juncOffset, end+juncOffset] )
                    else:
                        readExons.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                        for x in range(i+1,j):
                            readExons.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                        readExons.append( [offsets[j]+juncOffset, end+juncOffset] ) 

                    self.aligned.spliced.append(read.Read(self.aligned.getChromosome(readExons[0][0]), readExons, junc.xs, junc.NH))

                for p in paired:
                    start = p[0][0]
                    i = 0
                    while i < (len(mapping)-1) and mapping[i+1] <= start:
                        i += 1
                    start += offsets[i] - mapping[i]

                    end = p[0][1]
                    j = 0
                    while j < (len(mapping)-1) and mapping[j+1] < end:
                        j += 1
                    end += offsets[j] - mapping[j]
                    
                    readExonsA = []
                    if i == j:
                        readExonsA.append( [start+juncOffset, end+juncOffset] )
                    else:
                        readExonsA.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                        for x in range(i+1,j):
                            readExonsA.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                        readExonsA.append( [offsets[j]+juncOffset, end+juncOffset] )

                    start = p[1][0]
                    i = 0
                    while i < (len(mapping)-1) and mapping[i+1] <= start:
                        i += 1
                    start += offsets[i] - mapping[i]

                    end = p[1][1]
                    j = 0
                    while j < (len(mapping)-1) and mapping[j+1] < end:
                        j += 1
                    end += offsets[j] - mapping[j]
                    
                    readExonsB = []
                    if i == j:
                        readExonsB.append( [start+juncOffset, end+juncOffset] )
                    else:
                        readExonsB.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                        for x in range(i+1,j):
                            readExonsB.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                        readExonsB.append( [offsets[j]+juncOffset, end+juncOffset] )

                    self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(readExonsA[0][0]), readExonsA,  \
                                                             self.aligned.getChromosome(readExonsB[0][0]), readExonsB, junc.xs, junc.NH))

    def expandUnsplicedBinary(self, f):

        NHs = binaryIO.readVal(f, 2)
        self.sectionLen = binaryIO.readVal(f, 4)
        fragLenBytes = binaryIO.readVal(f, 1)
        readLenBytes = binaryIO.readVal(f, 1)

        for _ in range(NHs):
            NH = binaryIO.readVal(f, 2)
            coverage = []
            for chunk in range(len(self.unsplicedIndex[NH])):
                covLength = self.unsplicedIndex[NH][chunk]
                if covLength > 0:
                    s = f.read(covLength)
                    cov,_ = binaryIO.readCov(self.expandString(s))
                    coverage += self.RLEtoVector(cov)
                else:
                    if chunk == len(self.unsplicedIndex[NH])-1:
                        coverage += [0] * (self.aligned.exons[-1] - chunk*self.sectionLen)
                    else:
                        coverage += [0] * self.sectionLen

            for e in range(len(self.unsplicedExonsIndex[NH])):
                lenDists = self.expandString(f.read(self.unsplicedExonsIndex[NH][e]))
                startPos = 0
                i = 0

                while startPos < len(lenDists):
                    i += 1

                    unpairedLens, startPos = binaryIO.readLens(lenDists, readLenBytes, startPos)
                    pairs, startPos = binaryIO.readLens(lenDists, fragLenBytes, startPos)
                    lensLeft = dict()
                    lensRight = dict()
                    if len(pairs) > 0:
                        lensLeft, startPos = binaryIO.readLens(lenDists, readLenBytes, startPos)
                        if len(lensLeft) > 0:
                            lensRight, startPos = binaryIO.readLens(lenDists, readLenBytes, startPos)


                    if len(unpairedLens) > 0 or len(pairs) > 0:
                        exonStart = self.aligned.exons[e*self.exonChunkSize + i - 1]
                        exonEnd = self.aligned.exons[e*self.exonChunkSize + i]

                        debug = False

                        unpaired, paired = self.aligned.findReads(unpairedLens, pairs, lensLeft, lensRight, coverage[exonStart:exonEnd], None, debug)

                        for r in unpaired:
                            self.aligned.unspliced.append(read.Read(self.aligned.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))
                        for p in paired:
                            self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                                         self.aligned.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH))
    '''

    def expandUnsplicedBinary(self, f, length):
        ''' Expand a file containing compressed unspliced alignments
        '''

        s = self.expandString(f.read(length))

        start = 0
        num_NH, start = binaryIO.binaryToVal(s, 2, start)
        fragLenBytes, start = binaryIO.binaryToVal(s, 1, start)
        readLenBytes, start = binaryIO.binaryToVal(s, 1, start)

        for i in range(num_NH):
            NH, start = binaryIO.binaryToVal(s, 2, start)

            cov, start = binaryIO.readCov(s, start)
            cov = self.RLEtoVector(cov)

            for i in range(len(self.aligned.exons)-1):

                unpairedLens, start = binaryIO.readLens(s, readLenBytes, start)
                pairedLens, start = binaryIO.readLens(s, fragLenBytes, start)
                if len(pairedLens) > 0:
                    lensLeft, start = binaryIO.readLens(s, readLenBytes, start)
                    if len(lensLeft) > 0:
                        lensRight, start = binaryIO.readLens(s, readLenBytes, start)
                    else:
                        lensRight = dict()
                else:
                    lensLeft = dict()
                    lensRight = dict()

                if len(unpairedLens) > 0 or len(pairedLens) > 0:
                    exonStart = self.aligned.exons[i]
                    exonEnd = self.aligned.exons[i+1]

                    #print(cov[exonStart-self.aligned.exons[0]:exonEnd-self.aligned.exons[0]])
                    #print(unpairedLens)
                    #print(pairedLens)
                    #print('')

                    unpaired, paired = self.aligned.findReads(unpairedLens, pairedLens, lensLeft, lensRight, cov[exonStart-self.aligned.exons[0]:exonEnd-self.aligned.exons[0]], None, debug=False)

                    for r in unpaired:
                        self.aligned.unspliced.append(read.Read(self.aligned.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))
                    for p in paired:
                        self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                                     self.aligned.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH))

    def expandSplicedBinary(self, f, length):
        ''' Expand a file containing compressed spliced alignments
        '''

        s = self.expandString(f.read(length))

        startPos = 0
        readLenBytes, startPos = binaryIO.binaryToVal(s, 1, startPos)

        sorted_junctions, exonBytes, startPos = binaryIO.readJunctionsList(s, startPos)

        l = len(s)
        for key in sorted_junctions:
            exons = key[:-2]

            length = 0
            boundaries = []
            for i in range(len(exons)):
                e = exons[i]
                subexon_length = self.aligned.exons[e+1] - self.aligned.exons[e]

                length += subexon_length
                if i < len(exons)-1:
                    if i == 0:
                        boundaries.append(subexon_length)
                    else:
                        boundaries.append(boundaries[-1] + subexon_length)

            debug = False

            # Read the rest of the junction information
            junc, startPos = binaryIO.readJunction(s, junction.Junction(exons, length, boundaries), readLenBytes, startPos)
            junc.xs = key[-2]
            junc.NH = key[-1]

            #print(self.aligned.exons)
            #print(exons)
            #print(junc.coverage)
            #print(boundaries)
            junc.coverage = self.RLEtoVector(junc.coverage)

            #print(junc.unpairedLens)
            #print(junc.pairedLens)
            #print('')

            unpaired, paired = self.aligned.findReads(junc.unpairedLens, junc.pairedLens, junc.lensLeft, junc.lensRight, junc.coverage, junc.boundaries, debug)

            juncBounds = []
            for j in junc.exons:
                juncBounds.append([self.aligned.exons[j], self.aligned.exons[j+1]])

            # Offset of start of junction from beginning of chromosome
            juncOffset = juncBounds[0][0]

            # marks indices in junction coverage vector where exons are
            mapping = [0]
            # offsets of start of each exon in junction
            offsets = [0]
            for j in range(1,len(juncBounds)):
                mapping.append(mapping[-1] + juncBounds[j-1][1] - juncBounds[j-1][0])
                offsets.append(juncBounds[j][0] - juncBounds[0][0])

            for r in unpaired:
                start = r[0]
                i = 0
                while i < (len(mapping)-1) and mapping[i+1] <= start:
                    i += 1
                start += offsets[i] - mapping[i]

                end = r[1]
                j = 0
                while j < (len(mapping)-1) and mapping[j+1] < end:
                    j += 1
                end += offsets[j] - mapping[j]

                readExons = []
                if i == j:
                    readExons.append( [start+juncOffset, end+juncOffset] )
                else:
                    readExons.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                    for x in range(i+1,j):
                        readExons.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                    readExons.append( [offsets[j]+juncOffset, end+juncOffset] )

                self.aligned.spliced.append(read.Read(self.aligned.getChromosome(readExons[0][0]), readExons, junc.xs, junc.NH))

            for p in paired:
                start = p[0][0]
                i = 0
                while i < (len(mapping)-1) and mapping[i+1] <= start:
                    i += 1
                start += offsets[i] - mapping[i]

                end = p[0][1]
                j = 0
                while j < (len(mapping)-1) and mapping[j+1] < end:
                    j += 1
                end += offsets[j] - mapping[j]

                readExonsA = []
                if i == j:
                    readExonsA.append( [start+juncOffset, end+juncOffset] )
                else:
                    readExonsA.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                    for x in range(i+1,j):
                        readExonsA.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                    readExonsA.append( [offsets[j]+juncOffset, end+juncOffset] )

                start = p[1][0]
                i = 0
                while i < (len(mapping)-1) and mapping[i+1] <= start:
                    i += 1
                start += offsets[i] - mapping[i]

                end = p[1][1]
                j = 0
                while j < (len(mapping)-1) and mapping[j+1] < end:
                    j += 1
                end += offsets[j] - mapping[j]

                readExonsB = []
                if i == j:
                    readExonsB.append( [start+juncOffset, end+juncOffset] )
                else:
                    readExonsB.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                    for x in range(i+1,j):
                        readExonsB.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                    readExonsB.append( [offsets[j]+juncOffset, end+juncOffset] )

                self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(readExonsA[0][0]), readExonsA,  \
                                                         self.aligned.getChromosome(readExonsB[0][0]), readExonsB, junc.xs, junc.NH))


    def expandByCluster(self, f, out_name):
        clusters = binaryIO.readClusters(f)
        unspliced_index = binaryIO.readListFromFile(f)
        spliced_index = binaryIO.readListFromFile(f)

        for i in range(len(clusters)):
            self.aligned.exons = clusters[i]
            #print('Expanding (%d,%d)' % (self.aligned.exons[0], self.aligned.exons[-1]))

            self.expandUnsplicedBinary(f, unspliced_index[i])
            self.expandSplicedBinary(f, spliced_index[i])

            if i == 0:
                with open(out_name, 'w') as f2:
                    self.aligned.writeSAM(f2, header=True)
            else:
                with open(out_name, 'a') as f2:
                    self.aligned.writeSAM(f2, header=False)

            self.aligned.unspliced = []
            self.aligned.spliced = []
            self.aligned.paired = []


    def readIndexBinary(self, f):
        '''
            Read index information from beginning of file
        '''

        # Read an unzip index
        sLen = binaryIO.readVal(f, 4)
        s = self.expandString(f.read(sLen))
        startPos = 0

        self.exonChunkSize, startPos = binaryIO.binaryToVal(s, 4, startPos)
        self.junctionChunkSize, startPos = binaryIO.binaryToVal(s, 4, startPos)

        # Read junction names
        t = time.time()
        self.sortedJuncs, self.exonBytes, startPos = binaryIO.readJunctionsList(s, startPos)

        # Read junction chunk lengths
        self.junctionChunkLens, startPos = binaryIO.readList(s, startPos)

        self.unsplicedIndex = dict()
        self.unsplicedExonsIndex = dict()
        numNH, startPos = binaryIO.binaryToVal(s, 4, startPos)
        for _ in range(numNH):
            NH, startPos = binaryIO.binaryToVal(s, 4, startPos)

            # Read unspliced coverage index
            self.unsplicedIndex[NH], startPos = binaryIO.readList(s, startPos)

            # Read unspliced exons index
            self.unsplicedExonsIndex[NH], startPos = binaryIO.readList(s, startPos)

    def RLE(self, vector):
        rle = []

        val = vector[0]
        length = 0
        for v in vector:
            if v == val:
                length += 1
            else:
                if length == 1:
                    rle.append([val])
                else:
                    rle.append([val, length])
                val = v
                length = 1

        if length == 1:
            rle.append([val])
        else:
            rle.append([val, length])

        return rle

    def RLEtoVector(self, rle):
        ''' Convert a run length encoded vector to the original coverage vector '''

        vector = []
        for row in rle:
            vector += [row[0]] * row[1]
        return vector

    def expandString(self, s):
        ''' Use a predefined python library to expand the given string.
            Return the decompressed string '''

        if self.compressMethod == 0:
            return self.zlib.decompress(s)
        elif self.compressMethod == 1:
            return self.lzma.decompress(s)
        elif self.compressMethod == 2:
            return self.bz2.decompress(s)



    def getCoverage(self, compressedFilename, chrom, start=None, end=None):
        with open(compressedFilename, 'rb') as f:
            coverage = [0.0] * (end-start)
            chromosomes = binaryIO.readChroms(f)
            self.aligned = alignments.Alignments(chromosomes)
            if start == None or end == None:
                start = 0
                end = chromosomes[chrom]

            print('(%d,%d)' % (start, end))

            for k in sorted(chromosomes.keys()):
                if not k == chrom:
                    start += chromosomes[k]
                    end += chromosomes[k]
                else:
                    break

            clusters = binaryIO.readClusters(f)
            start_i, end_i = self.getRelevantClusters(clusters, start, end)
            if start_i == end_i:
                return [0] * (end-start)

            unspliced_index = binaryIO.readListFromFile(f)
            spliced_index = binaryIO.readListFromFile(f)

            skip = sum(unspliced_index[:start_i]) + sum(spliced_index[:start_i])
            f.seek(skip, 1)

            unspliced_t = 0.0
            spliced_t = 0.0

            for i in range(start_i, end_i):
                self.aligned.exons = clusters[i]

                start_t = time.time()
                coverage = self.getUnsplicedCoverage(f, unspliced_index[i], coverage, start, end)
                mid_t = time.time()
                coverage = self.getSplicedCoverage(f, spliced_index[i], coverage, start, end)
                end_t = time.time()
                unspliced_t += (mid_t - start_t)
                spliced_t += (end_t - mid_t)

        print('Unspliced time: %0.3fs' % (unspliced_t))
        print('Spliced time:   %0.3fs' % (spliced_t))
        print('')

        return coverage

    def getSplicedCoverage(self, f, length, coverage, range_start, range_end):
        s = self.expandString(f.read(length))

        startPos = 0
        readLenBytes, startPos = binaryIO.binaryToVal(s, 1, startPos)

        sorted_junctions, exonBytes, startPos = binaryIO.readJunctionsList(s, startPos)

        for key in sorted_junctions:
            exons = key[:-2]

            length = 0
            boundaries = []
            for i in range(len(exons)):
                e = exons[i]
                subexon_length = self.aligned.exons[e+1] - self.aligned.exons[e]

                length += subexon_length
                if i < len(exons)-1:
                    if i == 0:
                        boundaries.append(subexon_length)
                    else:
                        boundaries.append(boundaries[-1] + subexon_length)

            # Read the rest of the junction information
            junc, startPos = binaryIO.readJunction(s, junction.Junction(exons, length, boundaries), readLenBytes, startPos)
            junc.xs = key[-2]
            junc.NH = key[-1]

            junc_coverage = self.RLEtoVector(junc.coverage)
            for i in range(len(junc_coverage)):
                junc_coverage[i] /= float(junc.NH)

            juncBounds = []
            for j in junc.exons:
                juncBounds.append([self.aligned.exons[j], self.aligned.exons[j+1]])

            # marks indices in junction coverage vector where exons begin
            mapping = [0]
            for j in juncBounds:
                mapping.append(mapping[-1] + j[1] - j[0])

            genome_pos = self.aligned.exons[exons[0]] - range_start
            length = range_end - range_start
            map_id = 1
            for i in range(len(junc_coverage)):
                if i == mapping[map_id]:
                    genome_pos += juncBounds[map_id][0] - juncBounds[map_id-1][1]
                    map_id += 1

                if genome_pos >= 0:
                    if genome_pos >= length:
                        break
                    else:
                        coverage[genome_pos] += junc_coverage[i]

                genome_pos += 1



            '''
            # Expand coverage vector to include introns
            cov = []
            for i in range(len(junc.exons)-1):
                cov += junc_coverage[mapping[i]:mapping[i+1]] + [0] * (juncBounds[i+1][0]-juncBounds[i][1])
            cov += junc_coverage[mapping[-2]:mapping[-1]]

            cov_i = 0
            id_start = self.aligned.exons[exons[0]] - range_start
            if (id_start < 0):
                cov_i = -id_start
                id_start = 0

            id_end = min(self.aligned.exons[exons[-1]+1], range_end) - range_start

            for i in range(id_start, id_end):
                coverage[i] += cov[cov_i]
                cov_i += 1
            '''

        return coverage



    def getUnsplicedCoverage(self, f, length, coverage, range_start, range_end):
        s = self.expandString(f.read(length))

        start = 0
        num_NH, start = binaryIO.binaryToVal(s, 2, start)
        fragLenBytes, start = binaryIO.binaryToVal(s, 1, start)
        readLenBytes, start = binaryIO.binaryToVal(s, 1, start)

        for i in range(num_NH):
            NH, start = binaryIO.binaryToVal(s, 2, start)

            cov, start = binaryIO.readCov(s, start)
            cov = self.RLEtoVector(cov)

            cov_i = 0
            id_start = self.aligned.exons[0] - range_start
            if (id_start < 0):
                cov_i = -id_start
                id_start = 0

            id_end = min(self.aligned.exons[-1], range_end) - range_start

            for j in range(id_start, id_end):
                coverage[j] += cov[cov_i] / float(NH)
                cov_i += 1

            # Skip read length distributions
            for _ in range(len(self.aligned.exons)-1):
                _, start = binaryIO.skipLens(s, readLenBytes, start)
                paired, start = binaryIO.skipLens(s, fragLenBytes, start)
                if paired:
                    left, start = binaryIO.skipLens(s, readLenBytes, start)
                    if left:
                        _, start = binaryIO.skipLens(s, readLenBytes, start)

        return coverage

    def getRelevantClusters(self, c, start, end):
        '''

        :param clusters:
        :param start:
        :param end:
        :return: The start and end+1 index in clusters of all clusters that overlap [start,end)
        '''

        clusters = []
        for i in range(len(c)):
            clusters.append((c[i][0], c[i][-1]))

        for i in range(len(clusters)):
            c = clusters[i]
            if not len(c) == 2:
                print(i)
                print(c)
                exit()

        i = bisect.bisect_left(clusters, (start,start))

        if clusters[i-1][0] <= start and clusters[i-1][1] > start:
            start_i = i-1
        elif clusters[i][1] > start:
            start_i = i
        else:
            print('Error bisecting')
            print('(%d,%d)' % (start,end))
            print(clusters[i-3:i+3])
            exit()

        end_i = bisect.bisect_left(clusters, (end,end))

        return start_i, end_i


    def getReads(self, compressedFile, chrom, start=None, end=None):
        '''
        Return all reads for which at least one exon overlaps the given region
        :param compressedFile:
        :param chrom:
        :param start:
        :param end:
        :return:
        '''

        with open(compressedFile, 'rb') as f:
            chromosomes = binaryIO.readChroms(f)
            self.aligned = alignments.Alignments(chromosomes)
            if start == None or end == None:
                start = 0
                end = chromosomes[chrom]

            for k in sorted(chromosomes.keys()):
                if not k == chrom:
                    start += chromosomes[k]
                    end += chromosomes[k]
                else:
                    break

            clusters = binaryIO.readClusters(f)
            start_i, end_i = self.getRelevantClusters(clusters, start, end)
            if start_i == end_i:
                return [], []

            unspliced_index = binaryIO.readListFromFile(f)
            spliced_index = binaryIO.readListFromFile(f)

            skip = sum(unspliced_index[:start_i]) + sum(spliced_index[:start_i])
            f.seek(skip, 1)

            for i in range(start_i, end_i):
                self.aligned.exons = clusters[i]

                self.getUnsplicedReads(f, unspliced_index[i], start, end)
                self.getSplicedReads(f, spliced_index[i], start, end)

        return self.aligned.unspliced, self.aligned.paired

    def getSplicedReads(self, f, length, range_start, range_end):
        ''' Expand a file containing compressed spliced alignments
        '''

        s = self.expandString(f.read(length))

        startPos = 0
        readLenBytes, startPos = binaryIO.binaryToVal(s, 1, startPos)

        sorted_junctions, exonBytes, startPos = binaryIO.readJunctionsList(s, startPos)

        for key in sorted_junctions:
            exons = key[:-2]

            length = 0
            boundaries = []
            for i in range(len(exons)):
                e = exons[i]
                subexon_length = self.aligned.exons[e+1] - self.aligned.exons[e]

                length += subexon_length
                if i < len(exons)-1:
                    if i == 0:
                        boundaries.append(subexon_length)
                    else:
                        boundaries.append(boundaries[-1] + subexon_length)

            debug = False

            # Read the rest of the junction information
            junc, startPos = binaryIO.readJunction(s, junction.Junction(exons, length, boundaries), readLenBytes, startPos)
            junc.xs = key[-2]
            junc.NH = key[-1]

            junc.coverage = self.RLEtoVector(junc.coverage)

            unpaired, paired = self.aligned.findReads(junc.unpairedLens, junc.pairedLens, junc.lensLeft, junc.lensRight, junc.coverage, junc.boundaries, debug)

            juncBounds = []
            for j in junc.exons:
                juncBounds.append([self.aligned.exons[j], self.aligned.exons[j+1]])

            # Offset of start of junction from beginning of chromosome
            juncOffset = juncBounds[0][0]

            start_offset  = range_start - juncOffset
            end_offset = range_end - juncOffset

            # marks indices in junction coverage vector where exons are
            mapping = [0]
            # offsets of start of each exon in junction
            offsets = [0]
            for j in range(1,len(juncBounds)):
                mapping.append(mapping[-1] + juncBounds[j-1][1] - juncBounds[j-1][0])
                offsets.append(juncBounds[j][0] - juncBounds[0][0])

            for r in unpaired:
                start = r[0]
                i = 0
                while i < (len(mapping)-1) and mapping[i+1] <= start:
                    i += 1
                start += offsets[i] - mapping[i]

                end = r[1]
                j = 0
                while j < (len(mapping)-1) and mapping[j+1] < end:
                    j += 1
                end += offsets[j] - mapping[j]

                readExons = []
                if i == j:
                    readExons.append( [start+juncOffset, end+juncOffset] )
                else:
                    readExons.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                    for x in range(i+1,j):
                        readExons.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                    readExons.append( [offsets[j]+juncOffset, end+juncOffset] )

                if readExons[-1][1] > range_start and readExons[0][0] < range_end:
                    self.aligned.unspliced.append(read.Read(self.aligned.getChromosome(readExons[0][0]), readExons, junc.xs, junc.NH))

            for p in paired:
                start = p[0][0]
                i = 0
                while i < (len(mapping)-1) and mapping[i+1] <= start:
                    i += 1
                start += offsets[i] - mapping[i]

                end = p[0][1]
                j = 0
                while j < (len(mapping)-1) and mapping[j+1] < end:
                    j += 1
                end += offsets[j] - mapping[j]

                readExonsA = []
                if i == j:
                    readExonsA.append( [start+juncOffset, end+juncOffset] )
                else:
                    readExonsA.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                    for x in range(i+1,j):
                        readExonsA.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                    readExonsA.append( [offsets[j]+juncOffset, end+juncOffset] )

                start = p[1][0]
                i = 0
                while i < (len(mapping)-1) and mapping[i+1] <= start:
                    i += 1
                start += offsets[i] - mapping[i]

                end = p[1][1]
                j = 0
                while j < (len(mapping)-1) and mapping[j+1] < end:
                    j += 1
                end += offsets[j] - mapping[j]

                readExonsB = []
                if i == j:
                    readExonsB.append( [start+juncOffset, end+juncOffset] )
                else:
                    readExonsB.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                    for x in range(i+1,j):
                        readExonsB.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                    readExonsB.append( [offsets[j]+juncOffset, end+juncOffset] )

                if readExonsB[-1][1] > range_start and readExonsA[0][0] < range_end:
                    self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(readExonsA[0][0]), readExonsA,  \
                                                                     self.aligned.getChromosome(readExonsB[0][0]), readExonsB, junc.xs, junc.NH))


    def getUnsplicedReads(self, f, length, range_start, range_end):
        '''

        :param f:
        :param length:
        :param start:
        :param end:
        :return: Set of reads extracted from this chunk. Discard any reads completely outside [start,end)
        '''
        s = self.expandString(f.read(length))

        start = 0
        num_NH, start = binaryIO.binaryToVal(s, 2, start)
        fragLenBytes, start = binaryIO.binaryToVal(s, 1, start)
        readLenBytes, start = binaryIO.binaryToVal(s, 1, start)

        for i in range(num_NH):
            NH, start = binaryIO.binaryToVal(s, 2, start)

            cov, start = binaryIO.readCov(s, start)
            cov = self.RLEtoVector(cov)

            for i in range(len(self.aligned.exons)-1):

                unpairedLens, start = binaryIO.readLens(s, readLenBytes, start)
                pairedLens, start = binaryIO.readLens(s, fragLenBytes, start)
                if len(pairedLens) > 0:
                    lensLeft, start = binaryIO.readLens(s, readLenBytes, start)
                    if len(lensLeft) > 0:
                        lensRight, start = binaryIO.readLens(s, readLenBytes, start)
                    else:
                        lensRight = dict()
                else:
                    lensLeft = dict()
                    lensRight = dict()

                if len(unpairedLens) > 0 or len(pairedLens) > 0:
                    exonStart = self.aligned.exons[i]
                    exonEnd = self.aligned.exons[i+1]

                    start_offset  = range_start - exonStart
                    end_offset = range_end - exonStart

                    unpaired, paired = self.aligned.findReads(unpairedLens, pairedLens, lensLeft, lensRight, cov[exonStart-self.aligned.exons[0]:exonEnd-self.aligned.exons[0]], None, debug=False)

                    for r in unpaired:
                        if r[1] > start_offset and r[0] < end_offset:
                            self.aligned.unspliced.append(read.Read(self.aligned.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))
                    for p in paired:
                        if p[1][1] > start_offset and p[0][0] < end_offset:
                            self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                                             self.aligned.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH))

    def readOverlapsRegion(self, exons, start, end):
        '''

        :param exons: List of exon bounds
        :param start: Start of region of interest
        :param end: End of region of interest
        :return: True if any of the exons overlap the region of interest
        '''

        for e in exons:
            if e[0] >= end:
                # Passed region
                return False
            elif e[1] > start:
                return True
        return False


    '''
    TODO
    These functions are out of date and need to be revised/removed
    '''

    def expandSpliced(self, f):
        ''' Expand a file containing compressed spliced alignments
        '''

        #self.aligned.exons = None
        junc = None

        self.aligned.spliced = []
        for line in f:
            # Lines grouped by junction
            if line[0] == '>' or line[0] == '/':
                if not junc == None:
                    # process junction
                    unpaired, paired = self.aligned.findReads(junc.unpairedLens, junc.pairedLens, junc.lensLeft, junc.lensRight, junc.coverage)

                    juncBounds = []
                    for j in junctionExons:
                        juncBounds.append([self.aligned.exons[j], self.aligned.exons[j+1]])

                    # Offset of start of junction from beginning of chromosome
                    juncOffset = juncBounds[0][0]

                    # marks indices in junction coverage vector where exons are
                    mapping = [0]
                    # offsets of start of each exon in junction
                    offsets = [0]
                    for j in range(1,len(juncBounds)):
                        mapping.append(mapping[-1] + juncBounds[j-1][1] - juncBounds[j-1][0])
                        offsets.append(juncBounds[j][0] - juncBounds[0][0])

                    for r in unpaired:
                        start = r[0]
                        i = 0
                        while i < (len(mapping)-1) and mapping[i+1] <= start:
                            i += 1
                        start += offsets[i] - mapping[i]

                        end = r[1]
                        j = 0
                        while j < (len(mapping)-1) and mapping[j+1] < end:
                            j += 1
                        end += offsets[j] - mapping[j]

                        readExons = []
                        if i == j:
                            readExons.append( [start+juncOffset, end+juncOffset] )
                        else:
                            readExons.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                            for x in range(i+1,j):
                                readExons.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                            readExons.append( [offsets[j]+juncOffset, end+juncOffset] )

                        self.aligned.spliced.append(read.Read(self.aligned.getChromosome(readExons[0][0]), readExons, junc.xs, junc.NH))

                    for p in paired:
                        start = p[0][0]
                        i = 0
                        while i < (len(mapping)-1) and mapping[i+1] <= start:
                            i += 1
                        start += offsets[i] - mapping[i]

                        end = p[0][1]
                        j = 0
                        while j < (len(mapping)-1) and mapping[j+1] < end:
                            j += 1
                        end += offsets[j] - mapping[j]

                        readExonsA = []
                        if i == j:
                            readExonsA.append( [start+juncOffset, end+juncOffset] )
                        else:
                            readExonsA.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                            for x in range(i+1,j):
                                readExonsA.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                            readExonsA.append( [offsets[j]+juncOffset, end+juncOffset] )

                        start = p[1][0]
                        i = 0
                        while i < (len(mapping)-1) and mapping[i+1] <= start:
                            i += 1
                        start += offsets[i] - mapping[i]

                        end = p[1][1]
                        j = 0
                        while j < (len(mapping)-1) and mapping[j+1] < end:
                            j += 1
                        end += offsets[j] - mapping[j]

                        readExonsB = []
                        if i == j:
                            readExonsB.append( [start+juncOffset, end+juncOffset] )
                        else:
                            readExonsB.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                            for x in range(i+1,j):
                                readExonsB.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                            readExonsB.append( [offsets[j]+juncOffset, end+juncOffset] )

                        self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(readExonsA[0][0]), readExonsA,  \
                                                                 self.aligned.getChromosome(readExonsB[0][0]), readExonsB, junc.xs, junc.NH))

                if line[0] == '/':
                    return

                # Start of a new junction
                key = line[1:].rstrip().split('\t')
                junctionExons = [int(e) for e in key[:-2]]

                length = 0
                boundaries = []
                for i in range(len(junctionExons)):
                    e = junctionExons[i]
                    subexon_length = self.aligned.exons[e+1] - self.aligned.exons[e]

                    length += subexon_length
                    if i == 0:
                        boundaries.append(subexon_length)
                    elif i < len(junctionExons)-1:
                        boundaries.append(boundaries[-1] + subexon_length)

                junc = junction.Junction(junctionExons, length, boundaries)
                junc.unpairedLens = None
                junc.pairedLens = None
                junc.lensLeft = None
                junc.lensRight = None

                junc.xs = key[-2]
                junc.NH = int(key[-1])

            else:
                if junc.unpairedLens == None:
                    # First line after '>' contains unpaired read length distribution
                    junc.unpairedLens = dict()
                    for readLen in line.rstrip().split('\t'):
                        readLen = readLen.split(',')
                        junc.unpairedLens[int(readLen[0])] = int(readLen[1])
                    junc.coverage = []
                if junc.pairedLens == None:
                    # Second line after '>' contains paired read length distribution
                    junc.pairedLens = dict()
                    for readLen in line.rstrip().split('\t'):
                        readLen = readLen.split(',')
                        junc.pairedLens[int(readLen[0])] = int(readLen[1])
                    junc.coverage = []
                elif junc.lensLeft == None:
                    # Third line after '>' contains left read length distribution
                    junc.lensLeft = dict()
                    if len(line.rstrip()) > 0:
                        for length in line.rstrip().split('\t'):
                            length = length.split(',')
                            junc.lensLeft[int(length[0])] = int(length[1])
                elif junc.lensRight == None:
                    # Fourth line after '>' contains right read length distribution
                    junc.lensRight = dict()
                    if len(line.rstrip()) > 0:
                        for length in line.rstrip().split('\t'):
                            length = length.split(',')
                            junc.lensRight[int(length[0])] = int(length[1])
                else:
                    # Rest of lines contain run-length-encoded coverage vector
                    row = line.rstrip().split('\t')
                    if len(row) == 1:
                        length = 1
                    else:
                        length = int(row[1])
                    val = int(row[0])

                    junc.coverage += [val] * length
                    #junc.coverage += [[val, length]]

    def expandUnspliced(self, f):
        ''' Expand a file containing compressed unspliced alignments
        '''

        NH = 1
        unpairedLens = dict()
        pairedLens = dict()
        exonStart = 0
        exonEnd = 0
        for line in f:
            if line[0] == '#':
                if len(unpairedLens) > 0 or len(pairedLens) > 0:
                    unpaired, paired = self.aligned.findReads(unpairedLens, pairedLens, lensLeft, lensRight, coverage[exonStart:exonEnd])

                    for r in unpaired:
                        self.aligned.unspliced.append(read.Read(self.aligned.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))
                    for p in paired:
                        self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                                 self.aligned.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH))

                NH = int(line.rstrip()[1:])

                lineNum = 0
                unpairedLens = dict()
                pairedLens = dict()
                lensLeft = dict()
                lensRight = dict()
                coverage = []
                RLE_segment = True

            elif RLE_segment and not line[0] == '>':
                # Rest of lines contain run-length-encoded coverage vector
                row = line.rstrip().split("\t")

                if len(row) == 1:
                    length = 1
                else:
                    length = int(row[1])

                coverage += [int(row[0])] * length
                #coverage += [[int(row[0]), length]]

            elif line[0] == '>':
                RLE_segment = False
                if len(unpairedLens) > 0 or len(pairedLens) > 0:
                    # Process previous exon
                    unpaired, paired = self.aligned.findReads(unpairedLens, pairedLens, lensLeft, lensRight, coverage[exonStart:exonEnd])

                    for r in unpaired:
                        self.aligned.unspliced.append(read.Read(self.aligned.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))

                    for p in paired:
                        self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                                 self.aligned.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH ))
                #exonId += 1
                lineNum = 0
                exonId = int(line[1:].rstrip())
                bounds = line[1:].rstrip().split('-')
                exonStart = self.aligned.exons[exonId]
                exonEnd = self.aligned.exons[exonId+1]

            elif lineNum == 1:
                # First line after '>' contains unpaired lengths
                unpairedLens = dict()
                for row in line.rstrip().split('\t'):
                    if len(row) > 0:
                        readLen = row.split(',')
                        unpairedLens[int(readLen[0])] = int(readLen[1])
            elif lineNum == 2:
                # Second line after '>' contains paired lengths
                pairedLens = dict()
                for row in line.rstrip().split('\t'):
                    if len(row) > 0:
                        readLen = row.split(',')
                        pairedLens[int(readLen[0])] = int(readLen[1])
            elif lineNum == 3:
                # Second line after '>' contains left read length distribution
                lensLeft = dict()
                if len(line.rstrip()) > 0:
                    for length in line.rstrip().split('\t'):
                        length = length.split(',')
                        lensLeft[int(length[0])] = int(length[1])
            elif lineNum == 4:
                # Third line after '>' contains right read length distribution
                lensRight = dict()
                if len(line.rstrip()) > 0:
                    for length in line.rstrip().split('\t'):
                        length = length.split(',')
                        lensRight[int(length[0])] = int(length[1])
            lineNum += 1

        # Process the final exon
        if len(unpairedLens) > 0 or len(pairedLens) > 0:
            unpaired, paired = self.aligned.findReads(unpairedLens, pairedLens, lensLeft, lensRight, coverage[exonStart:exonEnd])

            for r in unpaired:
                self.aligned.unspliced.append(read.Read(self.aligned.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))

            for p in paired:
                self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                         self.aligned.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH))

    '''
    def getCoverage2(self, compressedFilename, chrom, start=None, end=None):
        self.aligned = None
        with open(compressedFilename, 'rb') as f:
            self.readIndexBinary(f)

            chromosomes = binaryIO.readChroms(f)
            self.aligned = alignments.Alignments(chromosomes)
            self.aligned.exons = binaryIO.readExons(f, self.exonBytes)

            if start == None or end == None:
                start = 0
                end = chromosomes[chrom]

            for k in sorted(chromosomes.keys()):
                if not k == chrom:
                    start += chromosomes[k]
                    end += chromosomes[k]
                else:
                    break

            coverage = [0.0] * (end-start)
            startTime = time.time()
            coverage = self.getUnsplicedCoverage2(f, coverage, start, end)
            #print(coverage[123823:123833])
            t = time.time() - startTime
            #print('  Unspliced time: %0.3f s' % t)
            startTime = time.time()
            coverage = self.getSplicedCoverage2(f, coverage, start, end)
            #print(coverage[123823:123833])
            t = time.time() - startTime
            #print('  Spliced time: %0.3f s' % t)

        return coverage

    def getSplicedCoverage2(self, filehandle, coverage, start, end):

        #print('Getting spliced coverage 2')

        readLenBytes = binaryIO.readVal(filehandle, 1)

        relevantExonsStart = 0
        while self.aligned.exons[relevantExonsStart+1] < start:
            relevantExonsStart += 1
        relevantExonsEnd = relevantExonsStart
        while self.aligned.exons[relevantExonsEnd] < end:
            relevantExonsEnd += 1

        relevant = False

        # Keep track of which chunk we are in, and what position in that chunk
        chunkId = 0
        chunkOffset = 0
        idInChunk = 0

        # Since the junctions are sorted, if all junctions are after the region of interest, we are done
        passedRange = False

        s = None
        lastJ = 0

        expandTime = 0
        skipTime = 0
        readTime = 0
        updateTime = 0
        t = time.time()

        for j in range(len(self.sortedJuncs)):
            junc = self.sortedJuncs[j]
            exons = [int(e) for e in junc[:-2]]
            NH = float(junc[-1])

            
            passed = True
            relevant = False
            for e in exons:
                if e < relevantExonsEnd:
                    passed = False
                    if e >= relevantExonsStart:
                        relevant = True
                        break
            if passed:
                break

            if relevant:
                length = 0
                for e in exons:
                    length += self.aligned.exons[e+1] - self.aligned.exons[e]

                t1 = time.time()
                if s == None:
                    if chunkOffset > 0:
                        filehandle.seek(chunkOffset, 1)
                        chunkOffset = 0
                    s = self.expandString(filehandle.read(self.junctionChunkLens[chunkId]))
                    startPos = 0

                t2 = time.time()
                expandTime += t2-t1

                for i in range(lastJ, j):
                    startPos = binaryIO.skipJunction(s, readLenBytes, startPos)

                t3 = time.time()
                skipTime += t3-t2

                junc, startPos = binaryIO.readJunction(s, junction.Junction(exons, length), readLenBytes, startPos)
                junc.coverage = self.RLEtoVector(junc.coverage)

                t4 = time.time()
                readTime += t4-t3

                
                # Update coverage vector
                currExon = 0
                exonLength = self.aligned.exons[exons[currExon]+1] - self.aligned.exons[exons[currExon]]
                offsetInExon = 0
                for i in range(length):
                    if offsetInExon == exonLength:
                        currExon += 1
                        exonLength = self.aligned.exons[exons[currExon]+1] - self.aligned.exons[exons[currExon]]
                        offsetInExon = 0

                    if junc.coverage[i] > 0:
                        pos = self.aligned.exons[exons[currExon]] + offsetInExon
                        if pos >= start and pos < end:
                            coverage[pos - start] += junc.coverage[i] / NH

                    offsetInExon += 1

                currExon = 0
                exonLength = self.aligned.exons[exons[currExon]+1] - self.aligned.exons[exons[currExon]]
                offsetInExon = 0

                for row in junc.coverage:
                    c = row[0]
                    l = row[1]
                    while l+offsetInExon >= exonLength:
                        if c > 0:
                            posEnd = min(self.aligned.exons[exons[currExon]+1]-start, end-start)
                            if posEnd > 0:
                                posStart = max(self.aligned.exons[exons[currExon]] + offsetInExon - start, 0)
                                for i in range(posStart, posEnd):
                                    coverage[i] += c / NH
                        l = l + offsetInExon - exonLength
                        currExon += 1

                        if currExon < len(exons):
                            if self.aligned.exons[exons[currExon]] >= end:
                                break
                            exonLength = self.aligned.exons[exons[currExon]+1] - self.aligned.exons[exons[currExon]]
                        offsetInExon = 0

                    if currExon >= len(exons) or self.aligned.exons[exons[currExon]] >= end:
                        break
                    
                    if l > 0:
                        if c > 0:
                            posEnd = min(self.aligned.exons[exons[currExon]] + offsetInExon + l - start, end-start)
                            if posEnd > 0:
                                posStart = max(self.aligned.exons[exons[currExon]] + offsetInExon - start, 0)
                                for i in range(posStart, posEnd):
                                    coverage[i] += c / NH

                                if posEnd == end-start:
                                    break
                        offsetInExon += l


                

                t5 = time.time()
                updateTime += t5-t4

                lastJ = j+1


            idInChunk += 1
            if idInChunk == self.junctionChunkSize:
                if s == None:
                    chunkOffset += self.junctionChunkLens[chunkId]
                else:
                    s = None
                chunkId += 1
                idInChunk = 0
                lastJ = j+1

        # Now skip to the end of the spliced section
        skip = 0
        while chunkId < len(self.junctionChunkLens):
            skip += self.junctionChunkLens[chunkId]
            chunkId += 1
        filehandle.seek(skip, 1)

        return coverage
    
    def getUnsplicedCoverage2(self, filehandle, coverage, start, end):
        #print('Getting unspliced coverage 2')
        NHs = binaryIO.readVal(filehandle, 2)
        self.sectionLen = binaryIO.readVal(filehandle, 4)
        fragLenBytes = binaryIO.readVal(filehandle, 1)
        readLenBytes = binaryIO.readVal(filehandle, 1)

        for _ in range(NHs):
            NH = binaryIO.readVal(filehandle, 2)
            offset = 0
            fileOffset = 0
            for i in range(len(self.unsplicedIndex[NH])):
                startOffset = offset
                if offset+self.sectionLen > start:
                    if fileOffset > 0:
                        filehandle.seek(fileOffset, 1)
                        fileOffset = 0

                    covLength = self.unsplicedIndex[NH][i]
                    if covLength > 0:
                        s = filehandle.read(covLength)
                        cov,_ = binaryIO.readCov(self.expandString(s))

                        for row in cov:
                            c = row[0]
                            l = row[1]
                            if c > 0 and offset+l > start:
                                if offset < start:
                                    startPos = 0
                                else:
                                    startPos = offset - start

                                if offset + l >= end:
                                    endPos = end - start
                                else:
                                    endPos = offset + l - start

                                for j in range(startPos, endPos):
                                    coverage[j] += c / float(NH)
                            offset += l
                            if offset >= end:
                                break
                    else:
                        offset += self.sectionLen
                else:
                    fileOffset += self.unsplicedIndex[NH][i]
                    offset += self.sectionLen
                if offset >= end:
                    break

            # Skip to end of coverage
            skip = 0
            for j in range(i+1, len(self.unsplicedIndex[NH])):
                skip += self.unsplicedIndex[NH][j]

            # Skip exon information
            for e in range(len(self.unsplicedExonsIndex[NH])):
                skip += self.unsplicedExonsIndex[NH][e]

            filehandle.seek(skip, 1)
        
        return coverage

    def getReads(self, compressedFile, chrom, start=None, end=None):
        with open(compressedFile, 'rb') as f:
            self.readIndexBinary(f)

            chromosomes = binaryIO.readChroms(f)
            self.aligned = alignments.Alignments(chromosomes)
            self.aligned.exons = binaryIO.readExons(f, self.exonBytes)

            if start == None or end == None:
                start = 0
                end = chromosomes[chrom]

            for k in sorted(chromosomes.keys()):
                if not k == chrom:
                    start += chromosomes[k]
                    end += chromosomes[k]
                else:
                    break

            reads = self.getSplicedReads(f, start, end)
            reads += self.getUnsplicedReads(f, start, end, self.exonBytes)

        return reads

    def getSplicedReads(self, filehandle, start, end):
        reads = []

        relevantExonsStart = 0
        while self.aligned.exons[relevantExonsStart+1] < start:
            relevantExonsStart += 1
        relevantExonsEnd = relevantExonsStart
        while self.aligned.exons[relevantExonsEnd] < end:
            relevantExonsEnd += 1


        junc = None
        relevant = False

        # Keep track of which chunk we are in, and what position in that chunk
        chunkId = 0
        chunkOffset = 0
        idInChunk = 0

        # For each junction in the chunk, 1 if it is relevant, 0 if not
        relevant = [0] * self.junctionChunkSize

        # Id in chunk of the last relevant junction
        # This way after we read the last relevant junction we can discard the rest without processing it
        lastRelevant = -1

        # Since the junctions are sorted, if all junctions are after the region of interest, we are done
        passedRange = False

        junctionExons = []
        junctionLengths = []

        for j in self.sortedJuncs:
            exons = [int(e) for e in j[:-2]]

            lastRelevant = -1
            passed = True
            length = 0
            for e in exons:
                length += self.aligned.exons[e+1] - self.aligned.exons[e]
                if e < relevantExonsEnd:
                    passed = False
                    if e >= relevantExonsStart:
                        relevant[idInChunk] = 1
                        lastRelevant = idInChunk
            if passed:
                passedRange = True

            junctionExons += [exons]
            junctionLengths += [length]

            NH = 1.0 / float(j[-1])

            idInChunk += 1
            if idInChunk == self.junctionChunkSize:
                # If any junctions in this chunk are relevant, read the chunk and process the relevant junctions
                if lastRelevant >= 0:
                    filehandle.seek(chunkOffset)
                    s = self.expandString(filehandle.read(self.junctionChunkLens[chunkId]))

                    for i in range(lastRelevant+1):
                        junc, s = binaryIO.readJunction(s, junction.Junction(junctionExons[i], juntionLengths[i]), readLenBytes)

                        if relevant[i] == 1:
                            # Process the junction
                            unpaired, paired = self.aligned.findReads(junc.unpairedLens, junc.pairedLens, junc.lensLeft, junc.lensRight, junc.coverage)

                            juncBounds = []
                            for j in junc.exons:
                                juncBounds.append([self.aligned.exons[j], self.aligned.exons[j+1]])

                            # Offset of start of junction from beginning of chromosome
                            juncOffset = juncBounds[0][0]

                            # marks indices in junction coverage vector where exons are
                            mapping = [0]
                            # offsets of start of each exon in junction
                            offsets = [0]
                            for j in range(1,len(juncBounds)):
                                mapping.append(mapping[-1] + juncBounds[j-1][1] - juncBounds[j-1][0])
                                offsets.append(juncBounds[j][0] - juncBounds[0][0])

                            for r in unpaired:
                                start = r[0]
                                i = 0
                                while i < (len(mapping)-1) and mapping[i+1] <= start:
                                    i += 1
                                start += offsets[i] - mapping[i]

                                end = r[1]
                                j = 0
                                while j < (len(mapping)-1) and mapping[j+1] < end:
                                    j += 1
                                end += offsets[j] - mapping[j]
                                
                                readExons = []
                                if i == j:
                                    readExons.append( [start+juncOffset, end+juncOffset] )
                                else:
                                    readExons.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                                    for x in range(i+1,j):
                                        readExons.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                                    readExons.append( [offsets[j]+juncOffset, end+juncOffset] ) 

                                if readExons[0][0] >= start and readExons[-1][1] <= end:
                                    reads.append(read.Read(self.aligned.getChromosome(readExons[0][0]), readExons, junc.xs, junc.NH))

                            for p in paired:
                                start = p[0][0]
                                i = 0
                                while i < (len(mapping)-1) and mapping[i+1] <= start:
                                    i += 1
                                start += offsets[i] - mapping[i]

                                end = p[0][1]
                                j = 0
                                while j < (len(mapping)-1) and mapping[j+1] < end:
                                    j += 1
                                end += offsets[j] - mapping[j]
                                
                                readExonsA = []
                                if i == j:
                                    readExonsA.append( [start+juncOffset, end+juncOffset] )
                                else:
                                    readExonsA.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                                    for x in range(i+1,j):
                                        readExonsA.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                                    readExonsA.append( [offsets[j]+juncOffset, end+juncOffset] )

                                start = p[1][0]
                                i = 0
                                while i < (len(mapping)-1) and mapping[i+1] <= start:
                                    i += 1
                                start += offsets[i] - mapping[i]

                                end = p[1][1]
                                j = 0
                                while j < (len(mapping)-1) and mapping[j+1] < end:
                                    j += 1
                                end += offsets[j] - mapping[j]
                                
                                readExonsB = []
                                if i == j:
                                    readExonsB.append( [start+juncOffset, end+juncOffset] )
                                else:
                                    readExonsB.append( [start+juncOffset, offsets[i]+mapping[i+1]-mapping[i]+juncOffset] )

                                    for x in range(i+1,j):
                                        readExonsB.append( [offsets[x]+juncOffset, offsets[x]+mapping[x+1]-mapping[x]+juncOffset] )

                                    readExonsB.append( [offsets[j]+juncOffset, end+juncOffset] )

                                if readExonsA[0][0] >= start and readExonsB[-1][1] <= end:
                                    reads.append(pairedread.PairedRead(self.aligned.getChromosome(readExonsA[0][0]), readExonsA,  \
                                                                       self.aligned.getChromosome(readExonsB[0][0]), readExonsB, junc.xs, junc.NH))

                relevant = [0] * self.junctionChunkSize
                lastRelevant = -1

                idInChunk = 0
                chunkOffset += self.junctionChunkLens[chunkId]
                chunkId += 1

                junctionExons = []
                junctionLengths = []

                if passedRange:
                    break

        # Now skip to the end of the spliced section
        skip = 0
        while chunkId < len(self.junctionChunkLens):
            skip += self.junctionChunkLens[chunkId]
            chunkId += 1
        filehandle.seek(skip, 1)

        return reads
    
    def getUnsplicedReads(self, filehandle, start, end):
        NHs = binaryIO.readVal(f, 2)
        self.sectionLen = binaryIO.readVal(f, 4)
        fragLenBytes = binaryIO.readVal(f, 1)
        readLenBytes = binaryIO.readVal(f, 1)

        reads = []

        for _ in range(NHs):
            NH = binaryIO.readVal(f, 2)
            offset = 0
            fileOffset = 0
            coverage = []
            for s in range(len(self.unsplicedIndex[NH])):
                offset += self.sectionLen
                if offset > start:
                    if fileOffset > 0:
                        filehandle.seek(fileOffset, 1)
                        fileOffset = 0
                        covStart = offset - self.sectionLen

                    covLength = self.unsplicedIndex[NH][s]
                    if covLength > 0:
                        s = filehandle.read(covLength)
                        cov,_ = binaryIO.readCov(self.expandString(s))

                    coverage += cov
                else:
                    fileOffset += self.unsplicedIndex[NH][s]

                if offset >= end:
                    covEnd = offset
                    break

            # Skip to end of coverage
            skip = 0
            for i in range(s+1, len(self.unsplicedIndex[NH])):
                skip += self.unsplicedIndex[NH]

            # Read exon information and find reads
            exonId = 0
            for e in range(len(self.unsplicedExonsIndex[NH])):
                if self.aligned.exons[exonId + self.exonChunkSize] <= start:
                    skip += self.unsplicedExonsIndex[NH][e]
                else:
                    if skip > 0:
                        filehandle.seek(skip, 1)
                        skip = 0

                    lenDists = self.expandString(f.read(self.unsplicedExonsIndex[NH][e]))
                    i = 0
                    while len(lenDists) > 0:
                        i += 1
                        unpairedLens, lenDists = binaryIO.readLens(lenDists, readLenBytes)
                        pairedLens, lenDists = binaryIO.readLens(lenDists, fragLenBytes)
                        lensLeft = dict()
                        lensRight = dict()
                        if len(pairedLens) > 0:
                            lensLeft, lenDists = binaryIO.readLens(lenDists, readLenBytes)
                            if len(lensLeft) > 0:
                                lensRight, lenDists = binaryIO.readLens(lenDists, readLenBytes)

                        if len(pairedLens) > 0 or len(unpairedLens) > 0:
                            exonStart = self.aligned.exons[e*self.exonChunkSize + i - 1]
                            exonEnd = self.aligned.exons[e*self.exonChunkSize + i]
                            if exonStart < end and exonEnd > start:
                                unpaired, paired = self.aligned.findReads(unpairedLens, pairedLens, lensLeft, lensRight, coverage[exonStart-covStart:exonEnd-covStart])

                                for r in unpaired:
                                    if r[0]+exonStart >= start and r[1]+exonStart <= end:
                                        reads.append(read.Read(self.aligned.getChromosome(r[0]+exonStart), [[r[0]+exonStart, r[1]+exonStart]], None, NH))
                                for p in paired:
                                    if p[0][0]+exonStart >= start and p[1][1]+exonStart <= end:
                                        reads.append(pairedread.PairedRead(self.aligned.getChromosome(p[0][0]+exonStart), [[p[0][0]+exonStart, p[0][1]+exonStart]],   \
                                                                         self.aligned.getChromosome(p[1][0]+exonStart), [[p[1][0]+exonStart, p[1][1]+exonStart]], None, NH))



                if self.aligned.exons[exonId] >= end:
                    break

            # Skip the rest of the exons
            skip = 0
            for i in range(e+1, len(self.unsplicedExonsIndex[NH])):
                skip += self.unsplicedExonsIndex[NH][i]
            filehandle.seek(skip, 1)
        
        return reads
    '''

    def getGenes(self, compressedFile, chrom, start=None, end=None, overlapRadius=50):
        '''
            Return all genes within the given range in the form [(x_0,y_0), (x_1,y_1), ...] with a tuple indicating the bounds of each gene
        '''

        genes = []
        with open(compressedFile, 'r') as f:
            chromosomes = self.readHeader(f)
            if not chrom in chromosomes:
                print('Error! Chromosome name not recognized!')
                print('Chromosomes: ' + ', '.join(chromosomes.keys()))
                exit()

            self.aligned = alignments.Alignments(chromosomes)
            self.aligned.exons = self.readExons(f)

            splicedIndex, unsplicedIndex = self.readIndex(f)
            startPos = f.tell()

            if start == None or end == None:
                start = 0
                end = chromosomes[chrom]

            chromOffset = 0
            for k in sorted(chromosomes.keys()):
                if not k == chrom:
                    chromOffset += chromosomes[k]
                else:
                    break
            start += chromOffset
            end += chromOffset

            genes = self.getSplicedGeneBounds(f, start, end, overlapRadius, splicedIndex, startPos)
            #print('Spliced genes:')
            #print([(g[0]-chromOffset, g[1]-chromOffset) for g in sorted(genes)[:10]])
            genes2 = self.getUnsplicedGeneBounds(f, start, end, overlapRadius, unsplicedIndex, startPos)
            #print('Unspliced genes:')
            #print([(g[0]-chromOffset, g[1]-chromOffset) for g in sorted(genes2)[:10]])
            genes += genes2

            genes.sort()

            i = 0
            while i < (len(genes)-1):
                while i < (len(genes)-1) and genes[i+1][0] - genes[i][1] <= overlapRadius:
                    genes[i] = (genes[i][0], max(genes[i][1], genes[i+1][1]))
                    del genes[i+1]
                i += 1

            i = 0
            while i < (len(genes)-1):
                while i < (len(genes)-1) and genes[i+1][0] - genes[i][1] <= overlapRadius:
                    genes[i] = (genes[i][0], max(genes[i][1], genes[i+1][1]))
                    del genes[i+1]
                i += 1

        for i in range(len(genes)):
            genes[i] = (genes[i][0]-chromOffset, genes[i][1]-chromOffset)
        return genes

    def getSplicedGeneBounds(self, filehandle, start, end, overlapRadius, splicedIndex=None, startPos=0):
        if not splicedIndex == None:
            numExons = len(self.aligned.exons)
            i = 0
            while i < numExons and self.aligned.exons[i+1] <= start:
                i += 1
            j = i
            while j < numExons and self.aligned.exons[j] < end:
                j += 1

            minPos = splicedIndex[i]
            for n in range(i,j):
                if splicedIndex[n] < minPos:
                    minPos = splicedIndex[n]
            filehandle.seek(startPos + minPos)

        junc = None
        relevant = False

        relevantExonsStart = 0
        while self.aligned.exons[relevantExonsStart+1] < start:
            relevantExonsStart += 1
        relevantExonsEnd = relevantExonsStart
        while self.aligned.exons[relevantExonsEnd] < end:
            relevantExonsEnd += 1

        geneBounds = []
        geneStart = None
        geneEnd = None

        relevant = False

        for line in filehandle:
            # Lines grouped by junction
            if line[0] == '/':
                # Store bounds from last junction
                if not geneStart == None:
                    if geneEnd == None:
                        print('Error! Start but no end!')
                        exit()
                    geneBounds.append((geneStart, geneEnd))
                return geneBounds

            elif line[0] == '>':
                # Store bounds from last junction
                if not geneStart == None:
                    if geneEnd == None:
                        print('Error! Start but no end!')
                        exit()
                    if geneStart >= start and geneEnd <= end:
                        geneBounds.append((geneStart, geneEnd))
                geneStart = None
                geneEnd = None

                # Start of a new junction
                key = line[1:].rstrip().split('\t')
                junctionExons = [int(e) for e in key[:-2]]

                relevant = True
                passedRange = True
                for e in junctionExons:
                    if e <= relevantExonsEnd:
                        passedRange = False
                    if e < relevantExonsStart or e > relevantExonsEnd:
                        relevant = False
                        break
                if passedRange:
                    break

                count = 0
                currExon = 0
                currExonLen = self.aligned.exons[junctionExons[0]+1] - self.aligned.exons[junctionExons[0]]
                offsetInExon = 0

            elif relevant:
                if count < 3:
                    count += 1
                elif count == 3:
                    row = line.rstrip().split('\t')
                    if len(row) == 1:
                        length = 1
                    else:
                        length = int(row[1])
                    val = int(row[0])

                    oldPos = self.aligned.exons[junctionExons[currExon]]
                    offsetInExon += length
                    while offsetInExon > currExonLen:
                        currExon += 1
                        offsetInExon -= currExonLen
                        currExonLen = self.aligned.exons[junctionExons[currExon]+1] - self.aligned.exons[junctionExons[currExon]]

                    if val > 0:
                        if oldPos >= start:
                            geneStart = oldPos
                            geneEnd = self.aligned.exons[junctionExons[currExon]] + offsetInExon
                        else:
                            relevant = False
                    else:
                        if self.aligned.exons[junctionExons[currExon]] + offsetInExon >= start:
                            geneStart = self.aligned.exons[junctionExons[currExon]] + offsetInExon
                        else:
                            relevant = False

                    count += 1
                else:
                    # Rest of lines contain run-length-encoded coverage vector
                    row = line.rstrip().split('\t')
                    if len(row) == 1:
                        length = 1
                    else:
                        length = int(row[1])

                    offsetInExon += length
                    while offsetInExon > currExonLen:
                        currExon += 1
                        offsetInExon -= currExonLen
                        currExonLen = self.aligned.exons[junctionExons[currExon]+1] - self.aligned.exons[junctionExons[currExon]]

                    val = int(row[0])
                    if val > 0:
                        geneEnd = self.aligned.exons[junctionExons[currExon]] + offsetInExon
        return geneBounds


    def getUnsplicedGeneBounds(self, filehandle, start, end, overlapRadius, unsplicedIndex, startPos=0):
        segmentLength = end-start        
        sectionLen = 100000
        breakpointId = int(start / sectionLen)

        geneBounds = []

        for k,v in unsplicedIndex.items():
            filehandle.seek(startPos + v[breakpointId])

            offset = breakpointId * sectionLen
            line = filehandle.readline().rstrip()

            inRange = False
            geneStart = None
            zeroLen = 0
            geneEnd = None

            firstRange = True

            while offset < end and not line[0] == '>':
                row = line.split('\t')
                if len(row) == 1:
                    length = 1
                else:
                    length = int(row[1])

                if inRange:
                    if int(row[0]) == 0:
                        zeroLen += length
                        if zeroLen >= overlapRadius:
                            if firstRange:
                                firstRange = False
                            if not geneStart == None:
                                if geneEnd == None:
                                    geneBounds.append((geneStart, offset))
                                else:
                                    geneBounds.append((geneStart, geneEnd))
                            zeroLen = 0
                            geneStart = None
                            geneEnd = None
                        elif geneEnd == None:
                            geneEnd = offset
                    else:
                        zeroLen = 0
                        geneEnd = None
                        if geneStart == None and not firstRange:
                            geneStart = offset
                else:
                    if offset+length > start:
                        inRange = True

                        if int(row[0]) == 0:
                            if (offset+length-start) >= overlapRadius:
                                firstRange = False
                            else:
                                zeroLen = offset+length-start

                offset += length
                line = filehandle.readline().rstrip()
        return geneBounds

    def readExons(self, filehandle):
        line = filehandle.readline().rstrip()
        return [int(e) for e in line.split('\t')]

    def readHeader(self, filehandle):
        chromInfo = filehandle.readline().rstrip().split('\t')
        chromosomes = dict()
        for c in chromInfo:
            ch = c.split(',')
            chromosomes[ch[0]] = int(ch[1])
        return chromosomes

    # def readIndex(self, filehandle):
    #     '''
    #         All index lines begin with #
    #     '''

    #     line = filehandle.readline().rstrip()
    #     splicedIndex = [int(i) for i in line.split('\t')[1:]]

    #     unsplicedIndex = dict()
    #     while True:
    #         offset = filehandle.tell()
    #         line = filehandle.readline().rstrip()
    #         if not line[0] == '#':
    #             # reset pointer
    #             filehandle.seek(offset)
    #             return splicedIndex, unsplicedIndex
    #         else:
    #             row = line[1:].split('\t')

    #             unsplicedIndex[int(row[0])] = [int(i) for i in row[1:]]

    #     return splicedIndex, unsplicedIndex

    # def skipIndex(self, filehandle):
    #     '''
    #         All index lines begin with #
    #     '''
    #     while True:
    #         offset = filehandle.tell()
    #         if not filehandle.readline()[0] == '#':
    #             # reset pointer
    #             filehandle.seek(offset)
    #             return

    # def readIndex(self, filehandle):
    #     '''
    #         All index lines begin with #
    #     '''

    #     line = filehandle.readline().decode('ascii').rstrip()
    #     self.exonChunkSize = int(line[1:])
    #     line = filehandle.readline().decode('ascii').rstrip()
    #     self.junctionsIndex = [(i.split(',')[0], int(i.split(',')[1])) for i in line[1:].split(' ')]

    #     self.unsplicedIndex = dict()
    #     self.unsplicedExonsIndex = dict()
    #     i = 0
    #     while True:
    #         offset = filehandle.tell()
    #         line = filehandle.readline().decode('ascii').rstrip()
    #         if not line[0] == '#':
    #             # reset pointer
    #             filehandle.seek(offset)
    #             return
    #         else:
    #             row = line[1:].split('\t')

    #             if i == 0:
    #                 self.unsplicedIndex[int(row[0])] = [int(i) for i in row[1:]]
    #             else:
    #                 self.unsplicedExonsIndex[int(row[0])] = [int(i) for i in row[1:]]

    #             i = (i+1) % 2



