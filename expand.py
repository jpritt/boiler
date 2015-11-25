#! /usr/bin/env python3
import alignments
import read
import pairedread
import junction
import time
import math
import binaryIO
import bisect
import resource

class Expander:
    aligned = None
    sectionLen = 100000    

    # 0 - zlib
    # 1 - lzma
    # 2 - bz2
    compressMethod = 0

    def __init__(self, force_xs=False):
        self.debug = False
        if self.compressMethod == 0:
            self.zlib = __import__('zlib')
        elif self.compressMethod == 1:
            self.lzma = __import__('lzma')
        elif self.compressMethod == 2:
            self.bz2 = __import__('bz2')

        self.readCounts = []
        self.readTimes = []
        self.pairCounts = []
        self.pairTimes = []
        self.pairCountsB = []
        self.pairTimesB = []

        self.force_xs = force_xs


    def expand(self, compressedFilename, uncompressedFilename, binary=False, debug=False):
        ''' Expand both spliced and unspliced alignments
        '''
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


    def expandByCluster(self, f, out_name):
        clusters = binaryIO.readClusters(f)
        spliced_index = binaryIO.readListFromFile(f)
        self.junctionChunkSize = binaryIO.readVal(f, 2)

        maxLen = 0

        for i in range(len(clusters)):
            self.aligned.exons = clusters[i]
            l = self.aligned.exons[-1] - self.aligned.exons[0]
            if l > maxLen:
                maxLen = l

            self.expandCluster(f, spliced_index[i])

            if i == 0:
                with open(out_name, 'w') as f2:
                    self.aligned.writeSAM(f2, True, self.force_xs)
            else:
                with open(out_name, 'a') as f2:
                    self.aligned.writeSAM(f2, False, self.force_xs)

            self.aligned.unpaired = []
            self.aligned.paired = []

    '''
    def expandCluster(self, f, length):
        cluster = self.expandString(f.read(length))
        readLenBytes, startPos = binaryIO.binaryToVal(cluster, 1, 0)
        sorted_junctions, exonBytes, startPos = binaryIO.readJunctionsList(cluster, startPos)

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
            junc, startPos = binaryIO.readJunction(cluster, junction.Junction(exons, length, boundaries), readLenBytes, startPos)
            junc.xs = key[-2]
            junc.NH = key[-1]

            junc.coverage = self.RLEtoVector(junc.coverage)

            self.expandJunc(junc)
    '''

    def expandCluster(self, f, length):
        index = self.expandString(f.read(length))
        chunkLens, startPos = binaryIO.readList(index, 0)
        sorted_junctions, exonBytes, startPos = binaryIO.readJunctionsList(index, startPos)
        readLenBytes, startPos = binaryIO.binaryToVal(index, 1, startPos)

        juncId = 0
        for chunkLen in chunkLens:
            chunk = self.expandString(f.read(chunkLen))
            pos = 0
            l = len(chunk)
            while pos < l:
                key = sorted_junctions[juncId]
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
                junc, pos = binaryIO.readJunction(chunk, junction.Junction(exons, length, boundaries), readLenBytes, pos)
                junc.xs = key[-2]
                junc.NH = key[-1]

                junc.coverage = self.RLEtoVector(junc.coverage)

                self.expandJunc(junc)

                juncId += 1

        #print('Found %d junctions, looking for %d' % (juncId, len(sorted_junctions)))
        if not juncId == len(sorted_junctions):
            exit()

    def expandJunc(self, junc):
        unpaired, paired, t1, t2 = self.aligned.findReads(junc.unpairedLens, junc.pairedLens, junc.lensLeft, junc.lensRight, junc.coverage, junc.boundaries, self.debug)


        numP = len(paired)
        numR = len(unpaired) + 2 * numP
        if numR >= len(self.readCounts):
            self.readCounts += [0] * (numR - len(self.readCounts) + 1)
            self.readTimes += [0] * (numR - len(self.readTimes) + 1)
        if numP >= len(self.pairCountsB):
            self.pairCountsB += [0] * (numP - len(self.pairCountsB) + 1)
            self.pairTimesB += [0] * (numP - len(self.pairTimesB) + 1)
        self.readCounts[numR] += 1
        self.readTimes[numR] += t1
        self.pairCountsB[numP] += 1
        self.pairTimesB[numP] += t2

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

            self.aligned.unpaired.append(read.Read(self.aligned.getChromosome(readExons[0][0]), readExons, junc.xs, junc.NH))

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

    def getGeneBounds(self, compressedFilename, chrom, start=None, end=None):
        with open(compressedFilename, 'rb') as f:
            chromosomes = binaryIO.readChroms(f)
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
            return [(c[0], c[-1]) for c in clusters[start_i:end_i]]

    def getCoverage(self, compressedFilename, chrom, start=None, end=None):
        print('Getting coverage in %s: %d - %d' % (chrom, start, end))
        with open(compressedFilename, 'rb') as f:
            chromosomes = binaryIO.readChroms(f)
            self.aligned = alignments.Alignments(chromosomes)
            if start == None or end == None:
                start = 0
                end = chromosomes[chrom]

            coverage = [0.0] * (end-start)

            for k in sorted(chromosomes.keys()):
                if not k == chrom:
                    start += chromosomes[k]
                    end += chromosomes[k]
                else:
                    break

            bundles = binaryIO.readClusters(f)
            start_i, end_i = self.getRelevantClusters(bundles, start, end)

            #print('%d -- %d' % (start_i, end_i))
            #for i in range(start_i, end_i):
            #    print(clusters[i])

            if start_i == end_i:
                return coverage

            spliced_index = binaryIO.readListFromFile(f)
            self.junctionChunkSize = binaryIO.readVal(f, 2)

            for i in range(start_i):
                index = self.expandString(f.read(spliced_index[i]))
                chunkLens, startPos = binaryIO.readList(index, 0)
                f.seek(sum(chunkLens),1)
            for i in range(start_i, end_i):
                self.aligned.exons = bundles[i]
                coverage = self.getBundleCoverage(f, spliced_index[i], coverage, start, end)

        return coverage

    def getBundleCoverage(self, f, length, coverage, range_start, range_end):
        index = self.expandString(f.read(length))
        chunkLens, startPos = binaryIO.readList(index, 0)
        sorted_junctions, exonBytes, startPos = binaryIO.readJunctionsList(index, startPos)
        readLenBytes, startPos = binaryIO.binaryToVal(index, 1, startPos)

        start_i, end_i = self.getRelevantExons(self.aligned.exons, range_start, range_end)

        skip = 0
        juncId = 0
        for i in range(len(chunkLens)):
            # Get the number of junctions in this chunk
            c_len = min(self.junctionChunkSize, len(sorted_junctions)-juncId)
            if i == len(chunkLens)-1:
                c_len = len(sorted_junctions) - i * self.junctionChunkSize

            # Only expand this chunk if one of the junctions overlaps the target region
            chunk = None

            startPos = 0

            for j in range(c_len):
                # Check if any of the subexons in this junction overlap the target region
                relevant = False
                for e in sorted_junctions[juncId+j][:-2]:
                    if e >= start_i and e < end_i:
                        relevant = True
                        break

                # If the junction does not overlap the target region and we have already expanded this chunk, skip it
                if not relevant:
                    if chunk:
                        startPos = binaryIO.skipJunction(chunk, readLenBytes, startPos)
                    continue

                # If the junction overlaps the target region and we have not expanded this chunk yet, expand it
                if not chunk:
                    if skip > 0:
                        f.seek(skip, 1)
                        skip = 0
                    chunk = self.expandString(f.read(chunkLens[i]))
                    startPos = 0

                    for _ in range(j):
                        startPos = binaryIO.skipJunction(chunk, readLenBytes, startPos)

                key = sorted_junctions[juncId + j]
                exons = key[:-2]

                # Boundaries contains the positions in the junction coverage vector where each new subexon begins
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
                junc, startPos = binaryIO.readJunction(chunk, junction.Junction(exons, length, boundaries), readLenBytes, startPos)
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
                if len(boundaries) == 0:
                    mapping = [0, juncBounds[-1][1] - juncBounds[-1][0]]
                else:
                    mapping = [0] + boundaries + [boundaries[-1] + juncBounds[-1][1] - juncBounds[-1][0]]

                #print(boundaries)
                #print(mapping)
                #print('')


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
                genome_pos = self.aligned.exons[exons[0]] - range_start
                length = range_end - range_start

                start = 0
                for b in range(len(boundaries)):
                    for i in range(start, boundaries[b]):
                        j = genome_pos+i
                        if j >= 0:
                            if j >= length:
                                break
                            else:
                                coverage[j] += junc_coverage[i]
                    genome_pos += juncBounds[b][0] - juncBounds[b-1][1]
                    start += boundaries[b]
                '''

            juncId += c_len

            if not chunk:
                skip += chunkLens[i]

        if skip > 0:
            f.seek(skip, 1)
        return coverage

    def getRelevantClusters(self, clusters, start, end):
        '''

        :param clusters:
        :param start:
        :param end:
        :return: The start and end+1 index in clusters of all clusters that overlap [start,end)
        '''

        cluster_bounds = [(c[0], c[-1]) for c in clusters]

        i = bisect.bisect_left(cluster_bounds, (start,start))

        if cluster_bounds[i-1][0] <= start and cluster_bounds[i-1][1] > start:
            start_i = i-1
        elif cluster_bounds[i][1] > start:
            start_i = i
        else:
            print('Error bisecting')
            print('(%d,%d)' % (start,end))
            print(cluster_bounds[i-3:i+3])
            exit()

        end_i = bisect.bisect_left(cluster_bounds, (end,end))

        return start_i, end_i

    def getRelevantExons(self, exons, start, end):
        start_i = bisect.bisect_right(exons, start) - 1
        end_i = bisect.bisect_left(exons, end)
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

                #self.getUnsplicedReads(f, unspliced_index[i], start, end)
                self.getSplicedReads(f, spliced_index[i], start, end)

        return self.aligned.unpaired, self.aligned.paired

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
                    self.aligned.unpaired.append(read.Read(self.aligned.getChromosome(readExons[0][0]), readExons, junc.xs, junc.NH))

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

