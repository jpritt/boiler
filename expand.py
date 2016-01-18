#! /usr/bin/env python3
import alignments
import read
import pairedread
import bucket
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
        self.bundles = binaryIO.readClusters(f)
        spliced_index = binaryIO.readListFromFile(f)

        self.expandCrossBundleBuckets(f)

        for i in range(len(self.bundles)):
            self.aligned.exons = self.bundles[i]

            #print('Expanding cluster')
            self.expandCluster(f, spliced_index[i])

            if i == 0:
                with open(out_name, 'w') as f2:
                    self.aligned.writeSAM(f2, True, self.force_xs)
            else:
                with open(out_name, 'a') as f2:
                    self.aligned.writeSAM(f2, False, self.force_xs)

            self.aligned.unpaired = []
            self.aligned.paired = []

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
            junc, startPos = binaryIO.readJunction(cluster, bucket.Bucket(exons, length, boundaries), readLenBytes, startPos)
            junc.strand = key[-2]
            junc.NH = key[-1]

            junc.coverage = self.RLEtoVector(junc.coverage)

            self.expandJunc(junc)

    def expandCrossBundleBuckets(self, filehandle):
        num_bundles = len(self.bundles)
        bundleIdBytes = binaryIO.findNumBytes(num_bundles)
        numBytes = binaryIO.readVal(filehandle, 1)
        length = binaryIO.readVal(filehandle, numBytes)

        if length > 0:
            index = self.expandString(filehandle.read(length))
            num_buckets, startPos = binaryIO.binaryToVal(index, 4, start=0)
            chunk_size, startPos = binaryIO.binaryToVal(index, 2, startPos)
            readLenBytes, startPos = binaryIO.binaryToVal(index, 1, startPos)
            buckets, startPos = binaryIO.readCrossBundleBucketNames(index, num_buckets, bundleIdBytes, startPos)
            chunk_lens, startPos = binaryIO.readList(index, startPos)

            i = 0
            for l in chunk_lens:
                chunk = self.expandString(filehandle.read(l))
                startPos = 0

                for i in range(i, min(i+chunk_size, num_buckets)):
                    b = buckets[i]
                    startPos = binaryIO.readCrossBundleBucket(chunk, b, readLenBytes, startPos)

                    b.coverage = self.RLEtoVector(b.coverage)

                    exonsA = self.bundles[b.bundleA]
                    exonsB = self.bundles[b.bundleB]
                    exon_bounds = [(exonsA[e], exonsA[e+1]) for e in b.exonIdsA] + [(exonsB[e], exonsB[e+1]) for e in b.exonIdsB]
                    b.exon_bounds = exon_bounds

                    boundaries = [exon_bounds[0][1]-exon_bounds[0][0]]
                    for n in range(1, len(exon_bounds)):
                        boundaries.append(boundaries[-1] + exon_bounds[n][1]-exon_bounds[n][0])
                    b.boundaries = boundaries
                    self.expandCrossBundleBucket(b)

                    # Delete this bucket to save space
                    buckets[i] = None
                i += 1

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
                junc.strand = key[-2]
                junc.NH = key[-1]

                junc.coverage = self.RLEtoVector(junc.coverage)

                self.expandJunc(junc)

                juncId += 1

        #print('Found %d junctions, looking for %d' % (juncId, len(sorted_junctions)))
        if not juncId == len(sorted_junctions):
            exit()
    '''

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

            self.aligned.unpaired.append(read.Read(self.aligned.getChromosome(readExons[0][0]), readExons, junc.strand, junc.NH))

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
                                                     self.aligned.getChromosome(readExonsB[0][0]), readExonsB, junc.strand, junc.NH))


    def expandCrossBundleBucket(self, bucket):
        if not sum([e[1]-e[0] for e in bucket.exon_bounds]) == bucket.length:
            print(bucket.exon_bounds)
            print(sum([e[1]-e[0] for e in bucket.exon_bounds]))
            print(bucket.coverage)
            print(len(bucket.coverage))
            print(bucket.length)
            exit()

        unpaired, paired, t1, t2 = self.aligned.findReads(dict(), bucket.pairedLens, bucket.lensLeft, bucket.lensRight, bucket.coverage, bucket.boundaries)

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

        # Offset of start of junction from beginning of chromosome
        bucket_offset = bucket.exon_bounds[0][0]

        # marks indices in junction coverage vector where exons are
        mapping = [0]
        # offsets of start of each exon in junction
        offsets = [0]
        for j in range(1,len(bucket.exon_bounds)):
            mapping.append(mapping[-1] + bucket.exon_bounds[j-1][1] - bucket.exon_bounds[j-1][0])
            offsets.append(bucket.exon_bounds[j][0] - bucket.exon_bounds[0][0])

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
                readExons.append( [start+bucket_offset, end+bucket_offset] )
            else:
                readExons.append( [start+bucket_offset, offsets[i]+mapping[i+1]-mapping[i]+bucket_offset] )

                for x in range(i+1,j):
                    readExons.append( [offsets[x]+bucket_offset, offsets[x]+mapping[x+1]-mapping[x]+bucket_offset] )

                readExons.append( [offsets[j]+bucket_offset, end+bucket_offset] )

            self.aligned.unpaired.append(read.Read(self.aligned.getChromosome(readExons[0][0]), readExons, bucket.strand, bucket.NH))

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
                readExonsA.append( [start+bucket_offset, end+bucket_offset] )
            else:
                readExonsA.append( [start+bucket_offset, offsets[i]+mapping[i+1]-mapping[i]+bucket_offset] )
                for x in range(i+1,j):
                    readExonsA.append( [offsets[x]+bucket_offset, offsets[x]+mapping[x+1]-mapping[x]+bucket_offset] )
                readExonsA.append( [offsets[j]+bucket_offset, end+bucket_offset] )

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
                readExonsB.append( [start+bucket_offset, end+bucket_offset] )
            else:
                readExonsB.append( [start+bucket_offset, offsets[i]+mapping[i+1]-mapping[i]+bucket_offset] )
                for x in range(i+1,j):
                    readExonsB.append( [offsets[x]+bucket_offset, offsets[x]+mapping[x+1]-mapping[x]+bucket_offset] )
                readExonsB.append( [offsets[j]+bucket_offset, end+bucket_offset] )

            self.aligned.paired.append(pairedread.PairedRead(self.aligned.getChromosome(readExonsA[0][0]), readExonsA,  \
                                                     self.aligned.getChromosome(readExonsB[0][0]), readExonsB, bucket.strand, bucket.NH))

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

    def getChromosomes(self, compressedFilename):
        with open(compressedFilename, 'rb') as f:
            return binaryIO.readChroms(f)

    def getGeneBounds(self, compressedFilename, chrom, start=None, end=None):
        with open(compressedFilename, 'rb') as f:
            chromosomes = binaryIO.readChroms(f)
            if start == None or end == None:
                start = 0
                end = chromosomes[chrom]

            offset = 0
            for k in sorted(chromosomes.keys()):
                if not k == chrom:
                    offset += chromosomes[k]
                else:
                    break
            start += offset
            end += offset

            clusters = binaryIO.readClusters(f)
            start_i, end_i = self.getRelevantClusters(clusters, start, end)
            return [(c[0]-offset, c[-1]-offset) for c in clusters[start_i:end_i]]

    def getCoverage(self, compressedFilename, chrom, start=None, end=None):
        #print('Getting coverage in %s: %d - %d' % (chrom, start, end))
        with open(compressedFilename, 'rb') as f:

            chromosomes = binaryIO.readChroms(f)
            self.aligned = alignments.Alignments(chromosomes)
            if start == None:
                start = 0
            else:
                start = max(start, 0)
            if end == None:
                end = chromosomes[chrom]
            else:
                end = min(end, chromosomes[chrom])

            start += self.aligned.chromOffsets[chrom]
            end += self.aligned.chromOffsets[chrom]

            coverage = [0.0] * (end-start)

            self.bundles = binaryIO.readClusters(f)
            start_i, end_i = self.getRelevantClusters(self.bundles, start, end)

            if start_i >= end_i:
                return coverage

            spliced_index = binaryIO.readListFromFile(f)

            st = time.time()
            coverage = self.getAllCrossBucketsCoverage(f, coverage, start_i, end_i, start, end)
            en = time.time()
            print('%f s to get cross buckets coverage' % (en-st))

            processT = 0.0

            f.seek(sum(spliced_index[:start_i]), 1)
            st = time.time()
            for i in range(start_i, end_i):
                self.aligned.exons = self.bundles[i]
                #print('Bundle %d - %d (%d)' % (bundles[i][0], bundles[i][-1], bundles[i][-1]-bundles[i][0]))
                coverage, t = self.getBundleCoverage(f, spliced_index[i], coverage, start, end)
                processT += t
            en = time.time()
            print('%f s to process bundles %d - %d' % (en-st, start_i, end_i))

        return coverage

    def getAllCrossBucketsCoverage(self, filehandle, coverage, start_i, end_i, range_start, range_end):
        num_bundles = len(self.bundles)
        bundleIdBytes = binaryIO.findNumBytes(num_bundles)
        numBytes = binaryIO.readVal(filehandle, 1)
        length = binaryIO.readVal(filehandle, numBytes)

        if length > 0:
            index = self.expandString(filehandle.read(length))
            num_buckets, startPos = binaryIO.binaryToVal(index, 4, start=0)
            chunk_size, startPos = binaryIO.binaryToVal(index, 2, startPos)
            readLenBytes, startPos = binaryIO.binaryToVal(index, 1, startPos)
            st = time.time()
            buckets, startPos = binaryIO.readCrossBundleBucketNames(index, num_buckets, bundleIdBytes, startPos)
            en = time.time()
            chunk_lens, startPos = binaryIO.readList(index, startPos)
            print('Parsing header time: %f s' % (en-st))

            print('%d cross-bundle buckets' % num_buckets)
            curr_bucket = 0
            skip = 0

            for l in chunk_lens:
                buckets_in_chunk = min(chunk_size, num_buckets-curr_bucket)
                relevant = [0] * buckets_in_chunk
                last_relevant = -1

                for i in range(buckets_in_chunk):
                    bundleA = buckets[i+curr_bucket].bundleA
                    bundleB = buckets[i+curr_bucket].bundleB
                    if (bundleA >= start_i and bundleA < end_i) or (bundleB >= start_i and bundleB < end_i):
                        relevant[i] = 1
                        last_relevant = i

                if last_relevant == -1:
                    skip += l
                else:
                    if skip > 0:
                        filehandle.seek(skip, 1)
                        skip = 0

                    chunk = self.expandString(filehandle.read(l))
                    startPos = 0

                    for i in range(last_relevant+1):
                        if relevant[i]:
                            b = buckets[i+curr_bucket]
                            startPos = binaryIO.readCrossBundleBucket(chunk, b, readLenBytes, startPos)
                            b.coverage = self.RLEtoVector(b.coverage)

                            exonsA = self.bundles[b.bundleA]
                            exonsB = self.bundles[b.bundleB]

                            # Is this necessary?
                            exon_bounds = [(exonsA[e], exonsA[e+1]) for e in b.exonIdsA] + [(exonsB[e], exonsB[e+1]) for e in b.exonIdsB]
                            boundaries = [0]
                            for n in range(len(exon_bounds)):
                                boundaries.append(boundaries[-1] + exon_bounds[n][1]-exon_bounds[n][0])

                            coverage = self.getBucketCoverage(b, coverage, range_start, range_end, exon_bounds, boundaries)
                        else:
                            startPos = binaryIO.skipCrossBundleBucket(chunk, readLenBytes, startPos)

                curr_bucket += buckets_in_chunk

            if skip > 0:
                filehandle.seek(skip, 1)

        return coverage

    def getBucketCoverage(self, bucket, coverage, range_start, range_end, subexon_bounds, boundaries):
        NH = float(bucket.NH)

        for i in range(len(boundaries)-1):
            if subexon_bounds[i][1] < range_start:
                continue
            elif subexon_bounds[i][0] >= range_end:
                break

            if subexon_bounds[i][0] <= range_start:
                subexon_start_id = range_start - subexon_bounds[i][0]
                range_start_id = 0
            else:
                subexon_start_id = 0
                range_start_id = subexon_bounds[i][0] - range_start

            if subexon_bounds[i][1] <= range_end:
                subexon_end_id = subexon_bounds[i][1] - subexon_bounds[i][0]
                range_end_id = subexon_bounds[i][1] - range_start
            else:
                subexon_end_id = range_end - subexon_bounds[i][0]
                range_end_id = range_end - range_start

            # Add junc_coverage[subexon_start_id:subexon_end_id] to coverage[range_start_id:range_end_id]
            range_id = range_start_id
            for j in range(boundaries[i]+subexon_start_id, boundaries[i]+subexon_end_id):
                c = bucket.coverage[j]
                if c > 0:
                    coverage[range_id] += c / NH
                range_id += 1

        return coverage

    def getBundleCoverage(self, f, length, coverage, range_start, range_end):
        bundle = self.expandString(f.read(length))
        startPos = 0
        readLenBytes, startPos = binaryIO.binaryToVal(bundle, 1, startPos)
        sorted_buckets, exonBytes, startPos = binaryIO.readJunctionsList(bundle, startPos)

        start_i, end_i = self.getRelevantExons(self.aligned.exons, range_start, range_end)

        t = 0.0

        for key in sorted_buckets:
            exons = key[:-2]

            relevant = False
            for e in exons:
                if e >= start_i and e < end_i:
                    relevant = True
                    break

            # If the junction does not overlap the target region, skip it
            if not relevant:
                startPos = binaryIO.skipJunction(bundle, readLenBytes, startPos)
                continue

            # Otherwise, expand this junction
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
            junc, startPos = binaryIO.readJunction(bundle, bucket.Bucket(exons, length, boundaries), readLenBytes, startPos)
            junc.NH = float(key[-1])
            junc.coverage = self.RLEtoVector(junc.coverage)

            t1 = time.time()

            subexon_bounds = []
            for j in junc.exons:
                subexon_bounds.append([self.aligned.exons[j], self.aligned.exons[j+1]])

            # marks indices in junction coverage vector where exons begin
            if len(boundaries) == 0:
                mapping = [0, subexon_bounds[-1][1] - subexon_bounds[-1][0]]
            else:
                mapping = [0] + boundaries + [boundaries[-1] + subexon_bounds[-1][1] - subexon_bounds[-1][0]]

            coverage = self.getBucketCoverage(junc, coverage, range_start, range_end, subexon_bounds, mapping)

            t2 = time.time()
            t += t2 - t1

        return coverage, t

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

    def getReads(self, compressedFilename, chrom, start=None, end=None):
        with open(compressedFilename, 'rb') as f:
            chromosomes = binaryIO.readChroms(f)
            self.aligned = alignments.Alignments(chromosomes)
            if start == None:
                start = 0
            else:
                start = max(start, 0)
            if end == None:
                end = chromosomes[chrom]
            else:
                end = min(end, chromosomes[chrom])

            start += self.aligned.chromOffsets[chrom]
            end += self.aligned.chromOffsets[chrom]

            self.bundles = binaryIO.readClusters(f)
            start_i, end_i = self.getRelevantClusters(self.bundles, start, end)

            if start_i >= end_i:
                return []

            spliced_index = binaryIO.readListFromFile(f)

            unpaired, paired = self.getAllCrossBucketsReads(f, start_i, end_i, start, end)

            f.seek(sum(spliced_index[:start_i]), 1)
            for i in range(start_i, end_i):
                self.aligned.exons = self.bundles[i]
                #print('Bundle %d - %d (%d)' % (bundles[i][0], bundles[i][-1], bundles[i][-1]-bundles[i][0]))
                self.getBundleReads(f, spliced_index[i], start, end, unpaired, paired)

        return unpaired, paired

    def getAllCrossBucketsReads(self, filehandle, start_i, end_i, range_start, range_end):
        num_bundles = len(self.bundles)
        bundleIdBytes = binaryIO.findNumBytes(num_bundles)
        numBytes = binaryIO.readVal(filehandle, 1)
        length = binaryIO.readVal(filehandle, numBytes)

        unpaired = []
        paired = []

        if length > 0:
            index = self.expandString(filehandle.read(length))
            num_buckets, startPos = binaryIO.binaryToVal(index, 4, start=0)
            chunk_size, startPos = binaryIO.binaryToVal(index, 2, startPos)
            readLenBytes, startPos = binaryIO.binaryToVal(index, 1, startPos)
            buckets, startPos = binaryIO.readCrossBundleBucketNames(index, num_buckets, bundleIdBytes, startPos)
            chunk_lens, startPos = binaryIO.readList(index, startPos)

            curr_bucket = 0
            skip = 0
            for l in chunk_lens:
                buckets_in_chunk = min(chunk_size, num_buckets-curr_bucket)
                relevant = [0] * buckets_in_chunk
                last_relevant = -1

                for i in range(buckets_in_chunk):
                    bundleA = buckets[i+curr_bucket].bundleA
                    bundleB = buckets[i+curr_bucket].bundleB
                    if (bundleA >= start_i and bundleA < end_i) or (bundleB >= start_i and bundleB < end_i):
                        relevant[i] = 1
                        last_relevant = i

                if last_relevant == -1:
                    skip += l
                else:
                    if skip > 0:
                        filehandle.seek(skip, 1)
                        skip = 0

                    chunk = self.expandString(filehandle.read(l))
                    startPos = 0

                    for i in range(last_relevant):

                        if relevant[i]:
                            b = buckets[i+curr_bucket]
                            startPos = binaryIO.readCrossBundleBucket(chunk, b, readLenBytes, startPos)
                            b.coverage = self.RLEtoVector(b.coverage)

                            exonsA = self.bundles[b.bundleA]
                            exonsB = self.bundles[b.bundleB]
                            exon_bounds = [(exonsA[e], exonsA[e+1]) for e in b.exonIdsA] + [(exonsB[e], exonsB[e+1]) for e in b.exonIdsB]
                            b.exon_bounds = exon_bounds

                            boundaries = [exon_bounds[0][1]-exon_bounds[0][0]]
                            for n in range(1, len(exon_bounds)):
                                boundaries.append(boundaries[-1] + exon_bounds[n][1]-exon_bounds[n][0])
                            b.boundaries = boundaries

                            self.getCrossBucketReads(b, range_start, range_end, unpaired, paired)

                        else:
                            startPos = binaryIO.skipCrossBundleBucket(chunk, readLenBytes, startPos)

                curr_bucket += buckets_in_chunk

            if skip > 0:
                filehandle.seek(skip, 1)

        return unpaired, paired

    def getCrossBucketReads(self, bucket, range_start, range_end, unpaired_reads, paired_reads):
        unpaired, paired, t1, t2 = self.aligned.findReads(dict(), bucket.pairedLens, bucket.lensLeft, bucket.lensRight, bucket.coverage, bucket.boundaries)

        # Offset of start of junction from beginning of chromosome
        bucket_offset = bucket.exon_bounds[0][0]

        # marks indices in junction coverage vector where exons are
        mapping = [0]
        # offsets of start of each exon in junction
        offsets = [0]
        for j in range(1,len(bucket.exon_bounds)):
            mapping.append(mapping[-1] + bucket.exon_bounds[j-1][1] - bucket.exon_bounds[j-1][0])
            offsets.append(bucket.exon_bounds[j][0] - bucket.exon_bounds[0][0])

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
                readExons.append( [start+bucket_offset, end+bucket_offset] )
            else:
                readExons.append( [start+bucket_offset, offsets[i]+mapping[i+1]-mapping[i]+bucket_offset] )

                for x in range(i+1,j):
                    readExons.append( [offsets[x]+bucket_offset, offsets[x]+mapping[x+1]-mapping[x]+bucket_offset] )

                readExons.append( [offsets[j]+bucket_offset, end+bucket_offset] )

            if self.readOverlapsRegion(readExons, range_start, range_end):
                unpaired_reads.append(read.Read(self.aligned.getChromosome(readExons[0][0]), readExons, bucket.strand, bucket.NH))

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
                readExonsA.append( [start+bucket_offset, end+bucket_offset] )
            else:
                readExonsA.append( [start+bucket_offset, offsets[i]+mapping[i+1]-mapping[i]+bucket_offset] )
                for x in range(i+1,j):
                    readExonsA.append( [offsets[x]+bucket_offset, offsets[x]+mapping[x+1]-mapping[x]+bucket_offset] )
                readExonsA.append( [offsets[j]+bucket_offset, end+bucket_offset] )

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
                readExonsB.append( [start+bucket_offset, end+bucket_offset] )
            else:
                readExonsB.append( [start+bucket_offset, offsets[i]+mapping[i+1]-mapping[i]+bucket_offset] )
                for x in range(i+1,j):
                    readExonsB.append( [offsets[x]+bucket_offset, offsets[x]+mapping[x+1]-mapping[x]+bucket_offset] )
                readExonsB.append( [offsets[j]+bucket_offset, end+bucket_offset] )

            if self.readOverlapsRegion(readExonsA, range_start, range_end) or self.readOverlapsRegion(readExonsB, range_start, range_end):
                paired_reads.append(pairedread.PairedRead(self.aligned.getChromosome(readExonsA[0][0]), readExonsA, self.aligned.getChromosome(readExonsB[0][0]), readExonsB, bucket.strand, bucket.NH))

    def getBundleReads(self, f, length, range_start, range_end, unpaired, paired):
        bundle = self.expandString(f.read(length))
        startPos = 0
        readLenBytes, startPos = binaryIO.binaryToVal(bundle, 1, startPos)
        sorted_buckets, exonBytes, startPos = binaryIO.readJunctionsList(bundle, startPos)

        start_i, end_i = self.getRelevantExons(self.aligned.exons, range_start, range_end)

        for key in sorted_buckets:
            exons = key[:-2]

            relevant = False
            for e in exons:
                if e >= start_i and e < end_i:
                    relevant = True
                    break

            # If the junction does not overlap the target region, skip it
            if not relevant:
                startPos = binaryIO.skipJunction(bundle, readLenBytes, startPos)
                continue

            # Otherwise, expand this junction
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
            b, startPos = binaryIO.readJunction(bundle, bucket.Bucket(exons, length, boundaries), readLenBytes, startPos)
            b.NH = float(key[-1])
            b.strand = key[-2]
            b.coverage = self.RLEtoVector(b.coverage)

            subexon_bounds = []
            for j in b.exons:
                subexon_bounds.append([self.aligned.exons[j], self.aligned.exons[j+1]])
            b.exon_bounds = subexon_bounds

            self.getBucketReads(b, range_start, range_end, unpaired, paired)

        return unpaired, paired

    def getBucketReads(self, b, range_start, range_end, unpaired_reads, paired_reads):
        unpaired, paired, t1, t2 = self.aligned.findReads(b.unpairedLens, b.pairedLens, b.lensLeft, b.lensRight, b.coverage, b.boundaries)

        juncBounds = []
        for j in b.exons:
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

            if self.readOverlapsRegion(readExons, range_start, range_end):
                unpaired_reads.append(read.Read(self.aligned.getChromosome(readExons[0][0]), readExons, b.strand, b.NH))

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

            if self.readOverlapsRegion(readExonsA, range_start, range_end) or self.readOverlapsRegion(readExonsB, range_start, range_end):
                paired_reads.append(pairedread.PairedRead(self.aligned.getChromosome(readExonsA[0][0]), readExonsA,  \
                                                     self.aligned.getChromosome(readExonsB[0][0]), readExonsB, b.strand, b.NH))


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

    def updateRLE(self, vector, start, end, val):
        '''
        Update the run-length encoded vector by adding val to each base in the range [start, end)
        :param vector:
        :param start:
        :param end:
        :param val:
        :return:
        '''

        length = end - start

        len_vec = len(vector)
        i = 0
        while i < len_vec:
            if start < vector[i][1]:
                break
            else:
                start -= vector[i][1]
                i += 1

        if i >= len_vec:
            return vector

        if start > 0:
            vector = vector[:i] + [[vector[i][0], start], [vector[i][0], vector[i][1]-start]] + vector[i+1:]
            i += 1
            len_vec += 1

        while i < len_vec:
            if length < vector[i][1]:
                break
            else:
                vector[i][0] += val
                length -= vector[i][1]
                i += 1

        if i < len_vec and length > 0:
            vector = vector[:i] + [[vector[i][0]+val, length], [vector[i][0], vector[i][1]-length]] + vector[i+1:]

        return vector
