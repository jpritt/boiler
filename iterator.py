import binaryIO

class BIterator:

    def __init__(self, compressed):
        '''
        :param compressed: Path to Boiler-compressed file
        '''

        try:
            self.f = open(compressed, 'rb')
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))

        # Read index
        self.chromosomes = binaryIO.readChroms(self.f)
        self.bundles = binaryIO.readClusters(self.f)
        self.spliced_index = binaryIO.readListFromFile(f)

        self.curr_cluster = 0
        self.num_clusters = len(self.clusters)

        # Read index for cross-bundle buckets
        num_bundles = len(self.bundles)
        bundleIdBytes = binaryIO.findNumBytes(num_bundles)
        numBytes = binaryIO.readVal(self.f, 1)
        length = binaryIO.readVal(self.f, numBytes)
        index = self.expandString(self.f.read(length))

        self.num_buckets, startPos = binaryIO.binaryToVal(index, 4, start=0)
        self.cross_bundle_chunk_size, startPos = binaryIO.binaryToVal(index, 2, startPos)
        self.readLenBytes, startPos = binaryIO.binaryToVal(index, 1, startPos)
        self.cross_bundle_buckets, startPos = binaryIO.readCrossBundleBucketNames(index, self.num_buckets, bundleIdBytes, startPos)
        self.cross_bundle_chunk_lens, startPos = binaryIO.readList(index, startPos)

        self.cross_bundle_start = self.f.tell()

    def get_coverage(self):
        if self.curr_cluster >= self.num_clusters:
            return None

        if self.cross_buckets_length > 0:
            coverage = self.getAllCrossBucketsCoverage()

    def get_alignments(self):
        if self.curr_cluster >= self.num_clusters:
            return None

    def next(self):
        self.curr_cluster += 1
        if self.curr_cluster >= self.num_clusters:
            return None

        else:
            return self.clusters[self.curr_cluster]


    def getAllCrossBucketsCoverage(self):
        self.f.seek(self.cross_bundle_start)

        curr_bucket = 0
        skip = 0

        for l in self.cross_bundle_chunk_lens:
            buckets_in_chunk = min(self.cross_bundle_chunk_size, self.num_buckets-curr_bucket)
            relevant = [0] * buckets_in_chunk
            last_relevant = -1

            for i in range(buckets_in_chunk):
                bundleA = self.cross_bundle_buckets[i+curr_bucket].bundleA
                bundleB = self.cross_bundle_buckets[i+curr_bucket].bundleB
                if (bundleA == self.curr_cluster) or (bundleB == self.curr_cluster):
                    relevant[i] = 1
                    last_relevant = i

            if last_relevant == -1:
                skip += l
            else:
                if skip > 0:
                    self.f.seek(skip, 1)
                    skip = 0

                chunk = self.expandString(self.f.read(l))
                startPos = 0

                for i in range(last_relevant+1):
                    if relevant[i]:
                        b = self.cross_bundle_buckets[i+curr_bucket]
                        startPos = binaryIO.readCrossBundleBucket(chunk, b, self.readLenBytes, startPos)
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
                        startPos = binaryIO.skipCrossBundleBucket(chunk, self.readLenBytes, startPos)

            curr_bucket += buckets_in_chunk

        if skip > 0:
            self.f.seek(skip, 1)

        return coverage

    def getBucketCoverage(self, bucket, coverage, subexon_bounds, boundaries):
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


    def expandString(self, s):
        ''' Use a predefined python library to expand the given string.
            Return the decompressed string '''

        if self.compressMethod == 0:
            return self.zlib.decompress(s)
        elif self.compressMethod == 1:
            return self.lzma.decompress(s)
        elif self.compressMethod == 2:
            return self.bz2.decompress(s)