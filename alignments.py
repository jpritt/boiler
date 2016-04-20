import bisect
import bucket
import cross_bundle_bucket
import pairedread
import copy
import random
import time

class Alignments:
    ''' A set of reads aligned to a genome '''

    def __init__(self, chromosomes, frag_len_cutoff=None, split_discordant=True):
        ''' Initialize a genome for alignments

            chromosomes: A dictionary with keys corresponding to chromosome names and values corresponding to lengths
        '''

        self.frag_len_cutoff = frag_len_cutoff
        self.split_discordant = split_discordant

        self.chromosomes = chromosomes
        self.chromosomeNames = sorted(chromosomes.keys())

        # Initialize exon breaks between all chromosomes
        self.exons = set()

        # Offset of each chromosome from the start of the genome
        self.chromOffsets = dict()

        # Bases are generally indexed from 1
        nextOffset = 1
        for c in self.chromosomeNames:
            self.chromOffsets[c] = nextOffset
            nextOffset += chromosomes[c]+1
            #self.exons += [nextOffset]


        # List of potential gene boundaries as tuples
        self.gene_bounds = []
        # If 2 reads are less than this far apart, combine them into a single gene
        self.overlap_radius = 50

        # paired reads for which the mate still needs to be found
        self.unmatched = dict()

        # Unmatched reads that span multiple bundles
        self.cross_bundle_reads = dict()

        # Unmatched reads *in the current bundle* that span multiple bundles
        self.curr_cross_bundle_reads = dict()

        # Matched mates spanning multiple bundles
        self.cross_bundle_pairs = []

        self.cross_bundle_buckets = dict()

        self.max_cross_bundle_read_len = 0

        self.unpaired = []
        self.paired = []

        self.read_timeA = 0.0
        self.read_countA = 0
        self.read_timeB = 0.0
        self.read_countB = 0

    def processRead(self, name, read, paired):
        ''' If read is unpaired, add it to the correct spliced or unspliced list of reads.
            If read is paired, find its pair or add it to a list to be found later. Once a pair of reads is found, add the combined read to the appropriate list of reads
        '''

        # Update read location for chromosome
        offset = self.chromOffsets[read.chrom]
        for i in range(len(read.exons)):
            read.exons[i] = [read.exons[i][0]+offset, read.exons[i][1]+offset]

        # update list of subexon bounds
        alignment = read.exons
        if len(alignment) > 1:
            for i in range(len(alignment)-1):
                self.exons.add(alignment[i][1])
                self.exons.add(alignment[i+1][0])

        # Update the boundaries of the current bundle
        #self.update_gene_bounds(read.exons[0][0], read.exons[-1][1])

        if not paired:
            # Update the boundaries of the current bundle
            self.update_gene_bounds(read.exons[0][0], read.exons[-1][1])

            self.add_unpaired(read)
        else:   # paired read
            # update pair location for chromsome
            read.pairOffset += self.chromOffsets[read.pairChrom]

            self.update_gene_bounds(read.exons[0][0], read.exons[-1][1])

            foundMate = False
            if read.pairOffset <= read.exons[0][0]:

                i = self.find_mate(read, name, self.unmatched)
                if i >= 0:
                    mate = self.unmatched[name][i]
                    if mate.exons[-1][1] > read.exons[-1][1]:
                        self.add_paired(read, mate)
                    else:
                        self.add_paired(mate, read)

                    foundMate = True
                    del self.unmatched[name][i]

                if not foundMate:
                    i = self.find_mate(read, name, self.cross_bundle_reads)
                    if i >= 0:
                        mate = self.cross_bundle_reads[name][i]
                        self.cross_bundle_pairs.append((mate, read))

                        foundMate = True
                        del self.cross_bundle_reads[name][i]
            if not foundMate:
                # Mate has not been processed yet

                if (read.pairOffset - read.exons[0][0]) < self.frag_len_cutoff:
                    self.update_gene_bounds(read.exons[0][0], read.pairOffset)
                

                if name in self.unmatched:
                    self.unmatched[name].append(read)
                else:
                    self.unmatched[name] = [read]
                
    def find_mate(self, read, name, unmatched):
        '''
        Search the list of unmatched reads for one matching the given name and location
        :param read: Read information including location and mate location
        :param name: Identifying name (first column of SAM information)
        :param unmatched: Dictionary of unmatched reads with key = name, value = list of reads
        :return: Index in unmatched of matching read, or -1 if no match found
        '''

        if not name in unmatched:
            return -1

        for i in range(len(unmatched[name])):
            match = unmatched[name][i]
            if read.pairOffset == match.exons[0][0] and match.pairOffset == read.exons[0][0] and (not self.split_discordant or not self.conflicts(read.exons, match.exons)):
                # Return index in unmatched dictionary of match
                return i

        # If no match found, return -1
        return -1


    def update_gene_bounds(self, start, end):
        '''
        Update the boundaries of the current bundle to include [start, end]
        '''
        if not self.gene_bounds:
            self.gene_bounds = [start, end]
        else:
            if start < self.gene_bounds[0]:
                self.gene_bounds[0] = start
            if end > self.gene_bounds[1]:
                self.gene_bounds[1] = end

    def conflicts(self, exonsA, exonsB):
        '''

        :param exonsA: List containing the exon bounds for gene A in the form [(x_0,y_0), (x_1,y_1),...]
        :param exonsB: List containing the exon bounds for gene B in the same form as exonsA
        :return: 1 if an exon in gene A overlaps an intron in gene B, 2 if vice versa, 3 if one gene range lies strictly inside the other, 0 otherwise.
        '''

        if (exonsA[0][0] < exonsB[0][0] and exonsA[-1][1] > exonsB[-1][1]) or (exonsB[0][0] < exonsA[0][0] and exonsB[-1][1] > exonsA[-1][1]):
            # One set of exons contains the other
            return 3

        for e in exonsB:
            if e[0] > exonsA[-1][0]:
                break

            for i in range(len(exonsA)-1):
                if e[0] >= exonsA[-i-1][0]:
                    break
                elif e[1] > exonsA[-i-2][1]:
                    # Exon in B overlaps an intron in A
                    return 1

        countA = len(exonsA)
        for i in range(countA):
            e = exonsA[countA-i-1]
            if e[1] < exonsB[0][1]:
                break

            for i in range(len(exonsB)-1):
                if e[1] <= exonsB[i][1]:
                    break
                elif e[1] > exonsB[i][1] and e[0] < exonsB[i+1][0]:
                    # Exon in A overlaps an intron in B
                    return 2

        return 0

    def add_unpaired(self, read):
        '''
        Add this read to the list of unpaired
        '''
        self.unpaired.append(read)

    def add_paired(self, read1, read2):
        '''
        Create a single paired read from these two reads and add it to the paired list
        '''


        strand = read1.strand or read2.strand
        NH = min(read1.NH, read2.NH)
        p = pairedread.PairedRead(read1.chrom, read1.exons, read2.chrom, read2.exons, strand, NH)

        if not self.split_discordant and read1.chrom == read2.chrom and self.conflicts(read1.exons, read2.exons):
            p.discordant = True
            #if not read1.name == read2.name:
            #    print('Names %s, %s do not match' % (read1.name, read2.name))
            #    exit()
            #print(read1.name)
            #exit()
        else:
            p.discordant = False

        self.paired.append(p)

    def finalizeUnmatched(self):
        '''
        Finalize unmatched (discordant) reads. We convert them to unpaired reads.
        '''
        for name,reads in self.unmatched.items():
            if reads:
                for r in reads:
                    if hasattr(r,'pairOffset') and r.pairOffset > self.gene_bounds[1]:
                        # Waiting for mate in a future bundle

                        if not hasattr(r, 'exonIds'):
                            r.exonIds, r.length = self.getExonIds(r.exons)
                            r.bucket_length = sum([self.exons[e+1]-self.exons[e] for e in r.exonIds])
                            r.startOffset = r.exons[0][0] - self.exons[r.exonIds[0]]
                            r.endOffset = self.exons[r.exonIds[-1]+1] - r.exons[-1][1]

                        if name in self.cross_bundle_reads:
                            self.cross_bundle_reads[name].append(r)
                        else:
                            self.cross_bundle_reads[name] = [r]
                    else:
                        self.add_unpaired(r)

        # Reset dictionary for next bundle
        self.unmatched = dict()

    def finalizeExons(self):
        '''
        Convert the set of exon boundaries to a list
        '''

        self.exons.add(self.gene_bounds[0])
        self.exons.add(self.gene_bounds[1])
        self.exons = sorted(list(self.exons))

    def finalize_cross_bundle_reads(self):
        '''
        Process the list of reads with mates outside this bundle
        '''
        # Finalize cross-bundle pairs that were discovered in this bundle
        for p in self.cross_bundle_pairs:
            if not hasattr(p[0], 'exonIds'):
                #p[0].bundle = bundle_id
                p[0].exonIds, p[0].length = self.getExonIds(p[0].exons)
                p[0].bucket_length = sum([self.exons[e+1]-self.exons[e] for e in p[0].exonIds])
                p[0].startOffset = p[0].exons[0][0] - self.exons[p[0].exonIds[0]]
                p[0].endOffset = self.exons[p[0].exonIds[-1]+1] - p[0].exons[-1][1]
            if not hasattr(p[1], 'exonIds'):
                #p[1].bundle = bundle_id
                p[1].exonIds, p[1].length = self.getExonIds(p[1].exons)
                p[1].bucket_length = sum([self.exons[e+1]-self.exons[e] for e in p[1].exonIds])
                p[1].startOffset = p[1].exons[0][0] - self.exons[p[1].exonIds[0]]
                p[1].endOffset = self.exons[p[1].exonIds[-1]+1] - p[1].exons[-1][1]


            # Update maximum length of a read in a cross-bucket bundle (needed to efficiently store length distribution)
            if p[0].length > self.max_cross_bundle_read_len:
                self.max_cross_bundle_read_len = p[0].length
            if p[1].length > self.max_cross_bundle_read_len:
                self.max_cross_bundle_read_len = p[1].length

            NH = min(p[0].NH, p[1].NH)
            strand = p[0].strand or p[1].strand
            if strand == '-':
                strand = -1
            elif strand == '+':
                strand = 1
            else:
                strand = 0
            key = str(p[0].bundle) + ' ' + ' '.join([str(e) for e in p[0].exonIds]) + ',' + str(p[1].bundle) + ' ' + ' '.join([str(e) for e in p[1].exonIds]) + ' ' + str(NH) + ' ' + str(strand)
            if key in self.cross_bundle_buckets:
                self.cross_bundle_buckets[key].add_pair(p[0], p[1])
            else:
                b = cross_bundle_bucket.CrossBundleBucket(p[0].bundle, p[0].exonIds, p[1].bundle, p[1].exonIds)
                b.set_length(p[0].bucket_length + p[1].bucket_length)
                b.NH = NH
                b.strand = strand
                b.add_pair(p[0], p[1])
                self.cross_bundle_buckets[key] = b

    def computeBuckets(self):
        # Compute coverage levels across every bucket
        partitions = dict()

        # Maximum fragment length
        max_len = 0

        for r in self.unpaired:
            partitions = self.add_unpaired_to_partition(r, partitions)
            if r.length > max_len:
                max_len = r.length
        for p in self.paired:
            partitions = self.add_paired_to_partition(p, partitions)
            if p.length > max_len:
                max_len = p.length

        return partitions, max_len

    def add_paired_to_partition(self, p, partitions):
        '''

        :param p: Paired read to add
        :param partitions: List of previously computed partitions
        :return: Updated partitions list
        '''

        self.finalize_paired_read(p)

        b = self.find_key(p, partitions)
        b.add_paired(p)

        return partitions

    def add_unpaired_to_partition(self, r, partitions):
        '''
        :param r: Read to add
        :param partitions: List of previously computed partitions
        :return: Updated partitions list
        '''

        self.finalize_unpaired_read(r)

        b = self.find_key(r, partitions)
        b.add_unpaired(r)

        return partitions

    def find_key(self, r, partitions):
        '''

        :param r: Either a Read or PairedRead object
        :param partitions: List of partitions
        :return: key and partition for this read
        '''

        exonIds = r.exonIds

        strand = '0'
        if r.strand == '-':
            strand = '-1'
        elif r.strand == '+':
            strand = '1'
        key = '\t'.join([str(e) for e in exonIds]) + '\t' + strand + '\t' + str(r.NH)

        if not key in partitions:
            covLength = 0
            for i in range(len(exonIds)):
                e = exonIds[i]
                subexon_length = self.exons[e+1] - self.exons[e]

                covLength += subexon_length

            partitions[key] = bucket.Bucket(exonIds, covLength)
            partitions[key].strand = r.strand
            partitions[key].NH = r.NH

            partitions[key].countBefore = 0

        j = partitions[key]

        return j

    def finalize_paired_read(self, p):
        # TODO: Take this out
        #if p.exonsB[0][0] < p.exonsA[0][0] or p.exonsA[-1][1] > p.exonsB[-1][1]:
        #    print('Weird reads:')
        #    print(p.exonsA)
        #    print(p.exonsB)
        #    print('')

        # Find what range of subexon ids contains this pair
        #for i in range(len(self.exons)-1):
        #    if self.exons[i+1] > p.exonsA[0][0]:
        #        break
        #for j in range(i, len(self.exons)):
        #    if self.exons[j] > p.exonsB[-1][1]:
        #        break

        # Exon ids spanned by read
        start_id = bisect.bisect_right(self.exons, p.exonsA[0][0]) - 1
        id = start_id
        p.startOffset = p.exonsA[0][0] - self.exons[id]
        p.lenLeft = 0
        exonIdsA = []
        for exon in p.exonsA:
            p.lenLeft += exon[1] - exon[0]
            while self.exons[id+1] <= exon[0]:
                id += 1
            while self.exons[id] < exon[1]:
                exonIdsA.append(id)
                id += 1

        id = start_id
        p.lenRight = 0

        exonIdsB = []
        for exon in p.exonsB:
            p.lenRight += exon[1] - exon[0]
            while self.exons[id+1] <= exon[0]:
                id += 1
            while self.exons[id] < exon[1]:
                exonIdsB.append(id)
                id += 1

        p.endOffset = self.exons[id] - p.exonsB[-1][1]

        if p.discordant:
            p.exonIds = exonIdsA + exonIdsB
        else:
            n = 0
            while n < len(exonIdsA) and exonIdsA[n] < exonIdsB[0]:
                n += 1
            p.exonIds = exonIdsA[:n] + exonIdsB

    def finalize_unpaired_read(self, r):
        # Exon ids spanned by read
        r.exonIds, r.length = self.getExonIds(r.exons)
        r.startOffset = r.exons[0][0] - self.exons[r.exonIds[0]]
        r.endOffset = self.exons[r.exonIds[-1]+1] - r.exons[-1][1]

    def getExonIds(self, exons):
        exonIds = []

        id = bisect.bisect_right(self.exons, exons[0][0]) - 1
        length = 0

        for exon in exons:
            length += exon[1] - exon[0]
            while self.exons[id+1] <= exon[0]:
                id += 1
            while self.exons[id] < exon[1]:
                exonIds.append(id)
                id += 1

        return exonIds, length

    def RLE(self, vector):
        rle = []

        val = vector[0]
        length = 0
        for v in vector:
            if v == val:
                length += 1
            else:
                rle.append([val, length])
                val = v
                length = 1

        rle.append([val, length])

        return rle

    def resetBundle(self):
        self.exons = set()
        self.gene_bounds = []
        self.unpaired = []
        self.paired = []
        self.unmatched = dict()
        self.cross_bundle_pairs = []
        self.curr_cross_bundle_reads = dict()



    def findReads(self, unpairedLens, pairedLens, lensLeft, lensRight, coverage, boundaries=None):
        ''' Find the set of reads that most closely matches the distribution of readLens and the coverage vector
        '''

        fragmentLens = copy.copy(unpairedLens)
        numFragments = 0
        for k,v in unpairedLens.items():
            numFragments += v
        for k,v in lensLeft.items():
            if k > 0:
                numFragments += v
                if k in fragmentLens:
                    fragmentLens[k] += v
                else:
                    fragmentLens[k] = v
        for k,v in lensRight.items():
            if k > 0:
                numFragments += v
                if k in fragmentLens:
                    fragmentLens[k] += v
                else:
                    fragmentLens[k] = v

        countPaired = 0
        for k,v in pairedLens.items():
            countPaired += v

        t1 = time.time()

        # If the number of potential start sites is equal to the number of reads, then we can solve it with high accuracy with v1. Otherwise use v3
        starts = 0

        lens_orig = dict()
        for k,v in fragmentLens.items():
            lens_orig[k] = v

        for i in range(len(coverage)):
            if i == 0 and coverage[i] > 0:
                starts += coverage[i]
            elif coverage[i] > coverage[i-1]:
                starts += coverage[i+1] - coverage[i]
        if starts == numFragments:
            ta = time.time()
            reads = self.findReadsInCoverage_v1(coverage, fragmentLens)
            tb = time.time()
            self.read_timeA += (tb - ta)
            self.read_countA += 1
        else:
            ta = time.time()
            reads = self.findReadsInCoverage_v3(coverage, fragmentLens)
            tb = time.time()
            self.read_timeB += (tb - ta)
            self.read_countB += 1
        t2 = time.time()

        if boundaries:
            unpaired, paired = self.findPairsWithBoundaries(reads, unpairedLens, pairedLens, boundaries)
        else:
            #unpaired, paired = self.findPairs(reads, pairedLens)
            unpaired, paired = self.findPairsRandom(reads, pairedLens)
        t3 = time.time()

        return unpaired, paired, t2-t1, t3-t2

    def printTime(self):
        print('Simple version:')
        print('  Total time: %f s' % self.read_timeA)
        print('  Num segments: %d' % self.read_countA)
        print('  Avg time: %f s' % (self.read_timeA / self.read_countA))
        print('')
        print('Full version:')
        print('  Total time: %f s' % self.read_timeB)
        print('  Num segments: %d' % self.read_countB)
        print('  Avg time: %f s' % (self.read_timeB / self.read_countB))
        print('')

    def findPairsGreedy(self, reads, pairedLens):
        pairedLensSorted = sorted(pairedLens.keys(), reverse=True)

        reads.sort()
        numReads = len(reads)

        # Keep track of which reads we have already assigned
        assigned = [0] * numReads

        # Map each distance to the paired reads that match it
        pairDists = dict()
        singleDists = dict()

        # Create a distance matrix between all pairs of reads
        dists = [0] * numReads
        for i in range(numReads):
            dists[i] = [0] * (i+1)
            for j in range(i):
                d = reads[i][1] - reads[j][0]
                dists[i][j] = d
                if d in pairDists:
                    pairDists[d] += 1
                else:
                    pairDists[d] = 1

            d = reads[i][1] - reads[i][0]
            dists[i][i] = d
            if d in singleDists:
                singleDists[d] += 1
            else:
                singleDists[d] = 1

        paired = []
        countPaired = 0
        for k,v in pairedLens.items():
            countPaired += v
        while countPaired > 0:
            bestFreq = None
            bestL = None
            foundPair = False

            # Look for unique paired read lengths
            for l in pairedLensSorted:
                if not (l in pairDists and pairDists[l] > 0 and pairedLens[l] > 0):
                    continue

                expected = pairedLens[l]
                freq = pairDists[l]
                if freq == 0:
                    continue
                if freq <= expected:
                    pairedLens[l] -= 1
                    i,j = self.findPairedDist(dists, assigned, numReads, l)
                    paired.append([reads[j], reads[i]])

                    assigned[i] = 1
                    assigned[j] = 1

                    if dists[i][i] > 0:
                        singleDists[dists[i][i]] -= 1
                        dists[i][i] = 0
                    if dists[j][j] > 0:
                        singleDists[dists[j][j]] -= 1
                        dists[j][j] = 0

                    for x in range(i):
                        if dists[i][x] > 0:
                            pairDists[dists[i][x]] -= 1
                            dists[i][x] = 0
                    for x in range(i+1,numReads):
                        if dists[x][i] > 0:
                            pairDists[dists[x][i]] -= 1
                            dists[x][i] = 0
                    for x in range(j):
                        if dists[j][x] > 0:
                            pairDists[dists[j][x]] -= 1
                            dists[j][x] = 0
                    for x in range(j+1,numReads):
                        if dists[x][j] > 0:
                            pairDists[dists[x][j]] -= 1
                            dists[x][j] = 0

                    foundPair = True
                    countPaired -= 1
                    break
                elif bestFreq == None or (freq-expected) < bestFreq:
                    bestFreq = freq - expected
                    bestL = l

            # No unique paired lengths, so choose one from the least frequent
            if not foundPair:
                if bestFreq == None:
                    break
                else:
                    pairedLens[bestL] -= 1
                    i,j = self.findPairedDist(dists, assigned, numReads, bestL)
                    paired.append([reads[j], reads[i]])

                    assigned[i] = 1
                    assigned[j] = 1

                    if dists[i][i] > 0:
                        singleDists[dists[i][i]] -= 1
                        dists[i][i] = 0
                    if dists[j][j] > 0:
                        singleDists[dists[j][j]] -= 1
                        dists[j][j] = 0

                    for x in range(i):
                        if dists[i][x] > 0:
                            pairDists[dists[i][x]] -= 1
                            dists[i][x] = 0
                    for x in range(i+1,numReads):
                        if dists[x][i] > 0:
                            pairDists[dists[x][i]] -= 1
                            dists[x][i] = 0
                    for x in range(j):
                        if dists[j][x] > 0:
                            pairDists[dists[j][x]] -= 1
                            dists[j][x] = 0
                    for x in range(j+1,numReads):
                        if dists[x][j] > 0:
                            pairDists[dists[x][j]] -= 1
                            dists[x][j] = 0

                    countPaired -= 1

        remaining = [0] * (numReads - sum(assigned))
        i = 0
        for j in range(numReads):
            if not assigned[j]:
                remaining[i] = reads[j]
                i += 1

        remainingPairedLens = dict()
        for k,v in pairedLens.items():
            if v > 0:
                remainingPairedLens[k] = v

        if countPaired > 0:
            newUnpaired, newPaired = self.findPairsGreedy2(remaining, remainingPairedLens)
            unpaired = newUnpaired
            paired += newPaired
        else:
            unpaired = remaining

        return unpaired, paired

    def findPairsGreedy2(self, reads, pairedLens):
        numReads = len(reads)
        reads.sort()
        pairedLensSorted = sorted(pairedLens.keys(), reverse=True)

        paired = []

        countPaired = 0
        for k,v in pairedLens.items():
            countPaired += v
        countUnpaired = numReads - 2 * countPaired

        # Create a distance matrix between all pairs of reads
        dists = [0] * numReads
        for i in range(numReads):
            dists[i] = [0] * i
            for j in range(i):
                d = reads[i][1] - reads[j][0]
                dists[i][j] = d
        assigned = [0] * numReads

        lenIndex = 0
        while countPaired > 0:
            targetL = pairedLensSorted[lenIndex]

            bestL = None
            bestDiff = None
            bestPos = None

            for i in range(numReads):
                for j in range(i,0,-1):
                    l = dists[i][j-1]
                    diff = abs(l - targetL)
                    if l > 0 and (bestDiff == None or diff < bestDiff):
                        bestDiff = diff
                        bestL = dists[i][j-1]
                        bestPos = (j-1, i)
                    elif l > targetL:
                        break

            if bestL == None:
                break
            else:
                pairedLens[targetL] -= 1
                if pairedLens[targetL] == 0:
                    lenIndex += 1

                i = bestPos[0]
                j = bestPos[1]
                paired.append([reads[i], reads[j]])
                assigned[i] = 1
                assigned[j] = 1

                for x in range(i):
                    dists[i][x] = 0
                for x in range(i+1,numReads):
                    dists[x][i] = 0
                for x in range(j):
                    dists[j][x] = 0
                for x in range(j+1,numReads):
                    dists[x][j] = 0

                countPaired -= 1

        unpaired = [0] * countUnpaired
        i = 0
        for j in range(numReads):
            if not assigned[j]:
                unpaired[i] = reads[j]
                i += 1

        return unpaired, paired

    def findPairedDist(self, dists, assigned, numReads, value):
        for i in range(numReads):
            if assigned[i]:
                continue
            for j in range(i):
                if assigned[j]:
                    continue
                if dists[i][j] == value:
                    return i,j

    def findPairsWithBoundaries(self, reads, unpaired_lens, paired_lens, boundaries):
        '''
        Use the fact that all reads must span all the subexons to improve our pairing

        :param reads:
        :param unpaired_lens:
        :param paired_lens:
        :param lens_left:
        :param lens_right:
        :param boundaries: Boundaries between subexons
        :return:
        '''

        numUnpaired = 0
        for k,v in unpaired_lens.items():
            numUnpaired += v
        numPaired = 0
        for k,v in paired_lens.items():
            numPaired += v

        left_reads = []
        spanning_reads = []
        right_reads = []

        for r in reads:
            if r[0] < boundaries[0]:
                if r[1] > boundaries[-1]:
                    spanning_reads.append(r)
                else:
                    left_reads.append(r)
            elif r[1] > boundaries[-1]:
                right_reads.append(r)
            else:
                spanning_reads.append(r)
                #print('Read does not overlap left or right...?')
                #exit()

        left_reads.sort()
        spanning_reads.sort()
        right_reads.sort()

        if left_reads:
            unique = [left_reads[0]]
            left_counts = [1]
            for r in left_reads[1:]:
                if r == unique[-1]:
                    left_counts[-1] += 1
                else:
                    unique.append(r)
                    left_counts.append(1)
            left_reads = unique
            left_bounds = [0, len(left_reads)]
        else:
            left_counts = []
            left_bounds = [0, 0]


        if spanning_reads:
            unique = [spanning_reads[0]]
            spanning_counts = [1]
            for r in spanning_reads[1:]:
                if r == unique[-1]:
                    spanning_counts[-1] += 1
                else:
                    unique.append(r)
                    spanning_counts.append(1)
            spanning_reads = unique
            spanning_bounds = [0, len(spanning_reads)]
        else:
            spanning_counts = []
            spanning_bounds = [0, 0]

        if right_reads:
            unique = [right_reads[0]]
            right_counts = [1]
            for r in right_reads[1:]:
                if r == unique[-1]:
                    right_counts[-1] += 1
                else:
                    unique.append(r)
                    right_counts.append(1)
            right_reads = unique
            right_bounds = [0, len(right_reads)]
        else:
            right_counts = []
            right_bounds = [0, 0]

        count_left = sum(left_counts)
        count_right = sum(right_counts)
        count_spanning = sum(spanning_counts)

        countUnpaired = 0
        for k,v in unpaired_lens.items():
            countUnpaired += v

        countPaired = 0
        for k,v in paired_lens.items():
            countPaired += v

        paired = []
        unpaired = None
        unmatched = []
        while countPaired > 0:
            if count_spanning == 0 and (count_left == 0 or count_right == 0):
                break

            if count_left == 0 and count_right == 0:
                if countPaired % 2 == 0:
                    p = self.findLeftPairRandom(spanning_reads, paired_lens, spanning_counts, spanning_bounds)
                else:
                    p = self.findRightPairRandom(spanning_reads, paired_lens, spanning_counts, spanning_bounds)
                if len(p) == 2:
                    paired.append([spanning_reads[p[0]][:], spanning_reads[p[1]][:]])
                    countPaired -= 1
                    count_spanning -= 2
                else:
                    unmatched.append(spanning_reads[p[0]][:])
                    count_spanning -= 1

            elif count_left >= count_right:
                p = self.findLeftPair(left_reads, spanning_reads, right_reads, left_counts, spanning_counts, right_counts, left_bounds, spanning_bounds, right_bounds, paired_lens)
                count_left -= 1
                if len(p) == 3:
                    if p[1] < 0:
                        paired.append([left_reads[p[0]][:], right_reads[p[2]][:]])
                        count_right -= 1
                    else:
                        paired.append([left_reads[p[0]][:], spanning_reads[p[1]][:]])
                        count_spanning -= 1
                    countPaired -= 1
                else:
                    unmatched.append(left_reads[p[0]][:])
            else:
                p = self.findRightPair(left_reads, spanning_reads, right_reads, left_counts, spanning_counts, right_counts, left_bounds, spanning_bounds, right_bounds, paired_lens)
                count_right -= 1
                if len(p) == 3:
                    if p[1] < 0:
                        paired.append([left_reads[p[0]][:], right_reads[p[2]][:]])
                        count_left -= 1
                    else:
                        paired.append([spanning_reads[p[1]][:], right_reads[p[2]][:]])
                        count_spanning -= 1
                    countPaired -= 1
                else:
                    unmatched.append(right_reads[p[0]][:])

        for i in range(left_bounds[0], left_bounds[1]):
            for _ in range(left_counts[i]):
                unmatched.append([left_reads[i][0], left_reads[i][1]])
        for i in range(spanning_bounds[0], spanning_bounds[1]):
            for _ in range(spanning_counts[i]):
                unmatched.append([spanning_reads[i][0], spanning_reads[i][1]])
        for i in range(right_bounds[0], right_bounds[1]):
            for _ in range(right_counts[i]):
                unmatched.append([right_reads[i][0], right_reads[i][1]])
        unmatched.sort()

        i = 0
        j = len(unmatched)
        for _ in range(countPaired):
            if i < j-1:
                paired.append([unmatched[i], unmatched[j-1]])
            i += 1
            j -= 1

        if i < j:
            if unpaired:
                unpaired += unmatched[i:j]
            else:
                unpaired = unmatched[i:j]
        elif not unpaired:
            unpaired = []

        return unpaired, paired

    def findLeftPair(self, left_reads, spanning_reads, right_reads, left_counts, spanning_counts, right_counts, left_bounds, spanning_bounds, right_bounds, paired_lens):
        '''
        :param left_reads:
        :param spanning_reads:
        :param right_reads:
        :param paired_lens:
        :return:
        '''

        i = left_bounds[0]
        start = left_reads[i][0]
        left_counts[left_bounds[0]] -= 1
        while left_bounds[0] < left_bounds[1] and left_counts[left_bounds[0]] == 0:
            left_bounds[0] += 1

        # Look for a match among right reads
        for j in range(right_bounds[0],right_bounds[1]):
            if right_counts[j] == 0:
                continue

            l = right_reads[j][1] - start
            if l in paired_lens:
                if paired_lens[l] > 1:
                    paired_lens[l] -= 1
                else:
                    del paired_lens[l]
                right_counts[j] -= 1
                while right_bounds[1] > right_bounds[0] and right_counts[right_bounds[1]-1] == 0:
                    right_bounds[1] -= 1
                while right_bounds[0] < right_bounds[1] and right_counts[right_bounds[0]] == 0:
                    right_bounds[0] += 1

                return [i, -1, j]

        # Look for a match among left reads
        for j in range(spanning_bounds[0],spanning_bounds[1]):
            if spanning_counts[j] == 0:
                continue

            l = spanning_reads[j][1] - start
            if l in paired_lens:
                if paired_lens[l] > 1:
                    paired_lens[l] -= 1
                else:
                    del paired_lens[l]
                spanning_counts[j] -= 1
                while spanning_bounds[1] > spanning_bounds[0] and spanning_counts[spanning_bounds[1]-1] == 0:
                    spanning_bounds[1] -= 1
                while spanning_bounds[0] < spanning_bounds[1] and spanning_counts[spanning_bounds[0]] == 0:
                    spanning_bounds[0] += 1

                return [i, j, -1]

        return [i]

    def findRightPair(self, left_reads, spanning_reads, right_reads, left_counts, spanning_counts, right_counts, left_bounds, spanning_bounds, right_bounds, paired_lens):
        '''
        :param left_reads:
        :param spanning_reads:
        :param right_reads:
        :param paired_lens:
        :return:
        '''

        j = right_bounds[1]-1
        end = right_reads[j][1]
        right_counts[right_bounds[1]-1] -= 1
        while right_bounds[0] < right_bounds[1] and right_counts[right_bounds[1]-1] == 0:
            right_bounds[1] -= 1

        # Look for a match among left reads
        for i in range(left_bounds[0],left_bounds[1]):
            if left_counts[i] == 0:
                continue

            l = end - left_reads[i][1]
            if l in paired_lens:
                if paired_lens[l] > 1:
                    paired_lens[l] -= 1
                else:
                    del paired_lens[l]
                left_counts[i] -= 1
                while left_bounds[1] > left_bounds[0] and left_counts[left_bounds[1]-1] == 0:
                    left_bounds[1] -= 1
                while left_bounds[0] < left_bounds[1] and left_counts[left_bounds[0]] == 0:
                    left_bounds[0] += 1

                return [i, -1, j]

        # Look for a match among left reads
        for i in range(spanning_bounds[0],spanning_bounds[1]):
            if spanning_counts[i] == 0:
                continue

            l = end - spanning_reads[i][0]
            if l in paired_lens:
                if paired_lens[l] > 1:
                    paired_lens[l] -= 1
                else:
                    del paired_lens[l]
                spanning_counts[i] -= 1
                while spanning_bounds[1] > spanning_bounds[0] and spanning_counts[spanning_bounds[1]-1] == 0:
                    spanning_bounds[1] -= 1
                while spanning_bounds[0] < spanning_bounds[1] and spanning_counts[spanning_bounds[0]] == 0:
                    spanning_bounds[0] += 1

                return [-1, i, j]

        return [j]


    def findPairsRandom(self, reads, paired_lens, debug=False):
        #print('Pairing %d reads' % len(reads))

        countPairs = 0
        for k,v in paired_lens.items():
            countPairs += v

        reads.sort()

        unique_reads = [reads[0]]
        read_counts = [1]
        for r in reads[1:]:
            if r == unique_reads[-1]:
                read_counts[-1] += 1
            else:
                unique_reads.append(r)
                read_counts.append(1)

        # Index of first and last reads in array that have not been used yet
        read_bounds = [0, len(read_counts)]

        paired = []
        unmatched = []
        unmatched_counts = []
        while countPairs > 0 and read_bounds[0] < read_bounds[1]:
            p = self.findLeftPairRandom(unique_reads, paired_lens, read_counts, read_bounds)
            if len(p) == 2:
                paired.append([unique_reads[p[0]][:], unique_reads[p[1]][:]])
                countPairs -= 1
            else:
                self.add_to_unmatched(unmatched, unmatched_counts, unique_reads[p[0]], 1)

            if countPairs == 0 or read_bounds[0] >= read_bounds[1]:
                break

            p = self.findRightPairRandom(unique_reads, paired_lens, read_counts, read_bounds)
            if debug:
                print(p)
            if len(p) == 2:
                paired.append([unique_reads[p[0]][:], unique_reads[p[1]][:]])
                countPairs -= 1
            else:
                self.add_to_unmatched(unmatched, unmatched_counts, unique_reads[p[0]], 1)

        # Add remaining reads to unmatched
        for i in range(read_bounds[0], read_bounds[1]):
            if read_counts[i] > 0:
                self.add_to_unmatched(unmatched, unmatched_counts, unique_reads[i], read_counts[i])

        num_remaining = sum(unmatched_counts)
        bounds = [0, len(unmatched)]

        paired_lens_sorted = sorted(paired_lens)

        while countPairs > 0 and num_remaining > 1:
            p = self.findClosestLeftPair(unmatched, unmatched_counts, bounds, paired_lens_sorted)
            paired.append(p)
            countPairs -= 1
            num_remaining -= 2

            if countPairs == 0 or num_remaining < 2:
                break

            p = self.findClosestRightPair(unmatched, unmatched_counts, bounds, paired_lens_sorted)
            paired.append(p)
            countPairs -= 1
            num_remaining -= 2

        unpaired = [0] * sum(unmatched_counts)
        id = 0
        for i in range(bounds[0], bounds[1]):
            for _ in range(unmatched_counts[i]):
                unpaired[id] = [unmatched[i][0], unmatched[i][1]]
                id += 1

        #print('Done!')
        return unpaired, paired

    def add_to_unmatched(self, unmatched, counts, read, num):
        i = bisect.bisect_left(unmatched, read)
        if i < len(unmatched) and unmatched[i] == read:
            counts[i] += num
        else:
            unmatched.insert(i, read)
            counts.insert(i, num)


    def findLeftPairRandom(self, reads, paired_lens, read_counts, read_bounds):
        i = read_bounds[0]
        read_counts[i] -= 1

        while read_bounds[0] < read_bounds[1] and read_counts[read_bounds[0]] == 0:
            read_bounds[0] += 1

        for j in range(read_bounds[1]-1, read_bounds[0]-1, -1):
            if read_counts[j] == 0:
                continue

            l = reads[j][1] - reads[i][0]
            if l in paired_lens:
                if paired_lens[l] > 1:
                    paired_lens[l] -= 1
                else:
                    del paired_lens[l]
                read_counts[j] -= 1
                while read_bounds[1] > read_bounds[0] and read_counts[read_bounds[1]-1] == 0:
                    read_bounds[1] -= 1
                while read_bounds[0] < read_bounds[1] and read_counts[read_bounds[0]] == 0:
                    read_bounds[0] += 1

                return [i, j]
        return [i]

    def findRightPairRandom(self, reads, paired_lens, read_counts, read_bounds):
        j = read_bounds[1]-1
        read_counts[j] -= 1

        while read_bounds[1] > read_bounds[0] and read_counts[read_bounds[1]-1] == 0:
            read_bounds[1] -= 1

        for i in range(read_bounds[0], read_bounds[1]):
            if read_counts[i] == 0:
                continue

            l = reads[j][1] - reads[i][0]
            if l in paired_lens:
                if paired_lens[l] > 1:
                    paired_lens[l] -= 1
                else:
                    del paired_lens[l]
                read_counts[i] -= 1
                while read_bounds[0] < read_bounds[1] and read_counts[read_bounds[0]] == 0:
                    read_bounds[0] += 1
                while read_bounds[1] > read_bounds[0] and read_counts[read_bounds[1]-1] == 0:
                    read_bounds[1] -= 1

                return [i,j]
        return [j]

    def findClosestLeftPair(self, reads, counts, bounds, paired_lens_sorted):
        i = bounds[0]
        start = reads[i][0]
        counts[i] -= 1
        while bounds[0] < bounds[1] and counts[bounds[0]] == 0:
            bounds[0] += 1

        num_lens = len(paired_lens_sorted)

        # Distance to closest match
        closestD = None
        # Index of closest match
        closestJ = None
        # Length of closest match
        closestL = None
        for j in range(bounds[0], bounds[1]):
            if counts[j] == 0:
                continue

            l = reads[j][1] - start
            id = bisect.bisect_left(paired_lens_sorted, l)
            if id < num_lens:
                d = abs(l - paired_lens_sorted[id])
                if d == 0:
                    closestD = 0
                    closestJ = j
                    closestL = paired_lens_sorted[id]
                    break
                elif closestD == None or d < closestD:
                    closestD = d
                    closestJ = j
                    closestL = paired_lens_sorted[id]
            if id > 0:
                d = abs(l - paired_lens_sorted[id-1])
                if closestD == None or d < closestD:
                    closestD = d
                    closestJ = j
                    closestL = paired_lens_sorted[id-1]

        if closestD == None:
            print('Error! Trying to pair only 1 read?')
            exit()

        pair = [reads[i][:], reads[closestJ][:]]
        counts[closestJ] -= 1
        while bounds[0] < bounds[1] and counts[bounds[0]] == 0:
            bounds[0] += 1
        while bounds[0] < bounds[1] and counts[bounds[1]-1] == 0:
            bounds[1] -= 1

        return pair

    def findClosestRightPair(self, reads, counts, bounds, paired_lens_sorted):
        j = bounds[1]-1
        end = reads[j][1]
        counts[j] -= 1
        while bounds[0] < bounds[1] and counts[bounds[1]-1] == 0:
            bounds[1] -= 1

        num_lens = len(paired_lens_sorted)

        # Distance to closest match
        closestD = None
        # Index of closest match
        closestI = None
        # Length of closest match
        closestL = None
        for i in range(bounds[1]-1, bounds[0]-1, -1):
            if counts[i] == 0:
                continue

            l = end - reads[i][1]
            id = bisect.bisect_left(paired_lens_sorted, l)
            if id < num_lens:
                d = abs(l - paired_lens_sorted[id])
                if d == 0:
                    closestD = 0
                    closestI = i
                    closestL = paired_lens_sorted[id]
                    break
                elif closestD == None or d < closestD:
                    closestD = d
                    closestI = i
                    closestL = paired_lens_sorted[id]
            if id > 0:
                d = abs(l - paired_lens_sorted[id-1])
                if closestD == None or d < closestD:
                    closestD = d
                    closestI = i
                    closestL = paired_lens_sorted[id-1]


        if closestD == None:
            print('Error! Trying to pair only 1 read?')
            exit()

        pair = [reads[closestI][:], reads[j][:]]
        counts[closestI] -= 1
        while bounds[0] < bounds[1] and counts[bounds[0]] == 0:
            bounds[0] += 1
        while bounds[0] < bounds[1] and counts[bounds[1]-1] == 0:
            bounds[1] -= 1

        return pair


    def findPairs(self, reads, pairedLens):
        if len(pairedLens) == 0:
            return reads, []

        length = 0
        for r in reads:
            if r[1] > length:
                length = r[1]

        paired = []
        unpaired = []

        reads.sort()

        # Sort read pairedLensSorted lengths
        pairedLensSorted = sorted(pairedLens, reverse=True)

        starts = [0] * (length+1)
        ends = [0] * (length+1)

        for r in reads:
            starts[r[0]] += 1
            ends[r[1]] += 1

        unmatched = []

        i = 0
        j = length

        while i < j and starts[i] <= 0:
            i += 1
        while j > i and ends[j] <= 0:
            j -= 1

        while i < j and len(pairedLensSorted) > 0:
            id1 = 0
            while not reads[id1][0] == i:
                id1 += 1

            startRead = reads[id1]

            starts[i] -= 1
            ends[reads[id1][1]] -= 1
            del reads[id1]

            foundRead = False
            for r in range(len(pairedLensSorted)):
                l = pairedLensSorted[r]
                if i+l <= length and ends[i+l] > 0:
                    pairedLens[l] -= 1
                    if pairedLens[l] == 0:
                        del pairedLensSorted[r]
                    foundRead = True

                    # Add paired read to list
                    id2 = 0
                    while not reads[id2][1] == i+l:
                        id2 += 1

                    starts[reads[id2][0]] -= 1
                    ends[i+l] -= 1


                    if reads[id2][0] < startRead[0]:
                        paired += [[reads[id2], startRead]]
                    else:
                        paired += [[startRead, reads[id2]]]

                    #paired += [[startRead, reads[id2]]]
                    del reads[id2]

                    break

            if not foundRead:
                unmatched += [startRead]

            while i < j and starts[i] <= 0:
                i += 1
            while j > i and ends[j] <= 0:
                j -= 1

            if j > i and len(pairedLensSorted) > 0:
                id1 = 0
                while not reads[id1][1] == j:
                    id1 += 1

                startRead = reads[id1]

                starts[reads[id1][0]] -= 1
                ends[j] -= 1
                del reads[id1]

                foundRead = False
                for r in range(len(pairedLensSorted)):
                    l = pairedLensSorted[r]
                    if j-l >= 0 and starts[j-l] > 0:
                        pairedLens[l] -= 1
                        if pairedLens[l] == 0:
                            del pairedLensSorted[r]
                        foundRead = True

                        # Add paired read to list
                        id2 = 0
                        while not reads[id2][0] == j-l:
                            id2 += 1

                        starts[j-l] -= 1
                        ends[reads[id2][1]] -= 1


                        if reads[id2][0] < startRead[0]:
                            paired += [[reads[id2], startRead]]
                        else:
                            paired += [[startRead, reads[id2]]]

                        #paired += [[reads[id2], startRead]]
                        del reads[id2]
                        break

                if not foundRead:
                    unmatched += [startRead]

                while i < j and starts[i] <= 0:
                    i += 1
                while j > i and ends[j] <= 0:
                        j -= 1

        # Pair up remaining reads until we meet the quota of paired-end reads

        paired_lens_new = dict()
        for k,v in pairedLens.items():
            if v > 0:
                paired_lens_new[k] = v

        u, p = self.findPairsGreedy(unmatched+reads, paired_lens_new)

        unpaired = u
        paired += p

        return unpaired, paired

    def findReadsInCoverage_v1(self, coverage, readLens, boundaries=None):
        ''' Given a coverage vector, return a set of reads the closely fits the coverage vector as well the distribution of read lengths.
            The algorithm creates new reads greedily, starting from both ends at once to ensure the ends of the vector are well marked.
            This algorithm is guaranteed to return a set of reads that covers every base at least to the corresponding depth of the coverage vector.
            In many cases the algorithm will overcompensate by creating extra reads to make sure every base in the coverage vector is covered.
            In such cases new reads have length equal to the median read length in the input distribution.

            coverage: Coverage vector containing the accumulation of many reads
            readLens: Dictionary containing the distribution of all read lengths in the coverage vector
            boundaries: Start and end points, of which all reads must cross at least 1
        '''

        lens = readLens.keys()

        # Find max and mode read lengths
        maxLen = max(lens)

        # Read lengths sorted by frequency, largest to smallest
        lensSorted = sorted(readLens, key=readLens.get, reverse=True)

        reads = []

        # start and end of coverage window
        # Keep finding reads from both ends until they meet in the middle
        start = 0
        while coverage[start] <= 0:
            start += 1
        end = len(coverage)
        while coverage[end-1] <= 0:
            end -= 1

        while end > start:

            # find a read from the beginning
            readStart = start
            readEnd = start

            closestEndpoint = None
            if boundaries and readStart >= boundaries[0]:
                minLen = boundaries[1]+1 - start
            else:
                minLen = 1
            for length in range(minLen, maxLen+1):
                if (readStart+length == end) or (readStart+length < end and coverage[readStart + length] < coverage[readStart + length - 1]):
                    if length in readLens:
                        readEnd = readStart + length
                        reads.append([readStart, readEnd])

                        readLens[length] -= 1

                        # reorder sorted lengths
                        for i in range(len(lensSorted)):
                            if lensSorted[i] == length:
                                break
                        j = i+1
                        while j < len(readLens) and readLens[lensSorted[j]] > readLens[lensSorted[i]]:
                            j += 1
                        if j > i+1:
                            lensSorted = lensSorted[:i] + lensSorted[i+1:j] + [lensSorted[i]] + lensSorted[j:]

                        if readLens[length] == 0:
                            del readLens[length]

                        break
                    else:
                        if closestEndpoint == None:
                            closestEndpoint = readStart + length

                # Don't extend into section where coverage is 0
                if (readStart+length) >= len(coverage) or coverage[readStart+length] == 0:
                    break
            if readEnd == readStart:
                if closestEndpoint == None:
                    lenId = 0
                    length = lensSorted[lenId]
                    readEnd = readStart + length

                    while readEnd > len(coverage) and lenId < (len(lensSorted)-1):
                        lenId += 1
                        length = lensSorted[lenId]
                        readEnd = readStart + length
                    if readEnd > len(coverage):
                        # No read lengths fit within the end of the vector
                        readEnd = len(coverage)

                    reads.append([readStart, readEnd])

                    if length in readLens:
                        readLens[length] -= 1

                        # reorder sorted lengths
                        for i in range(len(lensSorted)):
                            if lensSorted[i] == length:
                                break
                        j = i+1
                        while j < len(readLens) and readLens[lensSorted[j]] > readLens[lensSorted[i]]:
                            j += 1
                        if j > i+1:
                            lensSorted = lensSorted[:i] + lensSorted[i+1:j] + [lensSorted[i]] + lensSorted[j:]

                        if readLens[length] == 0:
                            del readLens[length]
                else:
                    readEnd = closestEndpoint
                    reads.append([readStart, readEnd])

            # Update coverage vector
            for i in range(readStart, readEnd):
                coverage[i] -= 1

            # update start
            while start < end and coverage[start] <= 0:
                start += 1
            while end > start and coverage[end-1] <= 0:
                end -= 1


            if end > start:
                readEnd = end
                readStart = end

                closestEndpoint = None
                if boundaries and readEnd <= boundaries[1]:
                    minLen = readEnd+1 - boundaries[0]
                else:
                    minLen = 1
                for length in range(minLen, maxLen+1):
                    if (readEnd-length == start) or (readEnd-length > start and coverage[readEnd - length] > coverage[readEnd - length - 1]):
                        if length in readLens:
                            readStart = readEnd - length
                            reads.append([readStart, readEnd])

                            readLens[length] -= 1


                            # reorder sorted lengths
                            for i in range(len(lensSorted)):
                                if lensSorted[i] == length:
                                    break
                            j = i+1

                            while j < len(readLens) and readLens[lensSorted[j]] > readLens[lensSorted[i]]:
                                j += 1
                            if j > i+1:
                                lensSorted = lensSorted[:i] + lensSorted[i+1:j] + [lensSorted[i]] + lensSorted[j:]


                            if readLens[length] == 0:
                                del readLens[length]

                            break
                        else:
                            if closestEndpoint == None:
                                closestEndpoint = readEnd - length

                if readStart == readEnd:
                    if closestEndpoint == None:
                        length = lensSorted[0]
                        readStart = readEnd - length
                        reads.append([readStart, readEnd])

                        if length in readLens:
                            readLens[length] -= 1

                            # reorder sorted lengths
                            for i in range(len(lensSorted)):
                                if lensSorted[i] == length:
                                    break

                            j = i+1
                            while j < len(readLens) and readLens[lensSorted[j]] > readLens[lensSorted[i]]:
                                j += 1
                            if j > i+1:
                                lensSorted = lensSorted[:i] + lensSorted[i+1:j] + [lensSorted[i]] + lensSorted[j:]

                            if readLens[length] == 0:
                                del readLens[length]
                    else:
                        readStart = closestEndpoint
                        reads.append([readStart, readEnd])

                for i in range(readStart, readEnd):
                    coverage[i] -= 1

                # update end
                while coverage[end-1] <= 0 and end > start:
                    end -= 1
                while coverage[start] <= 0 and start < end:
                    start += 1
        return reads

    def findReadsInCoverage_v3(self, cov, frag_lens):
        reads = []
        num_reads = 0

        start = 0
        end = len(cov)
        while start < end and cov[start] <= 0:
            start += 1
        while start < end and cov[end-1] <= 0:
            end -= 1

        min_len = None
        max_len = None
        mode_len = 0
        mode_freq = 0
        for l,f in frag_lens.items():
            if not min_len or l < min_len:
                min_len = l
            if not max_len or l > max_len:
                max_len = l
            if f > mode_freq:
                mode_freq = f
                mode_len = l

        while start < end:
            start, end, num_reads = self.findReadLeft(cov, start, end, reads, min_len, max_len, mode_len, num_reads)

            #if start < end:
            #    start, end = self.findReadRight(cov, start, end, reads, min_len, max_len, mode_len)

        reads = self.finalizeCovReads(reads, min_len, max_len, mode_len)

        return reads

    def finalizeCovReads(self, reads, min_len, max_len, mode_len):
        '''
        Split sections of end-to-end reads and return the list of all individual read boundaries

        :param reads:
        :param min_len:
        :param max_len:
        :param mode_len:
        :return:
        '''

        final = []
        for r in reads:
            if r[2] == 1:
                final.append([r[0], r[1]])
            else:
                pos = r[0]
                for i in range(r[2]):
                    total_len = r[1] - pos
                    min_left = total_len - max_len * (r[2]-i-1)
                    max_left = total_len - min_len * (r[2]-i-1)

                    l = max(min(mode_len, max_left), min_left)
                    final.append([pos, pos+l])
                    pos += l
        return final

    def findReadsEndingAt(self, reads, end):
        r = []
        for i in range(len(reads)-1, -1, -1):
            if reads[i][1] == end:
                r.append(i)
        return r

    def findReadsStartingAt(self, reads, start):
        r = []
        for i in range(len(reads)-1, -1, -1):
            if reads[i][0] == start:
                r.append(i)
        return r

    def findNextGap(self, cov, start, end, modeLen):
        gap = False
        for i in range(min(modeLen,end-start)):
            if cov[start+i] > 0:
                if gap:
                    return [gap_start, i]
            else:
                if not gap:
                    gap_start = i
                    gap = True
        if gap:
            return [gap_start, i+1]
        else:
            return None

    def extendReadRight(self, cov, start, stop, reads, max_len, num_reads):
        '''
        Extend the right end of a read ending at 'start' up to 'stop'
        '''

        if num_reads == 0:
            return False

        j = 0
        for i in range(num_reads-1, -1, -1):
            if reads[i][1] < stop:
                j = i
                break

        for i in range(j, -1, -1):
            r = reads[i]
            if r[1] < start:
                return False
            elif r[1] == start and (r[1] - r[0]) < (max_len * r[2]):
                r[1] = min(max_len*r[2] + r[0], stop)

                if r[1] < stop:
                    # Recalculate new position for read
                    for new_pos in range(j, i-1, -1):
                        if reads[new_pos][1] < r[1]:
                            break
                else:
                    new_pos = j

                if i < new_pos:
                    self.move_read(reads, i, new_pos)

                for i in range(start, r[1]):
                    cov[i] -= 1

                return True
        return False

    def shortenReadRight(self, cov, start, stop, reads, max_shorten, min_len, num_reads):
        '''
        Shorten the right end of a read ending at 'stop' down to 'start'
        '''

        if num_reads == 0:
            return False

        for i in range(num_reads-1, -1, -1):
            r = reads[i]
            if r[1] < stop:
                return False
            elif (start - r[0]) >= (min_len * r[2]):
                new_pos = 0
                for j in range(i, -1, -1):
                    if reads[j][1] < start:
                        new_pos = j+1
                        break

                for n in range(start, r[1]):
                    cov[n] += 1
                r[1] = start

                self.move_read(reads, i, new_pos)

                return True
        return False

    def add_read(self, reads, new_read, num_reads):
        if num_reads == 0:
            reads.append([new_read[0], new_read[1], 1])
            return 1

        # After adding new read, this should be its position in the array
        new_pos = -1
        for i in range(num_reads-1, -1, -1):
            if reads[i][1] < new_read[1]:
                new_pos = i
                break

        for i in range(new_pos, -1, -1):
            if reads[i][1] == new_read[0]:
                # We can extend this read rather than adding a new one
                reads[i][1] = new_read[1]
                reads[i][2] += 1

                # Move read from i to j
                if i < new_pos:
                    self.move_read(reads, i, new_pos)
                return num_reads

            elif reads[i][1] < new_read[0]:
                break

        # We have to add a new read
        r = [new_read[0], new_read[1], 1]
        reads.append(r)

        # Move read from end to j
        if new_pos+1 < num_reads:
            self.move_read(reads, num_reads, new_pos+1)

        return num_reads+1

    def move_read(self, reads, i, j):
        r = reads[i]
        if i < j:
            for x in range(i, j):
                reads[x] = reads[x+1]
            reads[j] = r
        elif i > j:
            for x in range(i, j, -1):
                reads[x] = reads[x-1]
            reads[j] = r

    def findReadLeft(self, cov, start, end, reads, min_len, max_len, mode_len, num_reads):
        gap = self.findNextGap(cov, start, end, mode_len)

        len_range = max_len - min_len

        if not gap:
            readEnd = min(start+mode_len,end)

            if readEnd-start < min_len:
                if not self.extendReadRight(cov, start, readEnd, reads, max_len, num_reads):
                    # We want to make sure that coverage ends at the right place
                    num_reads = self.add_read(reads, [readEnd - min_len, readEnd], num_reads)
                    for i in range(start, readEnd):
                        cov[i] -= 1
            else:
                #print('Adding read [%d, %d]' % (start, readEnd))
                num_reads = self.add_read(reads, [start, readEnd], num_reads)
                for i in range(start, readEnd):
                    cov[i] -= 1
        else:
            r = self.extendReadRight(cov, start, start+gap[0], reads, max_len, num_reads)
            if not r:
                r = self.shortenReadRight(cov, start+gap[0], start+gap[1], reads, len_range, min_len, num_reads)

            if not r:
                if gap[0] >= min_len:
                    readEnd = start+gap[0]
                    num_reads = self.add_read(reads, [start, readEnd], num_reads)
                elif gap[0] < min_len/2:
                    # Skip bases
                    readEnd = start+gap[0]
                else:
                    readEnd = min(start+mode_len, end)
                    num_reads = self.add_read(reads, [start, readEnd], num_reads)

                for i in range(start, readEnd):
                    cov[i] -= 1

        while start < end and cov[start] <= 0:
            start += 1
        while start < end and cov[end-1] <= 0:
            end -= 1

        return start, end, num_reads

    def getChromosome(self, index):
        ''' Return chromosome name containing the given index from the whole-genome vector
        '''

        for c in self.chromosomeNames:
            if self.chromosomes[c] > index:
                return c
            else:
                index -= self.chromosomes[c]

    def writeSAM(self, filehandle, unpaired, paired, header=True, force_xs=False, readId=0):
        ''' Write all alignments to a SAM file
        '''

        # write header
        if header:
            filehandle.write('@HD\tVN:1.0\tSO:unsorted\n')
            for k,v in self.chromosomes.items():
                filehandle.write('@SQ\tSN:' + str(k) + '\tLN:' + str(v) + '\n')

        for read in unpaired:
            exons = read.exons
            cigar = [str(exons[0][1] - exons[0][0]) + 'M']
            spliced = False
            for i in range(1, len(exons)):
                if exons[i][0] == exons[i-1][1]:
                    prevLen = int(cigar[-1][:-1])
                    cigar[-1] = str(prevLen + exons[i][1] - exons[i][0]) + 'M'
                else:
                    spliced = True
                    cigar += [str(exons[i][0] - exons[i-1][1]) + 'N']
                    cigar += [str(exons[i][1] - exons[i][0]) + 'M']
            cigar = ''.join(cigar)

            chrom = read.chrom
            offset = self.chromOffsets[chrom]

            if force_xs and spliced and not read.strand:
                #print('Assigning random XS value to spliced unpaired read')
                if random.randint(0,1) == 0:
                    read.strand = '+'
                else:
                    read.strand = '-'

            if read.strand:
                filehandle.write(chrom+':'+str(readId) + '\t0\t' + chrom + '\t' + str(exons[0][0]-offset) + '\t50\t' + cigar + '\t*\t0\t0\t*\t*\tXS:A:' + read.strand + '\tNH:i:' + str(int(read.NH)) + '\n')
            else:
                filehandle.write(chrom+':'+str(readId) + '\t0\t' + chrom + '\t' + str(exons[0][0]-offset) + '\t50\t' + cigar + '\t*\t0\t0\t*\t*\tNH:i:' + str(int(read.NH)) + '\n')
            readId += 1
        
        for pair in paired:
            exonsA = pair.exonsA
            cigarA = [str(exonsA[0][1] - exonsA[0][0]) + 'M']
            spliced = False
            for i in range(1, len(exonsA)):
                if exonsA[i][0] == exonsA[i-1][1]:
                    prevLen = int(cigarA[-1][:-1])
                    cigarA[-1] = str(prevLen + exonsA[i][1] - exonsA[i][0]) + 'M'
                else:
                    spliced = True
                    cigarA += [str(exonsA[i][0] - exonsA[i-1][1]) + 'N']
                    cigarA += [str(exonsA[i][1] - exonsA[i][0]) + 'M']
            cigarA = ''.join(cigarA)

            exonsB = pair.exonsB
            cigarB = [str(exonsB[0][1] - exonsB[0][0]) + 'M']
            for i in range(1, len(exonsB)):
                if exonsB[i][0] == exonsB[i-1][1]:
                    prevLen = int(cigarB[-1][:-1])
                    cigarB[-1] = str(prevLen + exonsB[i][1] - exonsB[i][0]) + 'M'
                else:
                    spliced = True
                    cigarB += [str(exonsB[i][0] - exonsB[i-1][1]) + 'N']
                    cigarB += [str(exonsB[i][1] - exonsB[i][0]) + 'M']
            cigarB = ''.join(cigarB)

            # Distance from start of first read to end of second read
            totalLen = max(exonsA[-1][1],exonsB[-1][1]) - exonsA[0][0]

            chromA = pair.chromA
            chromB = pair.chromB
            offsetA = self.chromOffsets[chromA]
            if not chromA == chromB:
                offsetB = self.chromOffsets[chromB]

            if force_xs and spliced and not pair.strand:
                #print('Assigning random XS value to spliced paired read')
                if random.randint(0,1) == 0:
                    pair.strand = '+'
                else:
                    pair.strand = '-'

            #if chromA == chromB and self.conflicts(exonsA, exonsB):
            #    if pair.strand:
            #        filehandle.write(chromA+':'+str(readId) + '\t0\t' + chromA + '\t' + str(exonsA[0][0]-offsetA) + '\t50\t' + cigarA + '\t*\t0\t0\t*\t*\tNH:i:' + str(pair.NH) + '\tXS:A:' + pair.strand + '\n')
            #        readId += 1
            #        filehandle.write(chromB+':'+str(readId) + '\t0\t' + chromB + '\t' + str(exonsB[0][0]-offsetA) + '\t50\t' + cigarB + '\t*\t0\t0\t*\t*\tNH:i:' + str(pair.NH) + '\tXS:A:' + pair.strand + '\n')
            #    else:
            #        filehandle.write(chromA+':'+str(readId) + '\t0\t' + chromA + '\t' + str(exonsA[0][0]-offsetA) + '\t50\t' + cigarA + '\t*\t0\t0\t*\t*\tNH:i:' + str(pair.NH) + '\n')
            #        readId += 1
            #        filehandle.write(chromB+':'+str(readId) + '\t0\t' + chromB + '\t' + str(exonsB[0][0]-offsetA) + '\t50\t' + cigarB + '\t*\t0\t0\t*\t*\tNH:i:' + str(pair.NH) + '\n')
            #else:
            if chromB == chromA:
                filehandle.write(chromA+':'+str(readId) + '\t81\t' + chromA + '\t' + str(exonsA[0][0]-offsetA) + '\t50\t' + cigarA + '\t=\t' + str(exonsB[0][0]-offsetA) + '\t' + str(totalLen) + '\t*\t*\tNH:i:' + str(int(pair.NH)))
            else:
                filehandle.write(chromA+':'+str(readId) + '\t81\t' + chromA + '\t' + str(exonsA[0][0]-offsetA) + '\t50\t' + cigarA + '\t' + chromB + '\t' + str(exonsB[0][0]-offsetB) + '\t0\t*\t*\tNH:i:' + str(int(pair.NH)))

            if pair.strand:# and 'N' in cigarA:
                filehandle.write('\tXS:A:' + pair.strand)
            filehandle.write('\n')

            if chromB == chromA:
                filehandle.write(chromA+':'+str(readId) + '\t161\t' + chromB + '\t' + str(exonsB[0][0]-offsetA) + '\t50\t' + cigarB + '\t=\t' + str(exonsA[0][0]-offsetA) + '\t' + str(-totalLen) + '\t*\t*\tNH:i:' + str(int(pair.NH)))
            else:
                filehandle.write(chromA+':'+str(readId) + '\t161\t' + chromB + '\t' + str(exonsB[0][0]-offsetB) + '\t50\t' + cigarB + '\t' + chromA + '\t' + str(exonsA[0][0]-offsetA) + '\t0\t*\t*\tNH:i:' + str(int(pair.NH)))
            if pair.strand:# and 'N' in cigarB:
                filehandle.write('\tXS:A:' + pair.strand)
            filehandle.write('\n')

            readId += 1

        return readId
