#! /usr/bin/env python3
import alignments
import re
import read
import pairedread
import bucket
import bisect
import binaryIO
import math
import readNode
import time
import os

class Compressor:
    aligned = None

    # 0 - zlib
    # 1 - lzma
    # 2 - bz2
    compressMethod = 0

    def __init__(self, force_xs, frag_len_cutoff):
        if self.compressMethod == 0:
            self.zlib = __import__('zlib')
        elif self.compressMethod == 1:
            self.lzma = __import__('lzma')
        elif self.compressMethod == 2:
            self.bz2 = __import__('bz2')

        self.force_xs = force_xs

        if frag_len_cutoff:
            print('Set fragment length cutoff to %d' % frag_len_cutoff)
        self.frag_len_cutoff = frag_len_cutoff

    def compress(self, samFilename, compressedFilename, min_filename=None, frag_len_z_cutoff=None, binary=False, debug=False):
        ''' Compresses the alignments to 2 files, one for unspliced and one for spliced

            file_prefix: Prefix for all output file names
        '''

        if not(frag_len_z_cutoff == None) and not self.frag_len_cutoff:
            self.firstPass(samFilename, frag_len_z_cutoff)
        elif not self.frag_len_cutoff:
            print('No fragment length cutoff')

        self.debug = debug

        # Read header
        header = ''
        with open(samFilename, 'r') as f:
            for line in f:
                if line[0] == '@':
                    header += line
                else:
                    break
        self.chromosomes = self.parseSAMHeader(header)
        self.aligned = alignments.Alignments(self.chromosomes, self.frag_len_cutoff, self.debug)

        self.compressByCluster(samFilename, compressedFilename, min_filename)

    def firstPass(self, samFilename, cutoff_z):
        '''
        Make a first pass over the file to determine the fragment length distribution and establish a cutoff for long fragments
        '''


        lens = dict()
        len_sum = 0
        N = 0
        with open(samFilename, 'r') as f:
            for line in f:
                row = line.rstrip().split('\t')
                if len(row) < 6:
                    continue

                if row[6] == '=':
                    frag_len = int(row[7]) - int(row[3])
                    if frag_len >= 0:
                        len_sum += frag_len
                        N += 1

                        if frag_len in lens:
                            lens[frag_len] += 1
                        else:
                            lens[frag_len] = 1

        # Calculate average and standard deviation
        avg = float(len_sum) / float(N)
        stdev = 0.0
        for length,freq in lens.items():
            stdev += ((length - avg) ** 2) * freq
        stdev = math.sqrt(stdev / N)

        #print('%d fragments' % N)
        #print('Average fragment length: %f' % avg)
        #print('Standard deviation: %f' % stdev)
        #print('')

        self.frag_len_cutoff = int(avg + cutoff_z * stdev)

        print('Set fragment length cutoff to z=%d (%d) based on length distribution' % (cutoff_z, self.frag_len_cutoff))
        count_longer = 0
        for l,f in lens.items():
            if l > self.frag_len_cutoff:
                count_longer += f
        print('%0.2f %% of pairs are longer than the cutoff' % (100.0 * float(count_longer) / float(N)))

        # Test different z scores
        #zs = [2, 2.5, 3, 3.5, 4, 5, 6, 7, 8]
        #cutoffs = [avg + c * stdev for c in zs]
        #num_c = len(cutoffs)
        #counts = [0] *  num_c

        #for length,freq in lens.items():
        #    for i in range(num_c):
        #        if length > cutoffs[i]:
        #            counts[i] += freq
        #        else:
        #            break

        #for i in range(num_c):
        #    print('z-score %0.1f (%d): %d / %d = %0.3f %%' % (zs[i], cutoffs[i], counts[i], N, 100.0*float(counts[i])/float(N)))
        #print('')


    def compressCluster(self, junctions, maxReadLen, filehandle):
        # Determine the number of bytes for read lengths
        readLenBytes = binaryIO.findNumBytes(maxReadLen)
        cluster = binaryIO.valToBinary(1, readLenBytes)
        cluster += binaryIO.writeJunctionsList(self.sortedJuncs, 2)

        junc_lens = []
        junc_string = b''
        for j in self.sortedJuncs:
            s = binaryIO.writeJunction(readLenBytes, junctions[j])
            junc_lens.append(len(s))
            junc_string += s

        #cluster += binaryIO.writeList(junc_lens)
        cluster += junc_string

        # Write to file
        start = filehandle.tell()
        filehandle.write(self.compressString(cluster))

        # return length of cluster in file
        return filehandle.tell() - start

    def compressCrossBundle(self, cross_bundle_buckets, maxReadLen, num_bundles, filehandle):
        '''
        Compress the bundle-spanning buckets
        '''

        readLenBytes = binaryIO.findNumBytes(maxReadLen)

        bundleIdBytes = binaryIO.findNumBytes(num_bundles)

        buckets_sorted = sorted(cross_bundle_buckets.keys())

        if len(buckets_sorted) > 0:
            print('%d cross-bundle buckets' % len(buckets_sorted))
            pos = filehandle.tell()

            s = b''
            #bucket_lens = []
            for b in buckets_sorted:
                start_len = len(s)
                s += binaryIO.writeCrossBundleBucket(bundleIdBytes, readLenBytes, cross_bundle_buckets[b])
                #bucket_lens.append(len(s)-start_len)
            #print('Bucket lengths: ' + str(bucket_lens))

            s = self.compressString(s)
            length = len(s)
            numBytes = binaryIO.findNumBytes(length)
            binaryIO.writeVal(filehandle, 1, readLenBytes)
            binaryIO.writeVal(filehandle, 1, numBytes)
            binaryIO.writeVal(filehandle, numBytes, length)
            filehandle.write(s)

            print('Compressed size: %d' % (filehandle.tell() - pos))
        else:
            binaryIO.writeVal(filehandle, 1, readLenBytes)
            binaryIO.writeVal(filehandle, 1, 1)
            binaryIO.writeVal(filehandle, 1, 0)


    def compressByCluster(self, input_name, compressed_name, intermediate_name=None):
        '''
        Read a sorted SAM file and compress in segments determined by clusters of reads

        :param filename:
        :return:
        '''

        # If coverage is 0 for at least this many bases end of a potential gene
        overlapRadius = 50

        if self.debug:
            start_time = time.time()

        spliced_index = []
        bundles = []

        first = True

        bundle_id = 0

        with open(input_name, 'r') as filehandle:
            for line in filehandle:
                row = line.strip().split('\t')

                # Check if header line
                if len(row) < 6 or row[5] == '*':
                    continue

                if not row[2] in self.chromosomes.keys():
                    print('Error! Chromosome ' + str(row[2]) + ' not found!')
                    exit()

                start = self.aligned.chromOffsets[row[2]] + int(row[3])

                if self.aligned.gene_bounds and start > (self.aligned.gene_bounds[-1] + overlapRadius):
                    # Compress most recent cluster
                    self.aligned.finalizeUnmatched()
                    self.aligned.finalizeExons()
                    self.aligned.finalize_cross_bundle_reads(bundle_id)
                    bundle_id += 1

                    bundles.append(self.aligned.exons)

                    if intermediate_name:
                        if first:
                            with open(intermediate_name, 'w') as f1:
                                self.aligned.writeSAM(f1, True, self.force_xs)
                            first = False
                        else:
                            with open(intermediate_name, 'a') as f1:
                                self.aligned.writeSAM(f1, False, self.force_xs)

                    junctions, maxReadLen = self.aligned.computeJunctions()
                    self.sortedJuncs = sorted(junctions.keys())

                    with open('temp.bin', 'ab') as f:
                        l = self.compressCluster(junctions, maxReadLen, f)
                        spliced_index.append(l)

                    # Start new cluster
                    self.aligned.resetCluster()
                    self.aligned.exons.add(start)

                exons = self.parseCigar(row[5], int(row[3]))

                # find XS and NH values
                xs = None
                NH = 1
                for r in row[11 : len(row)]:
                    if r[0:5] == 'XS:A:' or r[0:5] == 'XS:a:':
                        xs = r[5]
                    elif r[0:3] == 'NH:':
                        NH = int(r[5:])


                r = read.Read(row[2], exons, xs, NH)
                if row[6] == '*' or row[6] == row[2]:
                    self.aligned.processRead(r, row[0], paired=False)
                else:
                    if row[6] == '=':
                        r.pairChrom = row[2]
                    else:
                        r.pairChrom = row[6]
                    r.pairOffset = int(row[7])
                    self.aligned.processRead(r, row[0], paired=True)


            # Compress final cluster
            self.aligned.finalizeUnmatched()
            self.aligned.finalizeExons()
            self.aligned.finalize_cross_bundle_reads(bundle_id)
            bundle_id += 1

            bundles.append(self.aligned.exons)

            if intermediate_name:
                if first:
                    with open(intermediate_name, 'w') as f1:
                        self.aligned.writeSAM(f1, True, self.force_xs)
                    first = False
                else:
                    with open(intermediate_name, 'a') as f1:
                        self.aligned.writeSAM(f1, False, self.force_xs)

            junctions, maxReadLen = self.aligned.computeJunctions()
            self.sortedJuncs = sorted(junctions.keys())

            with open('temp.bin', 'ab') as f:
                l = self.compressCluster(junctions, maxReadLen, f)
                spliced_index.append(l)

        leftovers = 0
        for k,v in self.aligned.cross_bundle_reads.items():
            leftovers += len(v)
        print('%d cross-bundle reads unmatched' % leftovers)

        bundle_lens = [c[-1]-c[0] for c in bundles]
        print('Minimum bundle length: %d' % min(bundle_lens))
        print('Maximum bundle length: %d' % max(bundle_lens))
        print('Average bundle length: %d'% (sum(bundle_lens) / len(bundle_lens)))

        # Write index information and append spliced and unspliced files
        with open(compressed_name, 'wb') as f:
            s = binaryIO.writeChroms(self.aligned.chromosomes)
            s += binaryIO.writeClusters(bundles)
            s += binaryIO.writeList(spliced_index)
            f.write(s)

            self.compressCrossBundle(self.aligned.cross_bundle_buckets, self.aligned.max_cross_bundle_read_len, bundle_id, f)

            with open('temp.bin', 'rb') as f2:
                f.write(f2.read())

        os.remove('temp.bin')


        if self.debug:
            end_time = time.time()
            print('Compression time: %0.3fs' % (end_time - start_time))

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

    def parseSAMHeader(self, header):
        chromosomes = dict()
        for line in header.split('\n'):
            if line[0:3] == '@SQ':
                row = line.strip().split('\t')
                chromosomes[row[1][3:]] = int(row[2][3:])
        return chromosomes

    def compressString(self, s):
        ''' Use a predefined python library to compress the given string.
            Return the compressed string '''

        if self.compressMethod == 0:
            return self.zlib.compress(s)
        elif self.compressMethod == 1:
            return self.lzma.compress(s)
        elif self.compressMethod == 2:
            return self.bz2.compress(s)
