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
import preprocess

class Compressor:
    aligned = None

    # 0 - zlib
    # 1 - lzma
    # 2 - bz2
    compressMethod = 0

    def __init__(self, frag_len_cutoff):
        if self.compressMethod == 0:
            self.zlib = __import__('zlib')
        elif self.compressMethod == 1:
            self.lzma = __import__('lzma')
        elif self.compressMethod == 2:
            self.bz2 = __import__('bz2')

        if frag_len_cutoff:
            print('Set fragment length cutoff to %d' % frag_len_cutoff)
        self.frag_len_cutoff = frag_len_cutoff

    def compress(self, samFilename, compressedFilename, min_filename, frag_len_z_cutoff, split_diff_strands, split_discordant):
        ''' Compresses the alignments to 2 files, one for unspliced and one for spliced

            file_prefix: Prefix for all output file names
        '''

        self.p = preprocess.Preprocessor(samFilename, frag_len_z_cutoff, split_diff_strands)

        if not self.frag_len_cutoff:
            self.frag_len_cutoff = self.p.frag_len_cutoff
        print('Using fragment length cutoff of ' + str(self.frag_len_cutoff))

        if split_diff_strands:
            print('Splitting mates on different strands')
        else:
            print('Not splitting mates on different strands')

        if split_discordant:
            print('Splitting discordant')
        else:
            print('Not splitting discordant')

        # Reads on different strands that should be unpaired
        self.diff_strand_unpaired = self.p.unpaired
        del self.p

        # Read header
        header = ''
        with open(samFilename, 'r') as f:
            for line in f:
                if line[0] == '@':
                    header += line
                else:
                    break
        self.chromosomes = self.parseSAMHeader(header)
        self.aligned = alignments.Alignments(self.chromosomes, self.frag_len_cutoff, split_discordant)

        self.compressByBundle(samFilename, compressedFilename, min_filename)

    def compressByBundle(self, input_name, compressed_name, intermediate_name=None):
        '''
        Read a sorted SAM file and compress in segments determined by clusters of reads

        :param filename:
        :return:
        '''

        # If coverage is 0 for at least this many bases end of a potential gene
        overlapRadius = 50

        spliced_index = []
        bundles = []

        first = True

        bundle_id = 0

        read_id = 0

        diff_strand_unpaired_id = 0
        num_diff_strand_unpaired = len(self.diff_strand_unpaired)

        with open(input_name, 'r') as filehandle:
            id = 0
            for line in filehandle:
                # Check if header line
                if line[0] == '@':
                    continue

                row = line.strip().split('\t')

                if not row[2] in self.chromosomes.keys():
                    print('Error! Chromosome ' + str(row[2]) + ' not found!')
                    exit()

                # Starting position of this read
                start = self.aligned.chromOffsets[row[2]] + int(row[3])

                if self.aligned.gene_bounds and start > (self.aligned.gene_bounds[-1] + overlapRadius):
                    # Compress most recent bundle
                    self.aligned.finalizeExons()
                    self.aligned.finalizeUnmatched()
                    self.aligned.finalize_cross_bundle_reads()
                    bundle_id += 1

                    bundles.append(self.aligned.exons)

                    # Write to intermediate file
                    if intermediate_name:
                        if first:
                            # If it's the first bundle, write the header as well
                            with open(intermediate_name, 'w') as f1:
                                read_id = self.aligned.writeSAM(f1, self.aligned.unpaired, self.aligned.paired, True, False, read_id)
                        else:
                            with open(intermediate_name, 'a') as f1:
                                read_id = self.aligned.writeSAM(f1, self.aligned.unpaired, self.aligned.paired, False, False, read_id)

                    junctions, maxReadLen = self.aligned.computeBuckets()
                    self.sortedJuncs = sorted(junctions.keys())

                    # Compress bundle to temporary file
                    if first:
                        mode = 'wb'
                    else:
                        mode = 'ab'
                    with open('temp.bin', mode) as f:
                        l = self.compressBundle(junctions, maxReadLen, f)
                        spliced_index.append(l)

                    # Start new bundle
                    self.aligned.resetBundle()
                    self.aligned.exons.add(start)

                    first = False

                # Process read
                exons = self.parseCigar(row[5], int(row[3]))

                # find XS (strand) and NH values
                strand = None
                NH = 1
                for r in row[11 : len(row)]:
                    if r[0:5] == 'XS:A:' or r[0:5] == 'XS:a:':
                        strand = r[5]
                    elif r[0:3] == 'NH:':
                        NH = int(r[5:])

                r = read.Read(row[2], exons, strand, NH)
                #r.name = row[0]

                if row[6] == '*':
                    paired = False
                elif diff_strand_unpaired_id < num_diff_strand_unpaired and id == self.diff_strand_unpaired[diff_strand_unpaired_id]:
                    #if not row[6] == '*':
                    #    print('\t'.join(row))
                    paired = False
                    diff_strand_unpaired_id += 1
                else:
                    paired = True
                    r.bundle = bundle_id
                    r.pairOffset = int(row[7])
                    if row[6] == '=':
                        r.pairChrom = row[2]
                    else:
                        r.pairChrom = row[6]
                self.aligned.processRead(row[0], r, paired)

                id += 1

            # Compress final cluster
            self.aligned.finalizeExons()
            self.aligned.finalizeUnmatched()
            self.aligned.finalize_cross_bundle_reads()
            bundle_id += 1

            bundles.append(self.aligned.exons)

            # Write to intermediate file
            if intermediate_name:
                if first:
                    # If it's the first bundle, write the header as well
                    with open(intermediate_name, 'w') as f1:
                        read_id = self.aligned.writeSAM(f1, self.aligned.unpaired, self.aligned.paired, True, False, read_id)
                    first = False
                else:
                    with open(intermediate_name, 'a') as f1:
                        read_id = self.aligned.writeSAM(f1, self.aligned.unpaired, self.aligned.paired, False, False, read_id)

            junctions, maxReadLen = self.aligned.computeBuckets()
            self.sortedJuncs = sorted(junctions.keys())

            # Compress bundle to temporary file
            if first:
                mode = 'wb'
            else:
                mode = 'ab'
            with open('temp.bin', mode) as f:
                l = self.compressBundle(junctions, maxReadLen, f)
                spliced_index.append(l)

        leftovers = 0
        for k,v in self.aligned.cross_bundle_reads.items():
            #if len(v) > 0:
            #    print(k)
            #    print(v)
            #    exit()
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

            # Compress bundle-spanning buckets
            self.compressCrossBundle(self.aligned.cross_bundle_buckets, self.aligned.max_cross_bundle_read_len, bundle_id, f)

            # Move contents of temporary file to output file
            with open('temp.bin', 'rb') as f2:
                f.write(f2.read())

        os.remove('temp.bin')


    def compressBundle(self, junctions, maxReadLen, filehandle):
        # Determine the number of bytes for read lengths
        readLenBytes = binaryIO.findNumBytes(maxReadLen)
        cluster = binaryIO.valToBinary(1, readLenBytes)
        cluster += binaryIO.writeJunctionsList(self.sortedJuncs, 2)

        # TODO: No need for junc_lens?
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

            chunk_size = 20
            num_chunks = math.ceil(len(buckets_sorted) / chunk_size)
            chunk_lens = [0] * num_chunks

            index = binaryIO.valToBinary(4, len(buckets_sorted))
            index += binaryIO.valToBinary(2, chunk_size)
            index += binaryIO.valToBinary(1, readLenBytes)
            index += binaryIO.writeCrossBundleBucketNames(bundleIdBytes, cross_bundle_buckets, buckets_sorted)

            main = b''
            chunk = b''
            chunk_id = 0
            for i in range(len(buckets_sorted)):
                b = buckets_sorted[i]

                chunk += binaryIO.writeCrossBundleBucket(readLenBytes, cross_bundle_buckets[b])
                if (i+1) % chunk_size == 0:
                    compressed = self.compressString(chunk)
                    chunk_lens[chunk_id] = len(compressed)
                    chunk_id += 1
                    main += compressed
                    chunk = b''

            if len(chunk) > 0:
                compressed = self.compressString(chunk)
                chunk_lens[chunk_id] = len(compressed)
                main += compressed

            index += binaryIO.writeList(chunk_lens)

            index = self.compressString(index)
            length = len(index)
            numBytes = binaryIO.findNumBytes(length)
            binaryIO.writeVal(filehandle, 1, numBytes)
            binaryIO.writeVal(filehandle, numBytes, length)
            filehandle.write(index)
            filehandle.write(main)

            print('Compressed size: %d' % (filehandle.tell() - pos))
        else:
            binaryIO.writeVal(filehandle, 1, 1)
            binaryIO.writeVal(filehandle, 1, 0)

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
