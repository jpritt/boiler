# Contains functions for the first pass over the alignments file to determine the fragment length distribution, establish a cutoff for long fragments, and determine pairing
import math
import re

class Preprocessor:
    def __init__(self, samFilename, cutoff_z, split_diff_strand, split_discordant):
        self.num_reads = sum(1 for line in open(samFilename))
        self.frag_len_cutoff = None

        self.preprocess(samFilename, cutoff_z, split_diff_strand, split_discordant)

    def preprocess(self, samFilename, cutoff_z, split_diff_strand, split_discordant):
        '''
        Make a first pass over the file
        '''

        # For counting fragment lengths
        find_cutoff = False
        if not cutoff_z == None:
            find_cutoff = True
            self.lens = dict()
            self.len_sum = 0

        # For determining pairs
        self.pairing = [-1] * self.num_reads
        self.unmatched = dict()

        with open(samFilename, 'r') as f:
            id = 0
            for line in f:
                row = line.rstrip().split('\t')
                if len(row) < 6:
                    continue

                if find_cutoff and row[6] == '=':
                    self.process_frag_lens(int(row[7]) - int(row[3]))

                if not row[6] == '*':
                    self.process_pairs(id, row)

                id += 1

        if find_cutoff:
            self.calculate_cutoff(cutoff_z)
        self.reprocess_pairs(split_diff_strand, split_discordant)

    def process_pairs(self, id, row):
        name = row[0]
        chrom = row[2]
        pos = int(row[3])
        exons = self.parseCigar(row[5], pos)
        if row[6] == '=':
            mate_chrom = chrom
        else:
            mate_chrom = row[6]
        mate_pos = int(row[7])

        # -1 = negative strand, 1 = positive, 0 = unknown
        strand = 0
        for r in row[11:]:
            if r[:5] == 'XS:A:':
                if r[5] == '+':
                    strand = 1
                elif r[5] == '-':
                    strand = -1

        if name in self.unmatched:
            foundMatch = False

            mates = self.unmatched[name]
            for i in range(len(mates)):
                mate = mates[i]
                if chrom == mate[3] and mate_chrom == mate[1] and pos == mate[4] and mate_pos == mate[2] and (not chrom == mate_chrom or not self.conflicts(exons, mate[5])) and (strand == 0 or mate[6] == 0 or strand == mate[6]):
                    # Perfect match
                    self.pairing[id] = mate[0]
                    self.pairing[mate[0]] = id

                    # TODO: Delete list if empty?
                    del self.unmatched[name][i]
                    foundMatch = True
                    break
            if not foundMatch:
                self.unmatched[name].append((id, chrom, pos, mate_chrom, mate_pos, exons, strand))
        else:
            self.unmatched[name] = [(id, chrom, pos, mate_chrom, mate_pos, exons, strand)]

    def reprocess_pairs(self, split_diff_strand, split_discordant):
        '''
        Second pass looking for pairs, allow different strands and discordant if not required to split them
        :param split_diff_strand:
        :param split_discordant:
        :return:
        '''

        count = 0
        for name,reads in self.unmatched.items():
            count += len(reads)
        print('%d reads remaining before second pass' % count)

        count = 0
        for name,reads in self.unmatched.items():
            if len(reads) > 1:
                i = 0
                while i < len(reads)-1:
                    foundMatch = False
                    j = i+1
                    while j < len(reads):
                        if reads[i][1] == reads[j][3] and reads[i][3] == reads[j][1] and reads[i][2] == reads[j][4] and reads[i][4] == reads[j][2]:
                            if not split_discordant or not reads[i][1] == reads[i][3] or not self.conflicts(reads[i][5], reads[j][5]):
                                if not split_diff_strand or (reads[i][6] == 0 or reads[j][6] == 0 or reads[i][6] == reads[j][6]):
                                    self.pairing[reads[i][0]] = reads[j][0]
                                    self.pairing[reads[j][0]] = reads[i][0]

                                    count += 2

                                    del self.unmatched[name][j]
                                    del self.unmatched[name][i]

                                    foundMatch = True

                        if foundMatch:
                            break

                        j += 1
                    if not foundMatch:
                        i += 1

        print('%d reads matched up' % count)

        # All remaining reads are unmatched
        count = 0
        for name,reads in self.unmatched.items():
            count += len(reads)
        print('%d reads remaining' % count)

    def get_pair(self, i):
        return self.pairing[i]

    def process_frag_lens(self, frag_len):
        if frag_len >= 0:
            self.len_sum += frag_len

        if frag_len in self.lens:
            self.lens[frag_len] += 1
        else:
            self.lens[frag_len] = 1

    def calculate_cutoff(self, cutoff_z):
        if self.num_reads == 0:
            self.frag_len_cutoff = 0
            return

        # Calculate average and standard deviation
        avg = float(self.len_sum) / float(self.num_reads)
        stdev = 0.0
        for length,freq in self.lens.items():
            stdev += ((length - avg) ** 2) * freq
        stdev = math.sqrt(stdev / self.num_reads)
        self.frag_len_cutoff = int(avg + cutoff_z * stdev)

        print('Set fragment length cutoff to z=%f (%d) based on length distribution' % (cutoff_z, self.frag_len_cutoff))
        count_longer = 0
        for l,f in self.lens.items():
            if l > self.frag_len_cutoff:
                count_longer += f
        print('%0.2f %% of pairs are longer than the cutoff' % (100.0 * float(count_longer) / float(self.num_reads)))

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
