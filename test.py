#! /usr/bin/env python

import alignments
import time
import math


def findReadsInCoverage_v3(cov, frag_lens):
    reads = []

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
        reads, start, end = findReadLeft(cov, start, end, reads, min_len, max_len, mode_len)

        #if start < end:
        #    start, end = self.findReadRight(cov, start, end, reads, min_len, max_len, mode_len)

    reads = finalizeCovReads(reads, min_len, max_len, mode_len)

    return reads

def finalizeCovReads(reads, min_len, max_len, mode_len):
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

    print('')
    print('%d, %d, %d' % (min_len, max_len, mode_len))
    print(final)
    return final

def findReadsEndingAt(reads, end):
    r = []
    for i in range(len(reads)-1, -1, -1):
        if reads[i][1] == end:
            r.append(i)
    return r

def findReadsStartingAt(reads, start):
    r = []
    for i in range(len(reads)-1, -1, -1):
        if reads[i][0] == start:
            r.append(i)
    return r

def findNextGap(cov, start, end, modeLen):
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

def extendReadRight(cov, start, stop, reads, max_len):
    '''
    Extend the right end of a read ending at 'start' up to 'stop'
    '''
    for i in range(len(reads)):
        r = reads[i]
        if r[1] < start:
            return False
        elif r[1] == start:
            if (stop - r[0]) <= (max_len * r[2]):
                # Extend existing reads
                r[1] = stop
                reads = preserve_sorted(reads, i)
                for i in range(start, stop):
                    cov[i] -= 1
                return reads
            elif (stop - r[0]) >= (min_len * (r[2]+1)):
                # Add new read
                r[1] = stop
                r[2] += 1
                reads = preserve_sorted(reads, i)
                for i in range(start, stop):
                    cov[i] -= 1
                return reads
    return False

def shortenReadRight(cov, start, stop, reads, max_shorten, min_len):
    '''
    Shorten the right end of a read ending at 'stop' down to 'start'
    '''

    for i in range(len(reads)):
        r = reads[i]
        if r[1] < stop:
            return False
        elif (start - r[0]) >= (min_len * r[2]):
            r[1] = start
            reads = preserve_sorted(reads, i)
            for j in range(start, i):
                cov[j] += 1
            return reads
    return False

def addToMiddle(cov, start, stop, reads, min_len, max_len):
    for i in range(len(reads)):
        r = reads[i]
        if r[1] < start:
            return False
        elif stop-r[0] >= min_len and r[1]-start >= min_len:
            new_read = [r[0], stop, r[2]-1]

            if (new_read[1] - new_read[0]) > new_read[2]*max_len:
                new_read[2] = math.ceil((new_read[1] - new_read[0]) / max_len)
            elif (new_read[1] - new_read[0]) < new_read[2]*min_len:
                new_read[2] = math.floor((new_read[1] - new_read[0]) / min_len)
            r[0] = start
            r[2] = 1
            reads = add_read(reads, new_read)
            for j in range(start, stop):
                cov[j] -= 1
            return reads
    return False

def add_read(reads, new_read):
    # Try to find a string of reads ending where this one starts
    for i in range(len(reads)):
        r = reads[i]
        if r[1] == new_read[0]:
            r[1] = new_read[1]
            r[2] += 1

            reads = preserve_sorted(reads, i)

            return reads
    # If no reads end where this one starts, just add it to the list of reads
    return insert_in_order(reads, [new_read[0], new_read[1], 1])

def preserve_sorted(reads, updated_i):
    '''

    :param reads: List of reads
    :param updated_i: Last read updated

    Reorder reads so they are sorted by read end, in decreasing order
    '''

    new_i = updated_i
    while new_i > 0 and reads[new_i][1] > reads[new_i-1][1]:
        new_i -= 1

    num_reads = len(reads)
    while new_i < num_reads-1 and reads[new_i][1] < reads[new_i+1][1]:
        new_i += 1

    return reads[:new_i] + [reads[updated_i]] + reads[new_i:updated_i] + reads[updated_i+1:]

def insert_in_order(reads, r):
    '''

    :param reads:
    :param r:
    Insert r into reads such that the reads remain sorted by read end, in decreasing order
    '''

    num_reads = len(reads)
    i = 0
    while i < (num_reads-1) and r[1] < reads[i][1]:
        i += 1

    return reads[:i] + [r] + reads[i:]

def findReadLeft(cov, start, end, reads, min_len, max_len, mode_len):
    gap = findNextGap(cov, start, end, mode_len)

    len_range = max_len - min_len

    if not gap:
        readEnd = min(start+mode_len,end)

        r = None
        if readEnd-start < min_len:
            r = extendReadRight(cov, start, readEnd, reads, max_len)
            if r:
                reads = r
            else:
                r = addToMiddle(cov, start, readEnd, reads, min_len, max_len)
                if r:
                    reads = r

                else:
                    print('Skipping final read of length %d' % (readEnd-start))
                    for i in range(start, readEnd):
                        cov[i] -= 1
        else:
            reads = add_read(reads, [start, readEnd])
            for i in range(start, readEnd):
                cov[i] -= 1
    else:
        r = extendReadRight(cov, start, start+gap[0], reads, max_len)
        if r:
            reads = r

        if not r:
            r = addToMiddle(cov, start, start+gap[0], reads, min_len, max_len)
            if r:
                reads = r

        if not r:
            r = shortenReadRight(cov, start+gap[0], start+gap[1], reads, len_range, min_len)
            if r:
                reads = r

        if not r:
            if gap[0] >= min_len:
                readEnd = start+gap[0]
            else:
                readEnd = min(start+mode_len, end)

            reads = add_read(reads, [start, readEnd])
            for i in range(start, readEnd):
                cov[i] -= 1

    while start < end and cov[start] <= 0:
        start += 1
    while start < end and cov[end-1] <= 0:
        end -= 1

    #print(reads)

    return reads, start, end


def RLEtoVector(rle):
    ''' Convert a run length encoded vector to the original coverage vector '''

    vector = []
    for row in rle:
        vector += [row[0]] * row[1]
    return vector

#covRLE = [[3, 36], [1, 40], [0, 128], [1, 76], [0, 3]]
#covRLE = [[1, 2], [3, 2], [6, 3], [7, 1], [10, 1], [11, 3], [12, 2], [14, 4], [16, 1], [17, 1], [18, 3], [19, 2], [21, 3], [22, 1], [23, 1], [26, 1], [28, 3], [29, 1], [32, 2], [33, 1], [34, 2], [36, 1], [38, 1], [41, 1], [43, 3], [44, 1], [46, 4], [47, 1], [48, 1], [50, 4], [52, 1], [54, 1], [55, 1], [56, 2], [61, 2], [62, 1], [64, 1], [66, 3], [68, 2], [70, 2], [71, 1], [72, 1], [74, 1], [75, 1], [77, 1], [76, 1], [79, 1], [76, 1], [77, 1], [79, 1], [78, 1], [80, 2], [82, 1], [84, 1], [87, 2], [86, 1], [87, 1], [90, 1], [91, 1], [90, 1], [91, 1], [90, 1], [91, 1], [92, 1], [91, 1], [94, 1], [92, 2], [93, 1], [94, 2], [92, 1], [91, 1], [93, 1], [97, 1], [99, 2], [100, 1], [104, 3], [102, 1], [105, 1], [103, 2], [106, 1], [108, 2], [107, 1], [108, 1], [110, 1], [114, 1], [115, 1], [116, 1], [115, 1], [116, 1], [119, 3], [117, 1], [118, 1], [117, 2], [112, 2], [109, 1], [107, 1], [105, 2], [107, 1], [105, 2], [103, 2], [102, 1], [101, 1], [99, 1], [97, 1], [95, 1], [96, 1], [93, 2], [92, 1], [90, 2], [85, 1], [84, 1], [82, 1], [80, 1], [76, 1], [78, 1], [77, 1], [76, 1], [73, 1], [72, 1], [71, 1], [69, 2], [68, 1], [67, 2], [62, 3], [61, 1], [59, 1], [58, 1], [57, 1], [56, 1], [54, 1], [50, 1], [47, 1], [42, 1], [41, 1], [36, 1], [35, 3], [32, 1], [31, 1], [29, 1], [26, 1], [24, 1], [23, 1], [22, 1], [19, 1], [17, 1], [13, 1], [11, 1], [9, 1], [8, 1], [7, 1], [2, 2]]
covRLE = [[0, 3], [1, 50], [2, 12], [3, 14], [2, 50], [1, 12], [0, 1782731], [1, 4], [2, 4], [3, 1], [4, 3], [5, 6], [7, 4], [8, 3], [9, 1], [11, 4], [12, 1], [13, 3], [14, 1], [15, 2], [16, 2], [17, 6], [19, 1], [20, 1], [21, 1], [22, 3], [23, 1], [24, 1], [26, 1], [29, 1], [30, 1], [32, 1], [33, 2], [39, 1], [43, 2], [46, 1], [49, 1], [54, 1], [57, 1], [59, 1], [66, 1], [68, 1], [74, 1], [75, 1], [80, 1], [83, 1], [87, 1], [89, 2], [90, 1], [93, 1], [97, 1], [99, 1], [104, 1], [107, 1], [109, 1], [115, 1], [121, 1], [124, 1], [127, 1], [131, 1], [138, 1], [141, 1], [147, 1], [152, 1], [163, 1], [169, 1], [170, 1], [178, 1], [184, 1], [193, 1], [196, 1], [205, 1], [212, 1], [219, 1], [228, 1], [243, 1], [251, 1], [257, 1], [261, 1], [266, 1], [270, 1], [280, 1], [284, 1], [290, 1], [297, 1], [302, 1], [305, 1], [312, 1], [319, 1], [326, 1], [328, 1], [333, 1], [340, 1], [344, 1], [345, 1], [351, 1], [355, 1], [360, 1], [365, 1], [367, 1], [371, 1], [759, 1], [765, 1], [768, 1], [772, 1], [776, 1], [780, 1], [783, 1], [789, 1], [795, 1], [796, 1], [801, 1], [800, 1], [802, 1], [806, 1], [805, 1], [811, 1], [812, 1], [821, 1], [818, 1], [819, 2], [817, 1], [819, 1], [823, 2], [822, 1], [820, 1], [818, 1], [815, 2], [811, 1], [808, 1], [806, 1], [805, 2], [799, 1], [796, 1], [790, 1], [789, 1], [780, 1], [776, 1], [773, 1], [765, 1], [763, 1], [754, 1], [752, 1], [743, 1], [744, 1], [740, 1], [729, 1], [714, 1], [706, 1], [702, 1], [701, 1], [699, 1], [697, 1], [687, 1], [682, 1], [675, 1], [668, 1], [662, 1], [659, 1], [651, 1], [644, 1], [637, 1], [635, 1], [630, 1], [623, 1], [617, 1], [615, 1], [608, 1], [603, 1], [598, 1], [593, 1], [590, 1], [585, 1], [0, 225733]]
cov = RLEtoVector(covRLE)
#readLens = {76: 2, 36: 2}
#readLens = {32: 2, 52: 2, 23: 2, 42: 2, 76: 169, 46: 2}
readLens = {21: 2, 22: 4, 23: 4, 24: 2, 28: 4, 29: 8, 31: 2, 33: 4, 36: 2, 37: 2, 38: 4, 41: 2, 42: 4, 43: 2, 44: 2, 45: 4, 46: 2, 47: 2, 49: 4, 51: 3, 52: 3, 53: 6, 54: 2, 56: 4, 57: 4, 58: 2, 59: 10, 60: 7, 61: 8, 62: 6, 63: 6, 64: 5, 65: 4, 66: 8, 67: 4, 68: 6, 69: 10, 70: 9, 71: 4, 72: 5, 73: 6, 74: 4, 75: 9, 76: 788}


chroms = {'chr1': 0, 'chr7': 2434317170, 'chr9': 2739819855, 'chrY': 3036320417, 'chr16': 978694589, 'chr19': 1228321800, 'chr8': 2593455833, 'chr13': 653643779, 'chr17': 1069049342, 'chr12': 519791884, 'chr3': 1693110137, 'chr20': 1530650156, 'chr14': 768813657, 'chr5': 2082286843, 'chr18': 1150244552, 'chr11': 384785368, 'chrX': 2881049857, 'chr4': 1891132567, 'chr10': 249250621, 'chrM': 2881033286, 'chr6': 2263202103, 'chr22': 1641805571, 'chr21': 1593675676, 'chr15': 876163197, 'chr2': 1287450783}
#aligned = alignments.Alignments(chroms, False)

reads = findReadsInCoverage_v3(cov, readLens)

#print(reads)
#print('')

lens = dict()
for r in reads:
    l = r[1] - r[0]
    if l in lens:
        lens[l] += 1
    else:
        lens[l] = 1
print(sorted(readLens.items()))
print(sorted(lens.items()))