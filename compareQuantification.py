#! /usr/bin/env python
import sys
from scipy.stats.stats import pearsonr

'''
Compare 2 .fpkm_tracking files output by Cufflinks (run with the -G option)
'''

threshold = 0.001

names = []
fpkmsOrig = []
with open(sys.argv[1], 'r') as f:
    header = True
    label_id = 0
    fpkm_id = 0
    for line in f:
        row = line.rstrip().split('\t')

        if header:
            label_id = row.index('tracking_id')
            fpkm_id = row.index('FPKM')
            header = False
        else:
            fpkm = float(row[fpkm_id])
            names.append(row[label_id])
            fpkmsOrig.append(float(row[fpkm_id]))

fpkmsA = []
fpkmsB = []
with open(sys.argv[1], 'r') as f:
    header = True
    label_id = 0
    fpkm_id = 0
    line_id = 0
    for line in f:
        row = line.rstrip().split('\t')

        if header:
            label_id = row.index('tracking_id')
            fpkm_id = row.index('FPKM')
            header = False
        else:
            name = row[label_id]
            if not name == names[line_id]:
                print('Error! Files do not contain the same isoforms in the same order!')
            fpkm = float(row[fpkm_id])
            if fpkm > threshold or fpkmsOrig[line_id] > threshold:
                fpkmsA.append(fpkmsOrig[line_id])
                fpkmsB.append(fpkm)
            line_id += 1

print(fpkmsA[:20])
print(fpkmsB[:20])
print('Comparing %d transcripts' % len(fpkmsA))
p = pearsonr(fpkmsA, fpkmsB)
print('Pearson correlation coefficient:\t%f, %f' % (p[0], p[1]))
