#! /usr/bin/env python
import sys
from scipy.stats.stats import pearsonr
import math

'''
Compare 2 .fpkm_tracking files output by Cufflinks (run with the -G option)
'''

threshold = 0

fpkmsA = []
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
            fpkmsA.append((row[label_id], float(row[fpkm_id])))

fpkmsB = []
with open(sys.argv[2], 'r') as f:
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
            #print(row)
            fpkm = float(row[fpkm_id])
            fpkmsB.append((row[label_id], float(row[fpkm_id])))

if not len(fpkmsA) == len(fpkmsB):
    print('Lengths not equal!')
    print(len(fpkmsA))
    print(len(fpkmsB))
    #exit()
else:
   print('%d fpkms' % len(fpkmsA))

fpkmsA.sort()
fpkmsB.sort()
#for i in range(100):
#    print(fpkmsA[i][0] + '\t' + fpkmsB[i][0])
#print('')

namesA = [d[0] for d in fpkmsA]
valsA = [d[1] for d in fpkmsA]
namesB = [d[0] for d in fpkmsB]
valsB = [d[1] for d in fpkmsB]
matchesA = []
matchesB = []
unmatchedA = []
for i in range(len(namesA)):
    n = namesA[i]
    if n in namesB:
        matchesA.append(valsA[i])
        j = namesB.index(n)
        matchesB.append(valsB[j])
        #if valsA[i] < 300 and valsB[j] > 1000:
        #    print(n)
        #    print(valsA[i])
        #    print(valsB[j])
        #    exit()
        del namesB[j]
        del valsB[j]
        del fpkmsB[j]
    else:
        unmatchedA.append(fpkmsA[i])

#print('FBtr0100187' in fpkmsB)
print('%d total' % len(fpkmsA))
#print(len(unmatchedA))
#print(len(fpkmsB))
#print(unmatchedA[:10])
#print(fpkmsB[:10])

with open('quant_cufflinks.txt', 'w') as f:
    for i in range(len(matchesA)):
        f.write(str(matchesA[i]) + '\t' + str(matchesB[i]) + '\n')

#print('Comparing %d transcripts' % len(fpkmsA))
#p = pearsonr(matchesA, matchesB)
#print('Pearson correlation coefficient:\t%f, %f' % (p[0], p[1]))
#print('')

n = len(matchesA)
error = sum([(matchesA[i]-matchesB[i]) ** 2 for i in range(n)]) / n
print('RMSE:\t%f' % (math.sqrt(error)))
