#! /usr/bin/env python
import sys
from scipy.stats.stats import pearsonr
import math

'''
Compare 2 transcript files output by Stringtie (run with the -G and -e options)
'''

threshold = 0

fpkmsA = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        row = line.rstrip().split('\t')

        if row[2] == 'transcript':
            desc = row[8]
            i = desc.index('reference_id')
            start = desc.index('"', i+1)
            end = desc.index('"', start+1)
            name = desc[start+1:end]

            i = desc.index('FPKM')
            start = desc.index('"', i+1)
            end = desc.index('"', start+1)
            fpkm = float(desc[start+1:end])
            fpkmsA.append((name, fpkm))

fpkmsB = []
with open(sys.argv[2], 'r') as f:
    for line in f:
        row = line.rstrip().split('\t')

        if row[2] == 'transcript':
            desc = row[8]
            i = desc.index('reference_id')
            start = desc.index('"', i+1)
            end = desc.index('"', start+1)
            name = desc[start+1:end]

            i = desc.index('FPKM')
            start = desc.index('"', i+1)
            end = desc.index('"', start+1)
            fpkm = float(desc[start+1:end])
            fpkmsB.append((name, fpkm))

fpkmsA.sort()
fpkmsB.sort()

lenA = len(fpkmsA)
lenB = len(fpkmsB)

'''
if not len(fpkmsA) == len(fpkmsB):
    print('Lengths not equal!')
    print(len(fpkmsA))
    print(len(fpkmsB))
    #exit()
else:
   print('%d fpkms' % len(fpkmsA))
'''

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
        if valsA[i] > 1400 and valsA[i] < 1600:
            print(n + '\t' + str(valsA[i]) + '\t' + str(valsB[j]))
        del namesB[j]
        del valsB[j]
        del fpkmsB[j]
    else:
        matchesA.append(valsA[i])
        matchesB.append(0)
        unmatchedA.append(fpkmsA[i])

for v in valsB:
    matchesA.append(0)
    matchesB.append(v)

print('%d / %d unmatched from A' % (len(unmatchedA), lenA))
print('%d / %d unmatched from B' % (len(fpkmsB), lenB))
#print(unmatchedA[:10])
#print(fpkmsB[:10])

print('%d total in vector' % len(matchesA))

with open('quant_stringtie.txt', 'w') as f:
    for i in range(len(matchesA)):
        f.write(str(matchesA[i]) + '\t' + str(matchesB[i]) + '\n')

#print('Comparing %d transcripts' % len(fpkmsA))
#p = pearsonr(matchesA, matchesB)
#print('Pearson correlation coefficient:\t%f, %f' % (p[0], p[1]))

n = len(matchesA)
error = sum([(matchesA[i]-matchesB[i]) ** 2 for i in range(n)]) / n
print('RMSE:\t%f' % (math.sqrt(error)))
