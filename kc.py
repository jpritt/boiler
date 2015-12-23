#! /usr/bin/env python
import sys
import argparse
from itertools import izip
import time

'''
Calculate the KC score for the reference sequence and alignment given
'''

#k = 24
sumKmersWeighted = dict()

def buildProfile(transcripts, k):
    ''' Create frequency profile
    '''

    # frequency profile
    profile = dict()

    # denominator is the same for all kmers
    denom = 0.0
    for t in transcripts:
        denom += t.numKmers * t.expressed


    # loop through all possible kmers
    last_kmer = 'T' * k
    kmer = 'A' * k
    while kmer < last_kmer:
        numer = 0

        for t in transcripts:
            if kmer in t.kmers:
                #numer += t.kmers[kmer] * t.expressed
                numer += t.kmersWeighted[kmer]

        if numer > 0:
            profile[kmer] = float(numer) / denom

        kmer = increment_kmer(kmer)

    return profile

def computeWKRandICR(assembly, sequence, transcripts, dataLength, k):
    ''' Compute the WKR from the frequency profile and cufflinks GTF (assembly) file

        assembly: Cufflinks GTF file containing assembled transcript
        sequence: Genome sequence (read from fasta file)
        transcripts: List of all transcripts in reference with expression > 0
        dataLength: Length in nucleotides of input data
    '''
    with open(assembly, 'r') as f:
        # WKR denominator is the same for all kmers
        denom = 0.0
        for t in transcripts:
            denom += t.numKmers * t.expressed

        # Only process unique kmers
        kmersProcessed = set()

        WKR = 0.0

        # Length of assembly
        lenA = 0.0

        for line in f:
            #start_time = time.time()

            #n += 1
            #if (n%1000) == 0:
            #    print('%d / %d:' % (n, numLines))

            row = line.rstrip().split('\t')

            #split_time = time.time()
            #x = split_time-start_time
            #sys.stdout.write('Split (%f), ' % x)

            if len(row) >= 5 and row[2] == 'exon':
                chrom = row[0]
                start = int(row[3])
                end = int(row[4])

                lenA += end-start

                for i in xrange(start+k, end+1):
                    #loop_time = time.time()

                    kmer = sequence[chrom][i-k:i]

                    if not kmer in kmersProcessed:

                        '''
                        # WKR numerator
                        numer = 0

                        for t in transcripts:
                            if kmer in t.kmersSet:
                                numer += t.kmersWeighted[kmer]
                                #numer += t.kmers[kmer] * t.expressed

                        if numer > 0:
                            WKR += float(numer) / denom
                        '''

                        if kmer in sumKmersWeighted:
                            WKR += float(sumKmersWeighted[kmer]) / denom

                        kmersProcessed.add(kmer)
                #processed_time = time.time()
                #x = processed_time-split_time
                #sys.stdout.write('Processed %d (%d new) kmers (%f), ' % (countOld+countNew, countNew, x))
            #else:
            #    sys.stdout.write('Not an exon, ')
            #end_time = time.time()
            #x = end_time-start_time
            #print('\n  => %fs total' % x)

            #if n >= 100:
            #    exit()
        ICR = float(lenA) / float(dataLength)

        return WKR, ICR

def countDataLength(data):
    ''' Count total number of nucleotides in read data

        data:
    '''

    with open(data, 'r') as f:
        length = 0
        def grouped(iterable, n):
            "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
            return izip(*[iter(iterable)]*n)

        for w,x,y,z in grouped(f, 4):
           length += len(x.rstrip())
        #print 'Data length: ' + str(length)
        return length

def increment_kmer(kmer):
    ''' Return the next kmer incrementally
    '''
    newKmer = ''
    i = len(kmer)-1
    while kmer[i] == 'T':
        newKmer = 'A' + newKmer
        i -= 1
        if i < 0:
            return kmer
    if kmer[i] == 'A':
        newKmer = 'C' + newKmer
    elif kmer[i] == 'C':
        newKmer = 'G' + newKmer
    elif kmer[i] == 'G':
        newKmer = 'T' + newKmer
    #print kmer + ' --> ' + kmer[:i] + newKmer
    return kmer[:i] + newKmer

def parsePro(filename):
    ''' Return a dictionary with transcript id (e.g. 0300689) pointing to coverage level
    '''
    threshold = 0

    transcripts = dict()
    with open(filename, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            if len(row) < 8:
                continue

            tag = row[1]
            sep1 = row[0].find(':')
            sep2 = row[0].find('-', sep1)
            sep3 = row[0].find('W', sep2)
            chrom = row[0][:sep1]
            start = int(row[0][sep1+1:sep2])
            end = int(row[0][sep2+1:sep3])

            # expressed number
            number = float(row[5])

            if number > threshold:
                transcripts[tag] = Transcript(chrom, start, end, number, tag)
    return transcripts

def readGenome(genome):
    with open(genome, 'r') as f:
        chromosomes = dict()

        currName = ''
        currSeq = ''

        for l in f:
            if l[0] == '>':
                if len(currName) > 0:
                    chromosomes[currName] = currSeq

                currName = l[1:].rstrip()
                currSeq = ''
            else:
                currSeq += l.rstrip()
        # Add last chromosome
        if len(currName) > 0:
            chromosomes[currName] = currSeq

        return chromosomes

def updateTranscripts(transcripts, refGTF, sequence, k):
    with open(refGTF, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            if len(row) < 9:
                continue

            transcriptIdIndex = row[8].find('transcript_id')
            transcriptIdStart = row[8].find('"', transcriptIdIndex) + 1
            transcriptIdEnd = row[8].find('"', transcriptIdStart)
            transcriptId = row[8][transcriptIdStart : transcriptIdEnd]


            #if row[1] == 'protein_coding' and row[2] == 'exon' and transcriptId in transcripts:
            if row[2] == 'exon' and transcriptId in transcripts:
                transcripts[transcriptId].exons.append( (int(row[3]), int(row[4])) )

    transcripts = transcripts.values()
    #i = 0
    #numTranscripts = len(transcripts)
    for t in transcripts:
        #i += 1
        #if (i%1000) == 0:
        #    print 'Adding sequence %d / %d' % (i, numTranscripts)

        t.addSequence(sequence[t.chrom], k)
    return transcripts

class Transcript:
    '''
    '''
    def __init__(self, chrom, start, end, expressedNumber, name):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.expressed = expressedNumber
        self.name = name

        self.exons = []

    def addSequence(self, sequence, k):
        ''' Add nucleotide sequence covered by this transcript
        '''
        self.sequence = ''

        for e in self.exons:
            self.sequence += sequence[e[0]:e[1]]

        self.numKmers = len(self.sequence) - k + 1
        #self.kmers = dict()

        # number of each kmer times expression level of transcript
        #self.kmersWeighted = dict()

        for i in xrange(k, len(self.sequence)+1):
            kmer = self.sequence[i-k:i]
            '''
            if kmer in self.kmersWeighted:
                #self.kmers[kmer] += 1
                self.kmersWeighted[kmer] += self.expressed
            else:
                #self.kmers[kmer] = 1
                self.kmersWeighted[kmer] = self.expressed
            '''
            if kmer in sumKmersWeighted:
                sumKmersWeighted[kmer] += self.expressed
            else:
                sumKmersWeighted[kmer] = self.expressed

        #for kmer,value in self.kmers.items():
        #    self.kmersWeighted[kmer] = value * self.expressed

        #self.kmersSet = set(self.kmers.keys())
        #self.kmersSet = set(self.kmersWeighted.keys())

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--refGTF', type=str, required=True, 
        help='Full path of reference GTF file')
    parser.add_argument('--refPRO', type=str, required=True, 
        help='Full path of reference PRO file with transcript frequencies')
    parser.add_argument('--sequence', type=str, required=True, 
        help='Full path of fasta sequence file')
    parser.add_argument('--assembly', type=str, required=True, 
        help='Full path of assembled Cufflinks GTF file')
    parser.add_argument('--data', type=str, required=True,
        help='Full path of read data used for assembly')
    parser.add_argument('--kmer', type=str, required=True,
        help='k-mer length')
    
    
    args = parser.parse_args(sys.argv[1:])

    k = int(args.kmer)
    transcripts = parsePro(args.refPRO)
    sequence = readGenome(args.sequence)
    transcripts = updateTranscripts(transcripts, args.refGTF, sequence, k)

    dataLength = countDataLength(args.data)

    assembly = args.assembly
    WKR, ICR = computeWKRandICR(assembly, sequence, transcripts, dataLength, k)
    print('WKR = %f' % WKR)
