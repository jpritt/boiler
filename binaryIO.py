#! /usr/bin/env python3
import sys
import math
import bucket
import cross_bundle_bucket
import huffman
import time
from struct import *

# Reading/writing formats
formats = {1:'B', 2:'H', 4:'I', 8:'Q'}
byteVal = int(pow(2,8))

def valToBinary(numBytes, val):
    ''' Convert a value to a byte string which can then be written to a file
    '''

    return (val).to_bytes(numBytes, byteorder='big')

def binaryToVal(s, numBytes, start=0):
    ''' Convert a byte string from a file to an integer value
    '''

    return int.from_bytes(s[start:start+numBytes], byteorder='big'), numBytes+start

def writeVal(f, numBytes, val):
    f.write(valToBinary(numBytes, val))

def readVal(f, numBytes):
    v,_ = binaryToVal(f.read(numBytes), numBytes)
    return v

def writeChroms(chroms):
    names = []
    lens = []
    maxLen = 0
    for k,v in chroms.items():
        names.append(k)
        lens.append(v)
        if v > maxLen:
            maxLen = v

    # find length in bytes to fit all numbers
    numBytes = findNumBytes(maxLen)

    # Write keys
    s = bytes('\t'.join(names) + '\n', 'UTF-8')

    # Write values
    s += valToBinary(1, numBytes)
    for l in lens:
        s += valToBinary(numBytes, l)

    return s

def readChroms(f):
    line = f.readline().decode('ascii')
    names = line.rstrip().split('\t')

    numBytes = readVal(f, 1)
    chroms = dict()
    for n in names:
        chroms[n] = readVal(f, numBytes)
    return chroms

def writeDict(d):
    ''' Generic method to write a dictionary in tuple form with string keys and integer values '''

    keys = [v[0] for v in d]
    vals = [v[1] for v in d]

    # find length in bytes to fit all numbers
    maxVal = max(vals)
    valBytes = findNumBytes(maxVal)

    # Write keys
    s = bytes(','.join(keys) + '\n', 'UTF-8')

    # Add length of keys string to beginning
    s = valToBinary(4, len(s)) + s

    # Write values
    s += valToBinary(1, valBytes)
    for v in vals:
        s += valToBinary(valBytes, v)
    return s

def readDict(s, start=0):
    ''' Generic method to read a dictionary with string keys and integer values in the format written by writeDict() '''

    d = []

    # Read length of keys string
    keyLen, start = binaryToVal(s, 4, start)

    # Read keys
    keyStr = s[:keyLen].decode('ascii').rstrip()
    keys = keyStr.split(',')
    start += keyLen

    # Read number of bytes for each value
    valBytes, start = binaryToVal(s, 1, start)

    # Read values
    for k in keys:
        v, start = binaryToVal(s, valBytes, start)
        d.append((k,v))

    return d

def writeList(vals):
    ''' Generic method to write a list of integers '''

    maxVal = max(vals)

    # find length in bytes to fit all numbers
    valBytes = findNumBytes(maxVal)

    s = valToBinary(1, valBytes)
    s += valToBinary(4, len(vals))
    for v in vals:
        s += valToBinary(valBytes, v)
    return s

def readList(s, start=0):
    ''' Generic method to read a list of integers in the format written by writeList() '''
    valBytes, start = binaryToVal(s, 1, start)
    numVals, start = binaryToVal(s, 4, start)

    vals = [0] * numVals
    for i in range(numVals):
        vals[i], start = binaryToVal(s, valBytes, start)
    return vals, start

def readListFromFile(f):
    ''' Generic method to read a list of integers in the format written by writeList() '''
    valBytes = readVal(f, 1)
    numVals = readVal(f, 4)

    vals = [0] * numVals
    for i in range(numVals):
        vals[i] = readVal(f, valBytes)
    return vals

def writeExons(exons):
    ''' Write the list of exons. More specific version of writeList() above '''
    maxExon = exons[-1]

    # find length in bytes to fit all numbers
    exonBytes = findNumBytes(maxExon)

    # length in bytes to fit all exon ids
    numExonBytes = findNumBytes(len(exons))

    s = valToBinary(1, exonBytes)
    #s += valToBinary(1, numExonBytes)
    s += valToBinary(numExonBytes, len(exons))
    for e in exons:
        s += valToBinary(exonBytes, e)

    return s, numExonBytes

def readExons(f, numExonBytes):
    ''' Read the list of exons '''
    exonBytes = readVal(f, 1)
    numExons = readVal(f, numExonBytes)

    exons = []
    for _ in range(numExons):
        exons.append(readVal(f, exonBytes))
    return exons

def writeClusters(clusters):
    '''

    :param clusters: 2d list of splice sites, grouped by cluster
    :return:
    '''

    maxExon = clusters[-1][-1]

    # find length in bytes to fit all numbers
    exonBytes = findNumBytes(maxExon)

    # length in bytes to fit all exon ids
    max_exons = 0
    for c in clusters:
        if len(c) > max_exons:
            max_exons = len(c)
    exonIdBytes = findNumBytes(max_exons)

    s = valToBinary(1, exonBytes)
    s += valToBinary(1, exonIdBytes)
    s += valToBinary(exonBytes, len(clusters))
    for c in clusters:
        s += valToBinary(exonIdBytes, len(c))
        for e in c:
            s += valToBinary(exonBytes, e)

    return s

def readClusters(f):
    exonBytes = readVal(f, 1)
    exonIdBytes = readVal(f, 1)

    clusters = []
    num_c = readVal(f, exonBytes)

    for i in range(num_c):
        num_e = readVal(f, exonIdBytes)
        splice_sites = [0] * num_e
        for j in range(num_e):
            splice_sites[j] = readVal(f, exonBytes)
        clusters.append(splice_sites)

    return clusters

def writeJunctionsList(junctions, exonBytes):
    s = valToBinary(4, len(junctions))
    s += valToBinary(1, exonBytes)

    for j in junctions:
        j = j.split('\t')

        # First byte: XS value, 0 = None, -1 = -, 1 = +
        # Second byte: Number of exons in junction
        # These could probably be combined into a single byte, but I'm keeping them separate just in case
        #   A junction ever contains more than 64 subexons
        xs = int(j[-2])
        s += pack('b', xs)
        s += pack('B', len(j)-2)

        '''
        if j[-2] == '-':
            s += pack('b', 2-len(j))
        else:
            s += pack('b', len(j)-2)
        '''
        
        for e in j[:-2]:
            s += valToBinary(exonBytes, int(e))
        s += valToBinary(2, int(j[-1]))

    return s

def readJunctionsList(s, start=0):
    numJunctions, start = binaryToVal(s, 4, start)
    exonBytes, start = binaryToVal(s, 1, start)

    junctions = [[]] * numJunctions
    for j in range(numJunctions):
        # Read xs value
        v = unpack_from('b', s[start:start+1])[0]
        xs = None
        if v == -1:
            xs = '-'
        elif v == 1:
            xs = '+'

        # Read number of exons
        num_exons = unpack_from('B', s[start+1:start+2])[0]

        start += 2

        '''
        v = unpack_from('b', s[start:start+1])[0]

        start += 1

        # Read XS value
        if v < 0:
            xs = '-'
            numExons = -v
        else:
            xs = '+'
            numExons = v
        '''

        # Read exons
        juncExons = [0] * num_exons
        for e in range(num_exons):
            juncExons[e], start = binaryToVal(s, exonBytes, start)

        # Read NH value
        NH, start = binaryToVal(s, 2, start)

        # Create junction string
        junctions[j] = juncExons + [xs, NH]

    return junctions, exonBytes, start

def writeJunction(readLenBytes, junc):
    s = writeCov(junc.coverage)
    s += writeLens(readLenBytes, junc.unpairedLens)

    # Find max number of bytes needed to encode fragment lengths
    if len(junc.pairedLens) == 0:
        fragLenBytes = 1
    else:
        fragLenBytes = findNumBytes(max(junc.pairedLens))
    s += valToBinary(1, fragLenBytes)
    s += writeLens(fragLenBytes, junc.pairedLens)
    if len(junc.pairedLens) > 0:
        s += writeLens(readLenBytes, junc.lensLeft)
        if len(junc.lensLeft) > 0:
            s += writeLens(readLenBytes, junc.lensRight)

    return s

def readJunction(s, junc, readLenBytes, start=0):
    junc.coverage, start = readCov(s, start)
    junc.unpairedLens, start = readLens(s, readLenBytes, start)

    fragLenBytes, start = binaryToVal(s, 1, start)
    junc.pairedLens, start = readLens(s, fragLenBytes, start)

    if len(junc.pairedLens) > 0:
        junc.lensLeft, start = readLens(s, readLenBytes, start)
        if len(junc.lensLeft) > 0:
            junc.lensRight, start = readLens(s, readLenBytes, start)

    return junc, start

def writeCrossBundleBucket(bundleIdBytes, readLenBytes, bucket):
    s = valToBinary(bundleIdBytes, bucket.bundleA)
    s += writeList(bucket.exonIdsA)
    s += valToBinary(bundleIdBytes, bucket.bundleB)
    s += writeList(bucket.exonIdsB)
    s += pack('b', bucket.XS)
    s += valToBinary(1, bucket.NH)

    s += writeCov(bucket.coverage)

    # Find max number of bytes needed to encode fragment lengths
    if len(bucket.pairedLens) == 0:
        fragLenBytes = 1
    else:
        fragLenBytes = findNumBytes(max(bucket.pairedLens))
    s += valToBinary(1, fragLenBytes)
    s += writeLens(fragLenBytes, bucket.pairedLens)
    if len(bucket.pairedLens) > 0:
        s += writeLens(readLenBytes, bucket.lensLeft)
        if len(bucket.lensLeft) > 0:
            s += writeLens(readLenBytes, bucket.lensRight)

    return s

def readCrossBundleBucket(s, bundleIdBytes, readLenBytes, start=0):
    # Read bucket exon information
    bundleA, start = binaryToVal(s, bundleIdBytes, start)
    exonIdsA, start = readList(s, start)
    bundleB, start = binaryToVal(s, bundleIdBytes, start)
    exonIdsB, start = readList(s, start)

    v = unpack_from('b', s[start:start+1])[0]
    XS = None
    if v == -1:
        XS = '-'
    elif v == 1:
        XS = '+'
    start += 1

    NH, start = binaryToVal(s, 1, start)

    bucket = cross_bundle_bucket.CrossBundleBucket(bundleA, exonIdsA, bundleB, exonIdsB)
    bucket.XS = XS
    bucket.NH = NH

    coverage, start = readCov(s, start)
    length = 0
    for c in coverage:
        length += c[1]
    bucket.set_length(length)
    bucket.coverage = coverage

    fragLenBytes, start = binaryToVal(s, 1, start)
    bucket.pairedLens, start = readLens(s, fragLenBytes, start)

    if len(bucket.pairedLens) > 0:
        bucket.lensLeft, start = readLens(s, readLenBytes, start)
        if len(bucket.lensLeft) > 0:
            bucket.lensRight, start = readLens(s, readLenBytes, start)

    return bucket, start

def skipJunction(s, readLenBytes, start=0):
    _, start = readCov(s, start)
    _, start = readLens(s, readLenBytes, start)

    fragLenBytes, start = binaryToVal(s, 1, start)
    lens, start = readLens(s, fragLenBytes, start)

    if len(lens) > 0:
        lens, start = readLens(s, readLenBytes, start)
        if len(lens) > 0:
            _, start = readLens(s, readLenBytes, start)

    return start

'''
def skipJunction(s, readLenBytes, start=0):
    fragLenBytes, start = binaryToVal(s, 1, start)

    lens, start = readLens(s, readLenBytes, start)
    lens, start = readLens(s, fragLenBytes, start)
    if len(lens) > 0:
        lens, start = readLens(s, readLenBytes, start)
        if len(lens) > 0:
            _, start = readLens(s, readLenBytes, start)

    _, start = readCov(s, start)

    return start
'''

def writeLens(lenBytes, lens):
    # Write number of lengths
    s = valToBinary(2, len(lens))

    if len(lens) > 0:
        # Write number of bytes for each frequency
        freqBytes = findNumBytes(max(lens.values()))
        s += valToBinary(1, freqBytes)

        for l,n in lens.items():
            s += valToBinary(lenBytes, l)
            s += valToBinary(freqBytes, n)

    return s

def readLens(s, lenBytes, start=0):
    # Read number of lengths
    numLens, start = binaryToVal(s, 2, start)
    lens = dict()

    if numLens > 0:
        # Read number of bytes for each frequency
        freqBytes, start = binaryToVal(s, 1, start)

        for _ in range(numLens):
            l, start = binaryToVal(s, lenBytes, start)
            n, start = binaryToVal(s, freqBytes, start)
            lens[l] = n

    return lens, start

def skipLens(s, lenBytes, start=0):
    numLens, start = binaryToVal(s, 2, start)

    if numLens == 0:
        return False, start
    else:
        freqBytes, start = binaryToVal(s, 1, start)
        start += (lenBytes+freqBytes) * numLens
        return True, start

def writePairs(pairs, numBytes=3):
    '''
    Writes a list of index pairs as a simple list
    :param pairs: List of tuples corresponding to paired reads
    :param numBytes: Number of bytes to use for each read id
    :return: A binary string encoding the pairs as a simple concatenated list
    '''

    numBytes = 3

    s = valToBinary(numBytes, len(pairs))
    for p in pairs:
        s += valToBinary(numBytes, p[0]) + valToBinary(numBytes, p[1])
    return s

def readPairs(s, start=0, numBytes=3):
    '''
    :param s: Binary string to read from
    :param numBytes: Number of bytes used to encode each index
    :param start: Starting index to read from in s
    :return: A list of tuples containing paired indices, and the starting index for the remainder of the string
    '''

    numBytes = 3

    length, start = binaryToVal(s, numBytes, start)
    pairs = [None] * length
    for i in range(length):
        a, start = binaryToVal(s, numBytes, start)
        b, start = binaryToVal(s, numBytes, start)
        pairs[i] = (a,b)
    return pairs, start

def writeCov(cov):
    maxLen = 0
    maxCov = 0
    for c in cov:
        if c[0] > maxCov:
            maxCov = c[0]
        if c[1] > maxLen:
            maxLen = c[1]

    # find length in bytes to fit all lengths
    lenBytes = findNumBytes(maxLen)

    # find length in bytes to fit all covs
    covBytes = findNumBytes(maxCov)

    # Write size of length and cov in bytes
    s = valToBinary(1, lenBytes)
    s += valToBinary(1, covBytes)

    # Write the length of the vector
    s += valToBinary(4, len(cov))

    for c in cov:
        s += valToBinary(covBytes, c[0])
        s += valToBinary(lenBytes, c[1])

    return s

def readCov(s, start=0):
    # Read size of length and cov in bytes
    lenBytes, start = binaryToVal(s, 1, start)
    covBytes, start = binaryToVal(s, 1, start)

    # Read the length of the vector
    lenCov, start = binaryToVal(s, 4, start)

    cov = []
    for _ in range(lenCov):
        c, start = binaryToVal(s, covBytes, start)
        l, start = binaryToVal(s, lenBytes, start)
        cov.append([c,l])
        #cov += [c] * l

    return cov, start

def writeCovHuffman(cov, huffmanIndex):
    maxLen = 0
    for c in cov:
        if c[1] > maxLen:
            maxLen = c[1]

    # find length in bytes to fit all lengths
    lenBytes = findNumBytes(maxLen)
    #print('%d cov length bytes' % lenBytes)

    # Write size of length in bytes
    s = valToBinary(1, lenBytes)

    # Write the length of the vector
    s += valToBinary(2, len(cov))

    # Write lengths
    for c in cov:
        #print('Writing %d to %d bytes' % (c[1], lenBytes))
        s += valToBinary(lenBytes, c[1])

    # Write covs
    encoding = huffman.bitsToBytes(huffman.encode([c[0] for c in cov], huffmanIndex))
    #print('Writing %d to 4 bytes' % len(encoding))
    s += encoding
    #print('Finished writing')
    return s

def readCovHuffman(s, huffmanTree, start=0):
    # Read size of length and cov in bytes
    lenBytes, start = binaryToVal(s, 1, start)

    # Read the length of the vector
    lenCov, start = binaryToVal(s, 2, start)

    # Read lengths
    lengths = [0] * lenCov
    for i in range(lenCov):
        l, start = binaryToVal(s, lenBytes, start)
        lengths[i] = l

    # Read coverage 
    #l = binaryToVal(s, 4)
    #c, s = huffman.decode(huffman.bytesToBits(s, lenCov), huffmanTree, lenCov)

    c, start = huffman.decode(s, huffmanTree, lenCov, start)

    cov = []
    for i in range(lenCov):
        cov += [c[i]] * lengths[i]

    return cov, start


def writeHuffmanIndex(index):
    # Write number of values
    s = valToBinary(4, len(index))

    # Write the number of bytes used to write each value
    maxVal = max(index.keys())
    valBytes = findNumBytes(maxVal)
    s += valToBinary(1, valBytes)

    for val,bits in index.items():
        # Write the value
        s += valToBinary(valBytes, val)

        # Write the number of bits in the encoding
        s += valToBinary(1, len(bits))

        # Write the bits, padded by zeros
        s += huffman.bitsToBytes(bits)

    return s

def readHuffmanIndex(s, start=0):
    # Read number of values
    numVals, start = binaryToVal(s, 4, start)

    # Read the number of bytes used to write each value
    valBytes, start = binaryToVal(s, 1, start)

    index = dict()
    for _ in range(numVals):
        # Read the value
        val, start = binaryToVal(s, valBytes, start)

        # Read the number of bits in the encoding
        numBits, start = binaryToVal(s, 1, start)

        # Read the bits
        numBytes = math.ceil(numBits / 8)
        bits = huffman.bytesToBits(s[:numBytes])
        start += numBytes

        index[val] = bits[:numBits]
    return index, start

def RLE(vector):
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

def findNumBytes(val):
    ''' Return the number of bytes needed to encode the given value.
        Will only return 1, 2, 4, or 8 (the number of bytes that struct.pack supports)
    '''

    if val <= 0:
        return 1
    else:
        return math.ceil(math.log(val+1, 2) / 8.0)

