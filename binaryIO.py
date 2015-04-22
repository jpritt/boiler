#! /usr/bin/env python3
import sys
import math
import junction
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
    '''
    s = b''
    while numBytes > 0:
        v = val % byteVal
        s += pack('B', v)
        val = val >> 8
        numBytes -= 1

    if val > 0:
        print('Error! %d is too large to fit in %d bytes!' % (val, numBytes))
        exit()

    return s
    '''

def binaryToVal(s, numBytes, start=0):
    ''' Convert a byte string from a file to an integer value
    '''

    return int.from_bytes(s[start:start+numBytes], byteorder='big'), numBytes+start

    '''
    # Compute value for each character and add it to val
    val = 0
    #b = unpack('B'*numBytes, s[:numBytes])
    for i in range(numBytes):
        val += s[start+i] << (8*i)
    return val, start+numBytes
    '''

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

def writeJunctionsList(junctions, exonBytes):
    s = valToBinary(4, len(junctions))
    s += valToBinary(1, exonBytes)

    for j in junctions:
        j = j.split('\t')

        # First value: Number of exons in junction
        # Second value: 0 for negative XS, 1 for positive
        if j[-2] == '-':
            s += pack('b', 2-len(j))
        else:
            s += pack('b', len(j)-2)
        
        for e in j[:-2]:
            s += valToBinary(exonBytes, int(e))
        s += valToBinary(2, int(j[-1]))

    return s

def readJunctionsList(s, start=0):
    numJunctions, start = binaryToVal(s, 4, start)
    exonBytes, start = binaryToVal(s, 1, start)

    junctions = [[]] * numJunctions
    for j in range(numJunctions):
        # Read number of exons
        v = unpack_from('b', s[start:start+1])[0]
        start += 1
        numExons = abs(v)

        # Read XS value
        if v < 0:
            xs = '-'
        else:
            xs = '+'
        

        # Read exons
        juncExons = [0] * numExons
        for e in range(numExons):
            juncExons[e], start = binaryToVal(s, exonBytes, start)

        # Read NH value
        NH, start = binaryToVal(s, 2, start)

        # Create junction string
        junctions[j] = juncExons + [xs, NH]

    return junctions, exonBytes, start

def writeJunction(readLenBytes, junc, huffmanIndex=None):
    # More cost efficient to calculate/save the number of bytes needed for fragments for each junction
    maxFragLen = 0
    for l in junc.pairedLens.keys():
        if l > maxFragLen:
            maxFragLen = l
    fragLenBytes = findNumBytes(maxFragLen)
    s = valToBinary(1, fragLenBytes)

    s += writeLens(readLenBytes, junc.unpairedLens)
    s += writeLens(fragLenBytes, junc.pairedLens)
    if len(junc.pairedLens) > 0:
        s += writeLens(readLenBytes, junc.lensLeft)
        if len(junc.lensLeft) > 0:
            s += writeLens(readLenBytes, junc.lensRight)

    if not huffmanIndex == None:
        #s += writeCovHuffman(RLE(junc.coverage), huffmanIndex)
        s += writeCovHuffman(junc.coverage, huffmanIndex)
    else:
        #s += writeCov(RLE(junc.coverage))
        s += writeCov(junc.coverage)

    return s

def readJunction(s, junc, readLenBytes, start=0, huffmanTree=None):
    fragLenBytes, start = binaryToVal(s, 1, start)

    junc.unpairedLens, start = readLens(s, readLenBytes, start)
    junc.pairedLens, start = readLens(s, fragLenBytes, start)
    if len(junc.pairedLens) > 0:
        junc.lensLeft, start = readLens(s, readLenBytes, start)
        if len(junc.lensLeft) > 0:
            junc.lensRight, start = readLens(s, readLenBytes, start)

    if not huffmanTree == None:
        junc.coverage, start = readCovHuffman(s, huffmanTree, start)
    else:
        junc.coverage, start = readCov(s, start)

    return junc, start

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
    s += valToBinary(2, len(cov))

    for c in cov:
        s += valToBinary(covBytes, c[0])
        s += valToBinary(lenBytes, c[1])

    return s

def readCov(s, start=0):
    # Read size of length and cov in bytes
    lenBytes, start = binaryToVal(s, 1, start)
    covBytes, start = binaryToVal(s, 1, start)

    # Read the length of the vector
    lenCov, start = binaryToVal(s, 2, start)

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

