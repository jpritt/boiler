#! /usr/bin/env python3
import sys
import math
import junction
import huffman
from struct import *

# Reading/writing formats
formats = {1:'B', 2:'H', 4:'I', 8:'Q'}
byteVal = int(pow(2,8))

def valToBinary(numBytes, val):
    ''' Convert a value to a byte string which can then be written to a file
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

def binaryToVal(s, numBytes):
    ''' Convert a byte string from a file to an integer value
    '''

    # Compute value for each character and add it to val
    val = 0
    #b = unpack('B'*numBytes, s[:numBytes])
    for i in range(numBytes):
        val += s[i] << (8*i)
    return val, s[numBytes:]

def writeVal(f, numBytes, val):
    f.write(valToBinary(numBytes, val))

def readVal(f, numBytes):
    v,_ = binaryToVal(f.read(numBytes), numBytes)
    return v
    '''
    val = 0
    bs = unpack('B'*numBytes, f.read(numBytes))

    # Compute value for each character and add it to val
    for i in range(numBytes):
        val += bs[i] << (8*i)

    return val
    #return unpack(formats[numBytes], f.read(numBytes))[0]
    '''

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

def readDict(s):
    ''' Generic method to read a dictionary with string keys and integer values in the format written by writeDict() '''

    d = []

    # Read length of keys string
    keyLen, s = binaryToVal(s, 4)

    # Read keys
    keyStr = s[:keyLen].decode('ascii').rstrip()
    keys = keyStr.split(',')
    s = s[keyLen:]

    # Read number of bytes for each value
    valBytes, s = binaryToVal(s, 1)

    # Read values
    for k in keys:
        v, s = binaryToVal(s, valBytes)
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

def readList(s):
    ''' Generic method to read a list of integers in the format written by writeList() '''
    valBytes, s = binaryToVal(s, 1)
    numVals, s = binaryToVal(s, 4)

    vals = [0] * numVals
    for i in range(numVals):
        vals[i], s = binaryToVal(s, valBytes)
    return vals, s

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
    #numExonBytes = readVal(f, 1)
    numExons = readVal(f, numExonBytes)

    exons = []
    for _ in range(numExons):
        exons.append(readVal(f, exonBytes))
    return numExonBytes, exons

def writeJunctionsList(junctions, exonBytes):
    s = valToBinary(4, len(junctions))
    s += valToBinary(1, exonBytes)

    for j in junctions:
        #print(j)
        j = j.split('\t')

        # First value: Number of exons in junction, +/- depending on XS value
        if j[-2] == '-':
            s += pack('b', 2-len(j))
        else:
            s += pack('b', len(j)-2)
        for e in j[:-2]:
            s += valToBinary(exonBytes, int(e))
        s += valToBinary(2, int(j[-1]))

    return s

def readJunctionsList(s):
    numJunctions, s = binaryToVal(s, 4)
    exonBytes, s = binaryToVal(s, 1)

    junctions = [''] * numJunctions
    for j in range(numJunctions):
        # Read number of exons
        v = unpack_from('b', s)[0]
        numExons = abs(v)

        # Read XS value
        if v < 0:
            xs = '-'
        else:
            xs = '+'

        # Read exons
        s = s[1:]
        juncExons = [0] * numExons
        for e in range(numExons):
            juncExons[e], s = binaryToVal(s, exonBytes)

        # Read NH value
        NH, s = binaryToVal(s, 2)

        # Create junction string
        junctions[j] = '\t'.join([str(e) for e in juncExons]) + '\t' + xs + '\t' + str(NH)

    return junctions, exonBytes, s

def writeJunction(readLenBytes, junc, huffmanIndex=None):
    # More cost efficient to calculate/save the number of bytes needed for fragments for each junction
    maxFragLen = 0
    for l in junc.readLens.keys():
        if l > maxFragLen:
            maxFragLen = l
    fragLenBytes = findNumBytes(maxFragLen)
    s = valToBinary(1, fragLenBytes)

    s += writeLens(fragLenBytes, junc.readLens)
    if len(junc.readLens) > 0:
        s += writeLens(readLenBytes, junc.lensLeft)
        if len(junc.lensLeft) > 0:
            s += writeLens(readLenBytes, junc.lensRight)

    if not huffmanIndex == None:
        s += writeCovHuffman(RLE(junc.coverage), huffmanIndex)
    else:
        s += writeCov(RLE(junc.coverage))

    return s

def readJunction(s, junc, readLenBytes, huffmanTree=None):
    fragLenBytes, s = binaryToVal(s, 1)

    junc.readLens, s = readLens(s, fragLenBytes)
    if len(junc.readLens) > 0:
        junc.lensLeft, s = readLens(s, readLenBytes)
        if len(junc.lensLeft) > 0:
            junc.lensRight, s = readLens(s, readLenBytes)

    if not huffmanTree == None:
        junc.coverage, s = readCovHuffman(s, huffmanTree)
    else:
        junc.coverage, s = readCov(s)

    return junc, s

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

def readLens(s, lenBytes):
    # Read number of lengths
    numLens, s = binaryToVal(s, 2)
    lens = dict()

    if numLens > 0:
        # Read number of bytes for each frequency
        freqBytes, s = binaryToVal(s, 1)

        for _ in range(numLens):
            l, s = binaryToVal(s, lenBytes)
            n, s = binaryToVal(s, freqBytes)
            lens[l] = n

    return lens, s

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

def readCov(s):
    # Read size of length and cov in bytes
    lenBytes, s = binaryToVal(s, 1)
    covBytes, s = binaryToVal(s, 1)

    # Read the length of the vector
    lenCov, s = binaryToVal(s, 2)

    cov = []
    for _ in range(lenCov):
        c, s = binaryToVal(s, covBytes)
        l, s = binaryToVal(s, lenBytes)
        #cov.append([l,c])
        cov += [c] * l

    return cov, s

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

def readCovHuffman(s, huffmanTree):
    # Read size of length and cov in bytes
    lenBytes, s = binaryToVal(s, 1)

    # Read the length of the vector
    lenCov, s = binaryToVal(s, 2)

    # Read lengths
    lengths = [0] * lenCov
    for i in range(lenCov):
        l, s = binaryToVal(s, lenBytes)
        lengths[i] = l

    # Read coverage 
    #l = binaryToVal(s, 4)
    #c, s = huffman.decode(huffman.bytesToBits(s, lenCov), huffmanTree, lenCov)

    c, s = huffman.decode(s, huffmanTree, lenCov)

    cov = []
    for i in range(lenCov):
        cov += [c[i]] * lengths[i]

    return cov, s


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

def readHuffmanIndex(s):
    # Read number of values
    numVals, s = binaryToVal(s, 4)

    # Read the number of bytes used to write each value
    valBytes, s = binaryToVal(s, 1)

    index = dict()
    for _ in range(numVals):
        # Read the value
        val, s = binaryToVal(s, valBytes)

        # Read the number of bits in the encoding
        numBits, s = binaryToVal(s, 1)

        # Read the bits
        numBytes = math.ceil(numBits / 8)
        bits = huffman.bytesToBits(s[:numBytes])
        s = s[numBytes:]

        index[val] = bits[:numBits]
    return index, s

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

