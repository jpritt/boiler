import math
from struct import *

class HuffmanNode:
    def __init__(self, freq, val=None, child1=None, child2=None):
        self.freq = freq
        self.val = val
        self.child1 = child1
        self.child2 = child2
        self.parent = None

    def getCode(self):
        if self.parent == None:
            return []
        else:
            if self.parent.child1 == self:
                return self.parent.getCode() + [1]
            else:
                return self.parent.getCode() + [0]

    def readVal(self, code):
        if self.child1 == None:
            return self.val, code
        elif code[0] == 1:
            return self.child1.readVal(code[1:])
        else:
            return self.child2.readVal(code[1:])

class LinkedNode:
    def __init__(self, node):
        self.node = node
        self.prev = None
        self.next = None

def createHuffmanTree(data):
    ''' Data is a list of tuples containing values and their frequencies
    '''

    print('Creating huffman nodes')

    data = sorted(data, key=lambda d: d[1])

    # Create list of leaf nodes, one for each value
    first = None
    leaves = []
    for d in data:
        n = HuffmanNode(d[1], d[0])
        leaves.append(n)
        node = LinkedNode(n)
        if first == None:
            first = node
            last = node
        else:
            last.next = node
            last = node

    print('Building tree')
    # Iteratively connect the 2 least frequent nodes until all are connected
    while not first.next == None:
        n1 = first
        n2 = n1.next
        newNode = LinkedNode(HuffmanNode(n1.node.freq+n2.node.freq, None, n1.node, n2.node))
        n1.node.parent = newNode.node
        n2.node.parent = newNode.node

        first = n2.next
        if first == None:
            return newNode.node, leaves
        elif newNode.node.freq < first.node.freq:
            first.prev = newNode
            newNode.next = first
            first = newNode
        else:
            n = first

            while not n.next == None and newNode.node.freq > n.next.node.freq:
                n = n.next
            if not n.next == None:
                n.next.prev = newNode
            newNode.next = n.next
            newNode.prev = n
            n.next = newNode

def createHuffmanIndex(data):
    tree, leaves = createHuffmanTree(data)

    # Find the encoding for each value
    index = dict()
    codeLens = dict()
    for n in leaves:
        index[n.val] = n.getCode()
        length = len(n.getCode())
        if not length in codeLens:
            codeLens[length] = 1
        else:
            codeLens[length] += 1
    print('Encoding lengths')
    print(codeLens)
    #exit()
    return index

def createHuffmanTreeFromIndex(index):
    # First create an index mapping codes to values
    #print('Reversing index')
    reverseIndex = dict()
    for k,v in index.items():
        reverseIndex[''.join([str(i) for i in v])] = k

    #print('Constructing tree')
    # Since we don't care about frequencies, use freq to store the path to each node
    tree = HuffmanNode('')
    activeNodes = [tree]
    leaves = []
    i = 0
    while i < len(activeNodes):
        n = activeNodes[i]

        n.child1 = HuffmanNode(n.freq + '1')
        n.child1.parent = n
        if n.child1.freq in reverseIndex:
            n.child1.val = reverseIndex[n.child1.freq]
            leaves.append(n.child1)
        else:
            activeNodes.append(n.child1)

        n.child2 = HuffmanNode(n.freq + '0')
        n.child2.parent = n
        if n.child2.freq in reverseIndex:
            n.child2.val = reverseIndex[n.child2.freq]
            leaves.append(n.child2)
        else:
            activeNodes.append(n.child2)

        i += 1

    #print('Done')
    return tree, leaves

def encode(vals, index):
    ''' Return a list of bits encoding all the values in the list
    '''
    encoding = []
    for v in vals:
        encoding += index[v]
    return encoding

'''
def decode(encoding, tree, numVals):
    vals = []
    for _ in range(numVals):
        v, encoding = tree.readVal(encoding)
        vals.append(v)
    return vals
'''

def decode(s, tree, numVals):
    #print('Decoding')
    bitBuffer = bytesToBits(s[:4])
    i = 4
    vals = []
    for _ in range(numVals):
        # Keep the next 4 bytes as a buffer
        while i < len(s) and len(bitBuffer) < 32:
            bitBuffer += bytesToBits([s[i]])
            i += 1

        v, bitBuffer = tree.readVal(bitBuffer)
        vals.append(v)
    replace = math.floor(len(bitBuffer) / 8)
    i -= replace
    #print('Done')
    return vals, s[i:]

def bitsToBytes(bits):
    ''' Return a byte string containing the given list of bits, as well as the number of padding bits at the end of the last byte
    '''

    s = b''
    while len(bits) > 8:
        val = 0
        for i in range(8):
            val = (val << 1) + bits[i]
        s += pack('B', val)

        bits = bits[8:]

    # Write the last byte
    val = 0
    padding = 8 - len(bits)
    for i in range(len(bits)):
        val = (val << 1) + bits[i]
    val = val << padding
    s += pack('B', val)
    return s

def bytesToBits(bs):
    ''' Return a list of bits from the given byte string
    '''

    bits = [0] * (8*len(bs))
    for b in range(len(bs)):
        byte = bs[b]

        for i in range(8):
            bits[8*b + 7 - i] = (byte >> i) & 1

    return bits

'''
data = [('a1', 0.05), ('a2', 0.20), ('a3', 0.30), ('a4', 0.45)]
index = createHuffmanIndex(data)
c = bitsToBytes(encode(['a2', 'a4', 'a4', 'a2', 'a1', 'a3', 'a3', 'a1'], index))
c += b'hello'
tree, n = createHuffmanTreeFromIndex(index)
v, s = decode(c, tree, 8)
print(v)
print(s)
'''
