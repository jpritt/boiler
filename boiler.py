#! /usr/bin/env python3
import sys
from multiprocessing import Process
import os
import argparse
#import genome
import time
import expand
import compress


VERSION = '0.0.1'

def go(args):
    binary = False
    if args.binary:
        binary = True

    huffman = False
    if args.huffman:
        huffman = True

    compressedName = 'compressed/compressed_binary.txt'
    if args.output is not None:
        compressedName = args.output

    compressor = compress.Compressor()
    #print('Initial')
    #objgraph.show_growth(limit=3)
    startTime = time.time()
    compressor.compress(args.alignments, compressedName, args.intermediate, binary, huffman)
    endTime = time.time()
    compressor = None

    #print('Final')
    #objgraph.show_growth()
    print('Compression time: %0.3f s' % (endTime-startTime))
    #print('Alignments size: %f\n' % (asizeof.asizeof(compressor.aligned)/1000000))

    if args.expand_to is not None:
        expander = expand.Expander()
        startTime = time.time()
        expander.expand(compressedName, args.expand_to, binary, huffman)
        endTime = time.time()
        print('Decompression time: %0.3f s' % (endTime-startTime))
        #print('Alignments size: %f\n' % (asizeof.asizeof(expander.aligned)/1000000))

if __name__ == '__main__':

    if '--version' in sys.argv:
        print('Boiler v' + VERSION)
        sys.exit(0)

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--alignments', type=str, required=True, 
        help='Full path of SAM file containing aligned reads')
    parser.add_argument("--binary", help="Write in binary format",
        action="store_true")
    parser.add_argument("--huffman", help="Use huffman coding to compress coverage vectors",
        action="store_true")
    parser.add_argument("--expand-to", type=str, help="After compressing, decompress to this filename")
    parser.add_argument("--output", type=str, help="Compressed filename")
    parser.add_argument("--intermediate", type=str, help="Name of SAM file to write to after processing but before compressing")

    args = parser.parse_args(sys.argv[1:])

    go(args)

