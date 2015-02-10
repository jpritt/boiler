#! /usr/bin/env python
import sys
import argparse
#import genome
import time
import expand
import compress


if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--alignments', type=str, required=True, 
        help='Full path of SAM file containing aligned reads')
    parser.add_argument("--binary", help="Write in binary format",
        action="store_true")
    parser.add_argument('--protocol', type=int, required=False, help='Pickle protocol to use. Default = 2')
    
    args = parser.parse_args(sys.argv[1:])

    binary = False
    if args.binary:
        binary = True

    if args.protocol == None:
        protocol = 2
    else:
        protocol = args.protocol

    compressedName = 'compressed/compressed.txt'
    expandedName = 'expanded.sam'

    compressor = compress.Compressor()
    startTime = time.time()
    compressor.compress(args.alignments, compressedName, binary, protocol)
    endTime = time.time()
    print 'Compression time: %0.3f s' % (endTime-startTime)

    expander = expand.Expander()
    startTime = time.time()
    expander.expand(compressedName, expandedName, binary)
    endTime = time.time()
    print 'Decompression time: %0.3f s' % (endTime-startTime)
