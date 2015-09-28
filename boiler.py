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

    debug = False
    if args.debug:
        debug = True

    compressedName = 'compressed/compressed_binary.txt'
    if args.output is not None:
        compressedName = args.output

    #compressor = compress.Compressor()
    #startTime = time.time()
    #compressor.compress(args.alignments, compressedName, args.intermediate, binary, debug)
    #endTime = time.time()

    #print('Compression time: %0.3f s' % (endTime-startTime))

    if args.expand_to is not None:
        expander = expand.Expander()
        startTime = time.time()

        import cProfile
        import pstats
        import io
        pr = cProfile.Profile()
        pr.enable()

        expander.expand(compressedName, args.expand_to, binary, debug)

        pr.disable()
        s = io.StringIO()
        sortby = 'tottime'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats(30)
        print(s.getvalue())

        endTime = time.time()
        print('Decompression time: %0.3f s' % (endTime-startTime))

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
    parser.add_argument("--expand-to", type=str, help="After compressing, decompress to this filename")
    parser.add_argument("--output", type=str, help="Compressed filename")
    parser.add_argument("--intermediate", type=str, help="Name of SAM file to write to after processing but before compressing")
    parser.add_argument("--debug", help="Print debug information",
        action="store_true")

    args = parser.parse_args(sys.argv[1:])

    go(args)

