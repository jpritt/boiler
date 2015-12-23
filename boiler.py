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

    compressedName = 'compressed.bin'
    if args.output:
        compressedName = args.output

    if args.compress:
        if args.alignments is None:
            print('In order to run compression, you must specify a SAM file to compress with the --alignments flag')
            exit()

        print('Compressing to %s' % compressedName)

        compressor = compress.Compressor(args.force_xs, args.frag_len_cutoff)
        startTime = time.time()
        import cProfile
        import pstats
        import io
        pr = cProfile.Profile()
        pr.enable()
        compressor.compress(args.alignments, compressedName, args.intermediate, args.frag_len_z_cutoff, binary, debug)
        pr.disable()
        s = io.StringIO()
        sortby = 'tottime'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats(30)
        print(s.getvalue())
        endTime = time.time()

        print('Compression time: %0.3f s' % (endTime-startTime))
        print('')

    if args.decompress:
        if args.expand_to is None:
            print('In order to run decompression, you must specify a file to expand to with the --expand-to flag')
            exit()

        print('Decompressing from %s' % compressedName)

        expander = expand.Expander(args.force_xs)
        startTime = time.time()


        expander.expand(compressedName, args.expand_to, binary, debug)


        endTime = time.time()
        print('Decompression time: %0.3f s' % (endTime-startTime))
        print('')

if __name__ == '__main__':

    if '--version' in sys.argv:
        print('Boiler v' + VERSION)
        sys.exit(0)

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--compress', help="Run compression", action="store_true")
    parser.add_argument('--decompress', help="Run decompression", action="store_true")
    parser.add_argument('--alignments', type=str, help='Full path of SAM file containing aligned reads')
    parser.add_argument("--binary", help="Write in binary format", action="store_true")
    parser.add_argument("--expand-to", type=str, help="After compressing, decompress to this filename")
    parser.add_argument("--output", type=str, help="Compressed filename")
    parser.add_argument("--intermediate", type=str, help="Name of SAM file to write to after processing but before compressing")
    parser.add_argument("--debug", help="Print debug information", action="store_true")
    parser.add_argument("--force-xs", help="If we decompress a spliced read with no XS value, assign it a random one (so Cufflinks can run)", action="store_false")
    parser.add_argument("--frag-len-cutoff", type=int, help='Store any fragments longer than this in a bundle-spanning bucket')
    parser.add_argument("--frag-len-z-cutoff", type=float, help='Store any fragments above this z-score in a bundle-spanning bucket')

    args = parser.parse_args(sys.argv[1:])

    go(args)

