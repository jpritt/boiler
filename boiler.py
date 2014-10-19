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
    
    args = parser.parse_args(sys.argv[1:])


    compressedName = 'compressed/compressed.txt'
    expandedName = 'expanded.sam'

    compressor = compress.Compressor()
    compressor.compress(args.alignments, compressedName)

    expander = expand.Expander()
    expander.expand(compressedName, expandedName)
