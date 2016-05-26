#! /usr/bin/env python
import sys
import argparse
import time
import logging
import os

VERSION = '1.0.1'

def go(args):
    if sys.version_info < (3,0):
        print('Boiler requires Python version 3 or better to run')
        exit()

    if args.command == 'compress':
        import compress

        if args.preprocess:
            if args.preprocess.lower() == 'tophat':
                print('Preprocessing TopHat alignments')
                import inferXStags
                prefix = args.alignments[:args.alignments.index('.')] + '.processed'
                inferXStags.inferTags(args.alignments, prefix + '.sam')
                os.system('samtools view -bS ' + prefix + '.sam | samtools sort - ' + prefix)
                os.system('samtools view -h -o ' + prefix + '.sam ' + prefix + '.bam')
            elif args.preprocess.lower() == 'hisat':
                print('Preprocessing HISAT alignments')
                import removeUp
                import enumeratePairs
                import inferXStags
                prefix = args.alignments[:args.alignments.index('.')] + '.processed'
                removeUp.removeUnmapped(args.alignments, 'temp1.sam')
                enumeratePairs.processHISAT('temp1.sam', 'temp2.sam')
                os.system('rm temp1.sam')
                inferXStags.inferTags('temp2.sam', prefix + '.sam')
                os.system('rm temp2.sam')
                os.system('samtools view -bS ' + prefix + '.sam | samtools sort - ' + prefix)
                os.system('samtools view -h -o ' + prefix + '.sam ' + prefix + '.bam')
            args.alignments = prefix + '.sam'

        if args.verbose:
            print('Compressing')
            start = time.time()
        compressor = compress.Compressor(args.frag_len_cutoff)
        compressor.compress(args.alignments, args.compressed, args.gtf, None, args.frag_len_z_cutoff, args.split_diff_strands, args.split_discordant)
        if args.verbose:
            end = time.time()
            print('Compression took %0.3f s' % (end-start))

    elif args.command == 'query':
        import expand

        if args.output:
            #logging.info('Opening %s to write results' % args.output)
            try:
                f = open(args.output, 'w')
            except IOError:
                #logging.info('Couldn\'t open file %s for writing. Using standard out instead.' % args.output)
                f = sys.stdout
        else:
            #logging.warning('Writing results to standard out')
            f = sys.stdout

        expander = expand.Expander()
        if args.bundles:
            bundles = expander.getGeneBounds(args.compressed, args.chrom, args.start, args.end)
            for b in bundles:
                f.write(str(b[0])+'\t'+str(b[1])+'\n')
        if args.coverage:
            cov = expander.getCoverage(args.compressed, args.chrom, args.start, args.end)
            f.write(','.join([str(c) for c in cov]) + '\n')
        if args.reads:
            aligned, unpaired, paired = expander.getReads(args.compressed, args.chrom, args.start, args.end)
            aligned.writeSAM(f, unpaired, paired, False, False, 0)

    elif args.command == 'decompress':
        import expand

        if args.verbose:
            print('Decompressing')
            start = time.time()
        expander = expand.Expander(args.force_xs)
        expander.expand(args.compressed, args.expanded)
        if args.verbose:
            end = time.time()
            print('Decompression took %0.3f s' % (end-start))

if __name__ == '__main__':

    if '--version' in sys.argv:
        print('Boiler v' + VERSION)
        sys.exit(0)

    # Check for Python version
    version = sys.version_info[0]
    if version < 3:
        print('Python 3 is required to run Boiler')
        exit()

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter, prog='Boiler')

    subparsers = parser.add_subparsers(help='Commands', dest='command')

    parser_compress = subparsers.add_parser('compress', help="Compress a SAM file")
    parser_compress.add_argument("-c", "--frag-len-cutoff", type=int, help='Store any fragments longer than this in a bundle-spanning bucket')
    parser_compress.add_argument("-z", "--frag-len-z-cutoff", type=float, help='Store any fragments above this z-score in a bundle-spanning bucket')
    parser_compress.add_argument("-s", "--split-diff-strands", action="store_true", help='Split any pairs with different XS values')
    parser_compress.add_argument("-d", "--split-discordant", action="store_true", help='Treat discordant pairs as unpaired reads')
    parser_compress.add_argument("-p", "--preprocess", type=str, help="Set to 'tophat' to preprocess TopHat alignments, 'hisat' to preprocess HISAT alignments")
    parser_compress.add_argument("-g", "--gtf", type=str, help="Path to reference GTF to improve compression accuracy (this will result in larger file size)")
    parser_compress.add_argument("-v", "--verbose", help="Print timing information", action="store_true")
    parser_compress.add_argument("alignments", type=str, help='Full path of SAM file containing aligned reads')
    parser_compress.add_argument("compressed", type=str, nargs='?', default='compressed.bin', help="Compressed filename. Default: compressed.bin")

    parser_query = subparsers.add_parser('query', help="Query compressed file")
    group = parser_query.add_mutually_exclusive_group()
    group.add_argument('-b', '--bundles', help="Query bundles", action="store_true")
    group.add_argument('-c', '--coverage', help="Query coverage", action="store_true")
    group.add_argument('-r', '--reads', help="Query reads", action="store_true")
    parser_query.add_argument('--chrom', help="Chromosome to query", type=str, required=True)
    parser_query.add_argument('--start', help="Beginning of range to query", type=int)
    parser_query.add_argument('--end', help="End of range to query", type=int)
    parser_query.add_argument('compressed', help="Path to compressed file created by Boiler", type=str)
    parser_query.add_argument('output', nargs='?', default=None, help="File to write result to. Default: Standard out")

    parser_decompress = subparsers.add_parser('decompress', help="Decompress to a SAM file")
    parser_decompress.add_argument("-f", "--force-xs", help="If we decompress a spliced read with no XS value, assign it a random one (so Cufflinks can run)", action="store_true")
    parser_decompress.add_argument("-v", "--verbose", help="Print timing information", action="store_true")
    parser_decompress.add_argument("compressed", type=str, help="Compressed filename")
    parser_decompress.add_argument("expanded", type=str, nargs='?', default='expanded.sam', help="Write decompressed SAM to this filename. Default: expanded.sam")

    args = parser.parse_args(sys.argv[1:])

    go(args)

