#! /usr/bin/env python3
import sys
import argparse
import expand

def go(args):
    expander = expand.Expander()
    if args['bundles']:
        print('Querying bundles')
        print(expander.getGeneBounds(args['compressed'], args['chrom'], args['start'], args['end']))
    if args['coverage']:
        print('Querying coverage')
        print(expander.getCoverage(args['compressed'], args['chrom'], args['start'], args['end']))
    if args['reads']:
        print('Querying reads')
        print(expander.getReads(args['compressed'], args['chrom'], args['start'], args['end']))

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--compressed', help="Path to compressed file created by Boiler", type=str, required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--bundles', help="Query bundles", action="store_true")
    group.add_argument('--coverage', help="Query coverage", action="store_true")
    group.add_argument('--reads', help="Query reads", action="store_true")
    parser.add_argument('--chrom', help="Chromosome to query", type=str, required=True)
    parser.add_argument('--start', help="Beginning of range to query", type=int)
    parser.add_argument('--end', help="End of range to query", type=int)

    args = parser.parse_args(sys.argv[1:])

    go(vars(args))