#! /usr/bin/env python

def getCoverageFromRLE(RLE, start, end):
    ''' Return the coverage vector expanded from the RLE vector between start and end.
    '''
    coverage = []

    index = 0
    offset = RLE[index][1]
    while offset < start:
        index += 1
        offset += RLE[index][1]

    coverage += [RLE[index][0]] * (min(offset, end) - start)

    while offset < end:
        index += 1
        nextOffset = offset + RLE[index][1]

        if nextOffset <= end:
            coverage += [RLE[index][0]] * RLE[index][1]
        else:
            coverage += [RLE[index][0]] * (end - offset)

        offset = nextOffset

    if not len(coverage) == end-start:
        print 'Error!'

    print coverage