#! /usr/bin/env python3

import objgraph
import read

objgraph.show_growth(limit=4)
r1 = read.Read('2L', [100,102])
r2 = read.Read('2R', [10,12])
l = [r1, r2]
print('')
objgraph.show_growth()

r = l[1]
del[l[1]]
print(l)
print(r.exons)
print('')
objgraph.show_growth()

