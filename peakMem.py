#! /usr/bin/env python

import fileinput

max_mem = 0.0
i = 0
mem_id = 7
for line in fileinput.input():
    curr_id = 0
    if i == 0:
        i += 1
    else:
        row = line.rstrip().split(' ')
        for col in row:
            if len(col) > 0:
                if curr_id == mem_id:
                    mem = float(col)
                    if mem > max_mem:
                        max_mem = mem
                    break
                else:
                    curr_id += 1
print('Peak memory:\t%d' % max_mem)
