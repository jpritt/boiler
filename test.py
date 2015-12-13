#! /usr/bin/env python

import alignments
import time
import copy

def updateRLE(vector, start, end, val):
    '''
    Update the run-length encoded vector by adding val to each base in the range [start, end)
    :param vector:
    :param start:
    :param end:
    :param val:
    :return:
    '''

    length = end - start

    len_vec = len(vector)
    i = 0
    while i < len_vec:
        if start < vector[i][1]:
            break
        else:
            start -= vector[i][1]
            i += 1

    if i >= len_vec:
        return vector

    if start > 0:
        vector = vector[:i] + [[vector[i][0], start], [vector[i][0], vector[i][1]-start]] + vector[i+1:]
        i += 1
        len_vec += 1

    while i < len_vec:
        if length < vector[i][1]:
            break
        else:
            vector[i][0] += val
            length -= vector[i][1]
            i += 1

    if i < len_vec and length > 0:
        vector = vector[:i] + [[vector[i][0]+val, length], [vector[i][0], vector[i][1]-length]] + vector[i+1:]

    return vector


vec = [[0,100]]
vec = updateRLE(vec, 25, 75, 1)
print(vec)
vec = updateRLE(vec, 25, 100, 1)
print(vec)
vec = updateRLE(vec, 10, 50, 1)
print(vec)
vec = updateRLE(vec, 15, 80, 1)
print(vec)
vec = updateRLE(vec, 90, 110, 1)
print(vec)