#!usr/bin/python
########################################################################
###this script used for modify surface RNA name
import os
import sys
a = sys.argv[1]

with open(a) as handle:
    i1 = handle.readline().strip()
    print(i1)
    for line in handle:
        line = line.strip().split('\t')
        name = line[0].split('|')[-4]
        print('\t'.join([name]+line[1:]))
