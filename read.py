#!usr/bin/python

import os
f = os.listdir()
out = {}
f = [i for i in f if i.endswith('V2.xls')]
for i in f:
    with open(i) as handle:
        name = i.split('_')[0]
        for line in handle:
            line = line.strip().split('\t')
            if len(line) >2:
                cell = line[1].strip('"')
                number = line[2].strip('"')
                print('\t'.join([name,number,cell]))

