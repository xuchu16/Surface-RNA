#!usr/bin/python 

import os
dirs = os.getcwd()
f = os.listdir()
out = {}
sample = []
for i in f:
    d = dirs+'/'+i
    if os.path.isdir(d):
        sample.append(d.split('/')[-1])
        exp = d+'/'+'abundance.tsv'
        with open(exp) as handle:
            for line in handle:
                line = line.strip().split('\t')        
                name = line[0]
                TPM = line[-1]
                if name not in out:
                    out[name] = []
                    out[name].append(TPM)
                else:
                    out[name].append(TPM)
print('\t'.join(['Gene Name']+sample))
for i in out:
    print('\t'.join([i]+out[i]))
