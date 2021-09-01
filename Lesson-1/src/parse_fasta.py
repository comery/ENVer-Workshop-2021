#!/usr/bin/env python3

import sys
from icecream import ic

with open("../data/demo.fasta", 'r') as fh:
    name, seq = None, []
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                ic(name, seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        ic(name, seq)

f=open('../data/demo.fasta', 'r')
seq={}
for line in f:
        if line.startswith('>'):
                name=line.replace('>','').split()[0]
                seq[name]=''
        else:
                seq[name]+=line.replace('\n','').strip()
f.close()

ic(seq)
