#!/usr/bin/env python3
import os
import gzip
from icecream import ic
def smart_open(file, opera):
    if opera == 'r':
        if os.path.exists(file) ==False:
            print("Can not open file {}".format(file))
            exit()
        else:
            if file.endswith(".gz"):
                out = gzip.open(file, 'rt')
            else:
                out = open(file, 'r')
    elif opera == 'w':
        if file.endswith(".gz"):
            out = gzip.open(file, 'wt')
        else:
            out = open(file, 'w')
    return out


with smart_open("../data/demo.fasta", 'r') as fh:
    ic(fh)

with smart_open("../data/demo_1.fq.gz", 'r') as fq:
    ic(fq)


