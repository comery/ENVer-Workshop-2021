#!/usr/bin/env python3

import sys

def parse_fasta(fh):
    name, seq = None, []
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield name, ''.join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield name, ''.join(seq)


with open("../data/demo.fasta", 'r') as fh:
    for name, seq in parse_fasta(fh):
        print(f"{name}\t{seq}")
