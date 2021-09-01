#!/usr/bin/env python3

from icecream import ic

with open("../data/demo_1.fq", 'r') as fh:
    name, seq = None, []
    for line in fh:
        line = line.rstrip()
        if line.startswith("@"):
            if name:
                ic(name)
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        ic(name)

