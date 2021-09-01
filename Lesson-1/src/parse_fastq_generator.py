#!/usr/bin/env python3

import sys
from icecream import ic

def parse_pe_fastq(fq1, fq2):
    while True:
        name1 = fq1.readline().split(" ")[0].replace("/1", "")
        name2 = fq2.readline().split(" ")[0].replace("/2", "")
        if not name1 or not name2:
            break
        read1, nothing1, qual1 = fq1.readline()[:-1], fq1.readline(), fq1.readline()[:-1]
        read2, nothing2, qual2 = fq2.readline()[:-1], fq2.readline(), fq2.readline()[:-1]
        assert name1 == name2, 'fastq1, fastq2 is not paired-end'
        yield name1, read1, read2, qual1, qual2


if __name__ == '__main__':
    fq1 = open("../data/demo_1.fq", 'r')
    fq2 = open("../data/demo_2.fq", 'r')
    for name, read1, read2, qual1, qual2 in parse_pe_fastq(fq1, fq2):
        ic(name, read1, read2, qual1, qual2)



