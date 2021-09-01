#!/usr/bin/evn python3

import sys
from icecream import ic
if len(sys.argv) < 2:
    print(f"python3 {sys.argv[0]} <blastout>")
    exit()

def addtodict2(thedict,key_a,key_b,val):
    if key_a in thedict:
        thedict[key_a].update({key_b:val})
    else:
        thedict.update({key_a:{key_b:val}})


def parser_blast(file):
    pre = ""
    target_match = {}
    with open(file, 'r') as fh:
        while True:
            line = fh.readline().strip()
            if not line:
                break
            tmp = line.strip().split("\t")
            identity, score = float(tmp[2]), float(tmp[-1])
            if pre == "":
                target_match[line]= identity
                pre = tmp[0]
            elif tmp[0] == pre:
                target_match[line]= identity
            else:
                yield target_match
                target_match = {}
                target_match[line]= identity
                pre = tmp[0]
        yield target_match



def main():

    for arr in parser_blast(sys.argv[1]):
        b = sorted(arr.items(), key=lambda x: x[1], reverse=True)
        ic(b)

if __name__ == "__main__":
    main()


