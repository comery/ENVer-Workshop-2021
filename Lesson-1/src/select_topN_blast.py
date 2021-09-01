#!/usr/bin/evn python3

import sys

if len(sys.argv) < 4:
    sys.exit(f"python3 {sys.argv[0]} <blastout> <topmatch|int> <tohit|int> [desc.txt|option]")


def parser_blast(file):
    pre = ""
    target_match = []
    with open(file, 'r') as fh:
        while True:
            line = fh.readline().strip()
            if not line:
                break
            tmp = line.strip().split("\t")
            if pre == "":
                target_match = [line,]
                pre = tmp[0]
            elif tmp[0] == pre:
                target_match.append(line)
            else:
                yield target_match
                target_match = [line,]
                pre = tmp[0]
        yield target_match


def getinfo(desc):
    info = {}
    with open(desc, 'r') as fh:
        for i in fh:
            if i.startswith("#"):continue
            tmp = i.strip().split()
            info[tmp[0]] = " ".join(tmp[1:])
    return info


def main():
    match_lim = int(sys.argv[2])
    hit_lim = int(sys.argv[3])
    if len(sys.argv) > 4:
        descs = getinfo(sys.argv[4])

    for arr in parser_blast(sys.argv[1]):
       # print(len(arr))
        match = 0
        hit = {}
        pre = " "
        output = []
        for line in arr:
            tmp = line.split("\t")
            querry = tmp[1]
            if len(sys.argv) > 4:
                if querry in descs:
                    info = descs[querry]
                    tmp.append(info)
                else:
                    tmp.append("-")
            if pre == " ":
                match = 1
                pre = querry
                hit[pre] = 1
                output.append("\t".join(tmp))
            elif querry == pre and hit[querry] >= hit_lim:
                continue
            elif querry == pre:
                hit[querry] += 1
                output.append("\t".join(tmp))
            elif match >= match_lim:
                break
            else:
                match += 1
                pre = querry
                hit[pre] = 1
                output.append("\t".join(tmp))

        print("\n".join(output))

if __name__ == "__main__":
    main()


