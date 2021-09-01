#!/usr/bin/env python3
import sys
import os
import re
from icecream import ic
usage = """
this script is a liftover like tools, to convert genome position
from one genome assembly to another genome assembly based on mummer
alignment, which looks like:

============================================================
-- Alignments between Scaffold3 and scaffold1

-- BEGIN alignment [ +1 1 - 198330 | +1 1 - 198282 ]

                  10|       20|       30|       40|
1          ggttttctaaacaaaattcgtcttctaattttcccagaaacacggctat
1          ggtttt.taaacaaaattcgtcttctaattttcccagaaacacggctat
                 ^
you just need to give the *1delta and *.1coords files, script
will automatically generate the above files.
"""

if len(sys.argv) < 7:
    print(usage)
    sys.exit(f"python3 {sys.argv[0]} *.bed *.1coords *.1delta tmpDir out.bed unMap.bed")

def testcmd(command):
    import subprocess
    ret = subprocess.run(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8",timeout=1)
    if ret.returncode == 0:
        return True
    else:
        return False


class Alignment(object):
    """
    a customed object to describe the one-to-one alignment in mummer
    alignment looks like:
============================================================
-- Alignments between Scaffold3 and scaffold1

-- BEGIN alignment [ +1 1 - 198330 | +1 1 - 198282 ]

                  10|       20|       30|       40|
1          ggttttctaaacaaaattcgtcttctaattttcccagaaacacggctat
1          ggtttt.taaacaaaattcgtcttctaattttcccagaaacacggctat
    """
    def __init__(self, fh):
        self.blocks_ref = {}
        self.blocks_qry = {}
        self.block_order = 0
        self.block_infos = {}
        self.strand = {}
        block_order = 0
        switch = 0
        for i in fh:
            i = i.strip()
            if len(i) == 0: # skipping the blank lines
                continue
            if i.startswith("-- Alignments between"): # this is a Alignment object start
                alignment_header = i.split()
                self.ref = alignment_header[3]
                self.qry = alignment_header[5]
            elif i.startswith("-- BEGIN alignment"): # this is a block start
                block_order += 1 # use a number to index blocks
                switch = 0
                self.blocks_ref[block_order] = [] # the value of dict is a list, which will contain [ref_start, ref_end, qry_start, qry_end]
                self.blocks_qry[block_order] = []
                """ info
                [ +1 198454 - 740766 | +1 198284 - 740300 ]
                """
                info = i.split()
                # ref_start, ref_end, qry_start, qry_end
                self.block_infos[block_order] = (int(info[5]), int(info[7]), int(info[10]), int(info[12]))
                if info[9] == "+1": # strand is +
                    self.strand[block_order] = "+"
                else:
                    self.strand[block_order] = "-"
            elif re.match(r'[0-9]+\s+[natcg]+', i) and switch == 0:
                pos_and_seq = i.split()
                pos = int(pos_and_seq[0])
                seq = pos_and_seq[1]
                self.blocks_ref[block_order].append((pos, seq))
                switch = 1 # change the reading model, the next line is pos + seq for qry
            elif re.match(r'[0-9]+\s+[natcg]+', i) and switch == 1:
                pos_and_seq = i.split()
                pos = int(pos_and_seq[0])
                seq = pos_and_seq[1]
                self.blocks_qry[block_order].append((pos, seq))
                switch = 0 # return back


    def get_block_order(self):
        return sorted(self.block_infos.keys())

    def within_block(self, block_order, pos):
        """
        judge whether a position in this block
        """
        ref_s, ref_e, qry_s, qry_e = self.block_infos[block_order]
        if pos >= ref_s and pos <= ref_e:
            return True
        elif pos >= ref_e and pos <= ref_s:
            return True
        else:
            return False

    def locate_line(self, block_order, pos):
        """
        find the pos in which line. return the index(line number, 0-based)
        """
        lines = self.blocks_ref[block_order]
        for index, l in enumerate(lines):
            tmp_start, tmp_seq = l
            seq_len = len(tmp_seq)
            indel = tmp_seq.count(".")
            tmp_end = seq_len + tmp_start - indel - 1
            if pos >=tmp_start and pos <= tmp_end:
                return index

    def get_qry_base(self, block_order, pos):
        index = self.locate_line(block_order, pos)
        ic(index)
        target_line_ref = self.blocks_ref[block_order][index]
        target_line_qry = self.blocks_qry[block_order][index]
        strand = self.strand[block_order]
        ref_start, ref_seq = target_line_ref
        qry_start, qry_seq = target_line_qry
        if len(ref_seq) != len(qry_seq):
            print("bad line with different length")
            exit()
        r_skip = 0
        q_skip = 0
        for i in range(len(ref_seq)):
            rb = ref_seq[i]
            qb = qry_seq[i]

            if rb == ".":
                # rb is deletion, qb is 'a, t, c, g n';
                r_skip += 1

            elif qb == ".":
                # rb is 'a, t, c, g n'; qb is deletion
                q_skip += 1

            r_real_pos = ref_start + i - r_skip
            if strand == "+":
                q_real_pos = qry_start + i - q_skip
            else:
                q_real_pos = qry_start - i + q_skip
            if r_real_pos == pos:
                return q_real_pos
                break

def addtwodimdict(thedict, key_a, key_b, val):
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a:{key_b: val}})

def inDict(thedict, key_a, key_b):
    if key_a in thedict:
        if key_b in thedict[key_a]:
            return True
        else:
            return False
    else:
        return False

def main():
    Target_BED = {}

    if testcmd("show-aligns"):
        show_aligns = '/hwfssz5/ST_DIVERSITY/PUB/USER/yangchentao/software/All_kinds_align/mummer4/bin/show-aligns'
    else:
        show_aligns = 'show-aligns'

    # reading the bed of ref, which you want to convert to qry
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            if i.startswith("#"):continue
            tmp = i.strip().split()
            tmp[1], tmp[2] = int(tmp[1]), int(tmp[2])
            if tmp[0] not in Target_BED:
                Target_BED[tmp[0]] = [(tmp[1], tmp[2]),]
            else:
                Target_BED[tmp[0]].append((tmp[1], tmp[2]))
    # reading the matched pairs of ref and qry
    Target_PAIR = set()
    with open(sys.argv[2], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            if tmp[11] in Target_BED:
                Target_PAIR.add((tmp[11], tmp[12]))
    tmpDir = sys.argv[4]
    if os.path.exists(tmpDir) == False:
        os.mkdir(tmpDir)

    QRY_POS = {}
    # scanning all pairwise alignment, a is ref_scaf, b is qryry scaf
    for a, b in Target_PAIR:
        outfile = f"{tmpDir}/{a}_{b}.alns"
        # to generate mummer alignment files for each pair
        if os.path.exists(outfile) == False or os.path.getsize(outfile) == 0:
            cmd = f"{show_aligns} {sys.argv[3]} {a} {b} >{tmpDir}/{a}_{b}.alns"
            os.system(cmd)
        else:
            1
            #print(f"{outfile} exists")

        # parser alignment
        with open(outfile, 'r') as fh:
            alignment = Alignment(fh) # instantiation of object Alignment
            for order in alignment.get_block_order(): # loop the alignment blocks
                for p1, p2 in Target_BED[a]:
                    if p1 == p2: # for snp, speed it up
                        if alignment.within_block(order, p1): # p1 in this block, then find the coordinate site in qry
                            try:
                                p1_qpos = alignment.get_qry_base(order, p1)
                                addtwodimdict(QRY_POS, a, p1, (b,p1_qpos))
                                addtwodimdict(QRY_POS, a, p2, (b,p1_qpos))
                            except TypeError:
                                ic(p1, p2)
                    else:
                        if alignment.within_block(order, p1):
                            p1_qpos = alignment.get_qry_base(order, p1)
                            addtwodimdict(QRY_POS, a, p1, (b,p1_qpos))
                        if alignment.within_block(order, p2):
                            p2_qpos = alignment.get_qry_base(order, p2)
                            addtwodimdict(QRY_POS, a, p2, (b,p2_qpos))

    # output coordiate qryry position
    out = open(sys.argv[-2], 'w')
    unmap = open(sys.argv[-1], 'w')
    for ref_scaf in Target_BED:
        for p1, p2 in Target_BED[ref_scaf]:
            if inDict(QRY_POS, ref_scaf, p1):
                q1_scaf, q1_pos = QRY_POS[ref_scaf][p1]
            else:
                q1_scaf, q1_pos = ('-', '-')

            if inDict(QRY_POS, ref_scaf, p2):
                q2_scaf, q2_pos = QRY_POS[ref_scaf][p2]
            else:
                q2_scaf, q2_pos = ('-', '-')

            if q1_scaf == q2_scaf and q1_pos != '-' and q2_pos != '-':
                print(f"{ref_scaf}\t{p1}\t{p2}\t{q1_scaf}\t{q1_pos}\t{q2_pos}", file=out)
            else:
                print(f"{ref_scaf}\t{p1}\t{p2}\t-\t{q1_scaf}:{q1_pos}\t{q2_scaf}:{q2_pos}", file=unmap)

    out.close()
    unmap.close()

if __name__ == '__main__':
        main()
