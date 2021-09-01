import time

def revcom_base(sequence):
    # make a sequence complement #
    sequence = sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    return sequence.upper()[::-1]

def revcom_dict(s):
    complement = { 'A' : 'T', 'G' : 'C', 'C' : 'G', 'T' : 'A',
                 'a' : 'T', 'g' : 'C', 'c' : 'G', 't' : 'A'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def revcom_transtable(sequence):
    # make a sequence complement #
    # replace function of string is too low!
    sequence = sequence.upper()
    transtable = str.maketrans('ATCG', 'TAGC')
    sequence = sequence.translate(transtable)
    return sequence[::-1]


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def run_time(info):
    print("[INFO]: {0} total run time: {1:.2f}".format(info,
                                                       time.time() - t) + "s")

seqs = []
with open('../data/demo.fasta') as fp:
    for name, seq in read_fasta(fp):
        seqs.append(seq)

# method1
t = time.time()
for i in seqs:
    a = revcom_base(i)
run_time("method1")
t = time.time()

for i in seqs:
    b = revcom_dict(i)
run_time("method2")
t = time.time()

for i in seqs:
    c= revcom_transtable(i)
run_time("method3")
