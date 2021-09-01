from Bio import Entrez

# 参数设置
Entrez.email = "example@163.com"
Entrez.tool  = "exampleScript"

hd_search = Entrez.esearch(db="nucleotide", term="oct4", usehistory="y")
read_search = Entrez.read(hd_search)
webenv = read_search["WebEnv"]
query_key = read_search["QueryKey"]

# 使用历史记录特性来进行搜索。
# Entrez 将会提前进行缓冲，提高查询效率
step = 5
total = 10
with open("res/res_env_oct4.fasta", "w") as res_file:
    for start in range(0, total, step):
        end = min(total, start+step)
        print("Download record %i to %i" % (start+1, end))
        hd_fetch = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", retstart=start, retmax=step, webenv=webenv, query_key=query_key)
        records = hd_fetch.read()
        res_file.write(records)