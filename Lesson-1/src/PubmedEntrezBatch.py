from Bio import Entrez
# 参数设置
Entrez.email = "example@163.com"
Entrez.tool  = "exampleScript"
​
# 搜索
hd_esearch = Entrez.esearch(db="pubmed", term="Sus scrofa", reldate=365, ptyp="Review", usehistory="y")
read_esearch = Entrez.read(hd_esearch)
total = int(read_esearch["Count"])
webenv = read_esearch["WebEnv"]
query_key = read_esearch["QueryKey"]
​
# 这里演示设定total为 10
total = 10
step = 5
print("Result items: ", total)
with open("res/recent_review_sus_scrofa.txt", "w") as file:
    for start in range(0, total, step):
        print("Download record %i to %i" % (start + 1, int(start+step)))
        hd_efetch = Entrez.efetch(db="pubmed", retstart=start, retmax=step, webenv=webenv, query_key=query_key, rettype="medline", retmode="text")
        file.write(hd_efetch.read())