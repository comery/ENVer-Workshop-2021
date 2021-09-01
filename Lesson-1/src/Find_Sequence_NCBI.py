from Bio import Entrez,SeqIO
​
# 参数设置
Entrez.email = "example@163.com"
Entrez.tool  = "exampleScript"
​
# 查询 oct4 基因的在 Nucleotide 中的总数
hd_egquery = Entrez.egquery(term="oct4")
read_egquery = Entrez.read(hd_egquery)
total = 0
for ele in read_egquery["eGQueryResult"]:
    if ele["MenuName"] == "Nucleotide":
        total = ele["Count"]
​
# 得到查询 id 列表
hd_esearch = Entrez.esearch(db="nucleotide", term="oct4", retmax=total)
read_esearch = Entrez.read(hd_esearch)
# 这里我们只取前两个序列
ids = read_esearch["IdList"][:2]
​
# 用得到的 id 列表去下载每一条 fasta 文件，并合并，以便后续分析使用（比如进化树构建）
hd_efetch_fa = Entrez.efetch(db='nucleotide', id=ids, rettype='fasta')
read_efetch_fa = hd_efetch_fa.read()
with open("res/oct4.fasta","w") as file:
    file.write(read_efetch_fa)
# 同理你可以得到 xml 格式的序列信息
hd_efetch_xml = Entrez.efetch(db="nucleotide", id=ids, retmode="xml")
read_efetch_xml = Entrez.read(hd_efetch_xml)
print(read_efetch_xml)
hd_efetch_gb = Entrez.efetch(db="nuccore", id=ids, rettype="gb", retmode="text")
# 这里读取的是文本文件，保存为本地数据
read_efetch_gb = hd_efetch_gb.read()
with open("res/oct4.gb","w") as file:
    file.write(read_efetch_gb)
​
# 如果需要提取其中一些信息，可以按照以下步骤, 这里需要注意需要重新得到 efetch 句柄
hd_efetch_gb = Entrez.efetch(db="nuccore", id=ids, rettype="gb", retmode="text")
parse_efetch_gb = SeqIO.parse(hd_efetch_gb, "gb")
# 这里可以保存为 xls 或者 csv 格式
for ele in parse_efetch_gb:
    print(ele.name, ele.annotations['molecule_type'], ele.seq)