# =====查看物种谱系=====
from Bio import Entrez
​
# 参数设置
Entrez.email = "example@163.com"
Entrez.tool  = "exampleScript"
​
# 在 Taxonomy 库中搜索 Homo sapiens
hd_esearch = Entrez.esearch(db="Taxonomy", term="Homo sapiens")
read_esearch = Entrez.read(hd_esearch)
id = read_esearch["IdList"][0]
hd_eftech = Entrez.efetch(db="Taxonomy", id=id, retmode="xml")
read_eftech = Entrez.read(hd_eftech)
print(read_eftech[0]["Lineage"])