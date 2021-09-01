import sys
from Bio import Entrez
from Bio import Medline


if len(sys.argv) < 3:
    sys.exit(f"python3 {sys.argv[0]} term output")

def main():
    item = sys.argv[1]
    outfile = sy.argv[2]
    # 参数设置
    Entrez.email = "yangchentao@genomics.cn"
    Entrez.tool  = "exampleScript"

    # 用 esearch 在 pubmed 库中搜索关键字为 "mouse" 的文章
    # RetMax 这个参数为每次返回的最大个数，因此如果把Count的值赋给RetMax就会获取全部的mouse的文章，这里为实例设置为100
    hd_esearch = Entrez.esearch(db="pubmed", term=item, RetMax="100")
    read_esearch = Entrez.read(hd_esearch)
    idlist = read_esearch["IdList"]
    print ("Total: ", read_esearch["Count"])
    # 用 efetch下载
    hd_efetch = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text", )
    # 用 Medline 来解析
    parse_medline = Medline.parse(hd_efetch)
    with open(outfile, "w") as out:
        out.write("title\tauthors\tsource\tPubMed\n")
        for i, ele in enumerate(list(parse_medline)):
            title = ele['TI']
            author = ele['AU']
            source = ele['SO']
            pmid = ele['PMID']
            abstract = ele['AB']
            line = "\t".join([title, pmid, source, abstract, author])
            print(line, file=out)


if __name__ == '__main__':
    main()
