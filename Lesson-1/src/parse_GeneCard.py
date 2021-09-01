#https://www.genecards.org/cgi-bin/carddisp.pl?gene=TRIP6
import re
import os
import sys
import time
from bs4 import BeautifulSoup

if len(sys.argv) < 3:
    sys.exit("usage: python3 {} <gene.list> <inputdir>".format(sys.argv[0]))


def format_str(string):
    print(string.replace("\n", "").strip())

def makeup_md(info):
    for i in info.split("\n"):
        if len(i.strip()) > 0:
            print(i.strip())

def summary_info(soup):
    all_info = []
    for div in soup.find_all("div", class_="gc-subsection"):
        if div.find('h3') and "Summary" in div.find('h3').get_text():
            for i in div.descendants:
                if i.name == 'a' or i.name == 'p' or i.name == 'div':
                    if len(i.get_text()) > 0:
                        all_info.append(i.get_text().strip())
    tmp =  "\n".join(all_info)
    tmp = tmp.replace("GeneCards Summary for", "\nGeneCards Summary for")
    tmp = tmp.replace("UniProtKB/Swiss-Prot Summary for", "\nUniProtKB/Swiss-Prot Summary for")
    return tmp


def phenotype(soup):
    all_info = []
    for div in soup.find_all("div", class_="gc-subsection", id="function-phenotypes"):
        if div.find('h3') and "Phenotypes" in div.find('h3').get_text():
            for i in div.descendants:
                if i.name == 'a' or i.name == 'p' or i.name == 'div':
                    if len(i.get_text()) > 0:
                        all_info.append(i.get_text().strip())

    return "\n".join(all_info)

def translation(input_content, src, dest):
    # 导入模块
    from googletrans import Translator

    # 实例化翻译器，由于模块默认的服务url在国内无法使用，所以我们修改成国内可用的google翻译服务地址
    translator = Translator(service_urls=["translate.google.cn"])

    # 调用翻译函数，指定原语言的代码(en)，和目标语言的代码(zh-CN)
    result = translator.translate(input_content, src=src, dest=dest)

    # 原语言代码
    #print(result.src)
    # 目标语言代码
    #print(result.dest)
    # 要翻译的内容
    #print(result.origin)
    # 翻译后的内容
    #print(result.text)
    # 翻译后内容的发音
    #print(result.pronunciation)
    return result.origin


def main():
    inputdir = sys.argv[2]
    if os.path.exists(inputdir) == False:
        sys.exit("can not find {}".format(inputdir))

    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            gene = i.strip()
            print("#### " + gene )
            infile = inputdir + "/" + gene + ".html"
            if os.path.exists(infile) == False:
                print("Can not find html of {} in {}".format(gene, inputdir))
            else:
                fh = open(infile, 'r')
                html = fh.read()
                fh.close()
                soup = BeautifulSoup(html, features="html5lib")
                # summary info
                print("\n### Summary")
                summaries = summary_info(soup)
                makeup_md(summaries)
                #translation(summaries, src='en', dest='zh-CN')
                #tmp_tab = [gene, summaries, translation]
                #print("\t".join(tmp_tab))
                print(summaries)

                # phenotype
                #print("\n### Phenotypes")
                #phenotypes = phenotype(soup)
                #makeup_md(phenotypes)

            print("\n")




if __name__ == '__main__':
    main()

