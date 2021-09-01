#https://www.genecards.org/cgi-bin/carddisp.pl?gene=TRIP6
import re
import os
import sys
import time
import urllib
import requests
from urllib import request
if len(sys.argv) < 3:
    sys.exit("usage: python3 {} <gene.list> <outdir>".format(sys.argv[0]))


def getContentGeneCard(gene):
    api = "https://www.genecards.org/cgi-bin/carddisp.pl?"
    item = "gene={}&keywords={}".format(gene, gene)
    url = api + item

    headers = {
        # paste your specific information here, then you can run this script
    }

    params = (
        ('gene', gene),
        ('keywords', gene),
    )

    response = requests.get('https://www.genecards.org/cgi-bin/carddisp.pl', headers=headers, params=params)

    #NB. Original query string below. It seems impossible to parse and
    #reproduce query strings 100% accurately so the one below is given
    #in case the reproduced version is not "correct".
    # response = requests.get('https://www.genecards.org/cgi-bin/carddisp.pl?gene=RTN4R&keywords=RTN4R', headers=headers)

    print(url)
    content = response.text
    return content

def write_html(html, outfile):
    with open(outfile, 'w') as fh:
        fh.write(html)

def main():
    outdir = sys.argv[2]
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            gene = i.strip()
            outfile = outdir + "/" + gene + ".html"
            if os.path.exists(outfile) == False or os.path.getsize(outfile) == 0:
                html = getContentGeneCard(gene)
                write_html(html, outfile)
                time.sleep(4)
            else:
                print("{} html has existed".format(gene))

if __name__ == '__main__':
    main()

