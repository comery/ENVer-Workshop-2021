import sys
import re
import jieba
from wordcloud import WordCloud
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print("Usage: python3 " + sys.argv[0] + " *.txt\n")
    exit(0)

def count(string):
    ok = True
    for i in string.split("\n"):
        if i.count(",") >5:
            ok = False
            break
        elif i.count("(") > 3:
            ok = False
            break
    return ok

tmp = ""
text = ""
with open(sys.argv[1],'r') as fh:
    for i in fh.readlines():
        i = i.strip()
        if len(i) == 0:
            if len(tmp.split("\n")) >4 and tmp[:6] != "Author" and count(tmp) == True:
                text += tmp
                tmp = ""
            else:
                tmp = ""
        else:
            tmp += i + "\n"

# width,height,margin可以设置图片属性

#wordcloud = WordCloud(font_path = r'D:\Fonts\simkai.ttf').generate(f)
# 你可以通过font_path参数来设置字体集
#background_color参数为设置背景颜色,默认颜色为黑色

wordlist_after_jieba = jieba.cut(text, cut_all = True)
wl_space_split = " ".join(wordlist_after_jieba)
# this size is suit for PPT page.
wordcloud = WordCloud(background_color="Black",width=3385, height=1905, margin=2).generate(wl_space_split)

plt.imshow(wordcloud)
plt.axis("off")
plt.show()

wordcloud.to_file('test1.png')
