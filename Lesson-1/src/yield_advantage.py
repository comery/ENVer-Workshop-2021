
for n in range(1000):
    a=n
    print(a) # 相当于 return a
print("*" * 100)

# 生成器实现
def foo(num):
    print("starting...")
    while num<1000:
        num=num+1
        yield num

for n in foo(0):
    print(n)
