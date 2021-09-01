# 一个普通函数：
def foo():
    print('Starting.....')
# 调用函数，直接执行语句
g = foo()
print("*" * 100)

# 包含 yield 关键字，就变成了生成器函数
# 调用函数并不会执行语句
def foo():
    print('Starting.....')
    while True:
        res = yield 4
        print("res:", res)


# 逐行运行代码观察效果
g = foo()
print("第一次调用执行结果：")
print(next(g))
print("*" * 100)

print("第二次调用执行结果：")
print(next(g))
print("*" * 100)
