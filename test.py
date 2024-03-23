import numpy as np
from scipy.integrate import laplace_transform
 
# 定义原始函数
def f(t):
    return t**2 + 3*np.exp(-0.5*t)
 
# 对原始函数进行拉普拉斯变换
result = laplace_transform(f, t, s)
print("拉普拉斯变换结果：", result)