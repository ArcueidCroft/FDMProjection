# * 用隐式差分法求解方程
# * Write by ArcueidCroft@SCU.math
# * Email: zqy_edu@163.com
# * Last edit: 2024.3.22 22:55

using LinearAlgebra
using SparseArrays
import Plots as plt

force(x, t) = 1

function Mainf()
    # * 初始化参数
    nx::Int = 100
    nt::Int = 100
    bwidth::Int = 2 # 边界宽度
    ht = 0.01
    hx = 0.01
    u = zeros(nt, nx)

    # * 边界控制和系数计算
    row_Q = nx-2*bwidth
    p = ht^2/(2*hx^4)
    
    # TODO 边界条件

    # * 计算u
    for k = 2:nt-1
        # * 组装Q矩阵
        # ** 首先组装非边界
        # ** 对bwidth单独组装
        Q = sparse(zeros(row_Q, row_Q))
        for row = 3:row_Q-2
            Q[row, row+2] = -p
            Q[row, row+1] = 4*p
            Q[row, row]   = (1+ht)-6*p
            Q[row, row-1] = 4*p
            Q[row, row-2] = -p
        end
        Q[1:2, 1:5] = Q[3:4, 3:3+5-1]
        Q[row_Q-1:row_Q, row_Q-1-2:end] = Q[3:4, 1:4]
        
        # * 组装f向量
        fk = zeros(row_Q)
        for row = 1:row
            fk[row] = ht^2*force(row*nx, k*nt)
        end

        # * 组装右端uk
        # ** 首先组装Qt矩阵，Qt矩阵组装思路和Q相同
        # ** 再生成uk
        Qt = sparse(zeros(row_Q, row_Q))
        for row = 3:row_Q
            Qt[row, row+2] = p
            Qt[row, row+1] = -4p
            Qt[row, row]   = (2+ht)+6*p
            Qt[row, row-1] = -4p
            Qt[row, row-2] = p
        end
        Qt[1:2, 1:5] = Qt[3:4, 3:3+5-1]
        Qt[row_Q-1:row_Q, row_Q-1-2:end] = Qt[3:4, 1:4]
        uk = Qt*u[nt, 3:end-2]

        # * 组装ukk:=u^{k-1}_{j}
        ukk = u[nt-1, 3:end-2]

        # * 求解生成u^{k+1}
        u[nt+1, :] = Q\(fk+uk+ukk)
    end

    plt.heatmap(u)
end