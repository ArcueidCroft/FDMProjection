# * 用隐式差分法求解方程
# * Write by ArcueidCroft@SCU.math
# * Email: zqy_edu@163.com
# * Last edit: 2024.3.22 22:55

using LinearAlgebra
using SparseArrays

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
    row_Q = nx-bwidth
    p = ht^2/(2*hx^4)
    
    # TODO 边界条件

    # * 计算u
    for k = 2:nt-1
        # * 组装Q矩阵
        Q = sparse(zeros(row_Q, row_Q))
        for row = 3:row_Q
            Q[row, row+2] = -p
            Q[row, row+1] = 4*p
            Q[row, row]   = (1+ht)-6*p
            Q[row, row-1] = 4*p
            Q[row, row-2] = -p
        end

        # * 组装f向量
        fk = zeros(row_Q)
        for row = 1:row
            fk[row] = ht^2*force(row*nx, k*nt)
        end

        # * 组装右端uk
        # ** 首先组装Qt矩阵
        # ** 再生成uk
        Qt = sparse(zeros(row_Q, row_Q))
        for row = 3:row_Q
            Qt[row, row+2] = p
            Qt[row, row+1] = -4p
            Qt[row, row]   = (2+ht)+6*p
            Qt[row, row-1] = -4p
            Qt[row, row-2] = p
        end
        uk = Qt*u[nt, :]

        # * 组装ukk:=u^{k-1}_{j}
        ukk = u[nt-1, :]

        # * 求解生成u^{k+1}
        u[nt+1, :] = Q\(fk+uk+ukk)
    end
end