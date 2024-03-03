# 用差分法求解方程
include("force.jl")
using Plots, LinearAlgebra, SparseArrays
using BenchmarkTools
GC.enable(false)
# ! 还可以继续优化
# TODO 将前面的循环改成广播
# TODO 稀疏数组A √
# TODO 避免直接对大数组切片，用@view

# 网格划分
function Mainf()
    ll::Float64 = 2.93 # 燃料棒长度 
    tt::Float64 = 10 # 终止时间
    nx::Int = 31; 
    nt::Int = 1000001 
    h_x::Float64 = ll/(nx-1);
    h_t::Float64 = tt/(nt-1)
    x = (0:h_x:ll)';t = (0:h_t:tt)' # x与t方向的网格划分
    u::Matrix{Float64} = zeros(nx,nt) # 网格点形成的函数值数值解矩阵
    v::Matrix{Float64} = zeros(nx,nt) # 网格点形成的一阶导数值数值解矩阵
    u_exact::Matrix{Float64} = zeros(nx,nt) # 网格点形成的精确解矩阵
    F::Matrix{Float64} = zeros(nx,nt) # 右端项矩阵

    # 形成精确解矩阵
    d_exact = (x,t) -> t^2*sin(pi*x/ll)

    # TODO 可以把这个改为广播的形式
    @inbounds for i = 1:nx
        for j = 1:nt
            u_exact[i,j] = d_exact((i-1)*h_x,(j-1)*h_t)
        end
    end

    # 形成右端项矩阵
    @inbounds for i = 1:nx
        for j = 1:nt
            F[i,j] = Force((i-1)*h_x,(j-1)*h_t,ll)
        end
    end

    # 形成差分矩阵
    A = sparse(
        diagm(-6*ones(nx))+
        diagm(-2 => -1*ones(nx-2)) + diagm(2 => -1*ones(nx-2))+
        diagm(-1 => 4*ones(nx-1)) + diagm(1 => 4*ones(nx-1))
        )
    B = A[3:nx-2,:]

    # 根据差分法求解方程
    @inbounds for j = 2:nt
        v[3:nx-2,j] = (1-h_t) * v[3:nx-2,j-1] + h_t/(h_x^4) * B * u[:,j-1] + h_t * F[3:nx-2,j]
        u[3:nx-2,j] = h_t * v[3:nx-2,j-1] + u[3:nx-2,j-1]
        u[2,j] = (2/5) * (2 * u[3,j] - 0.5 * u[4,j])
        u[nx-1,j] = (2/5) * (2 * u[nx-2,j] - 0.5 * u[nx-3,j])
    end

    # * @view版本
    # ! bug 还需修改
    # @inbounds for j = 2:nt
    #     fill!(@view(v[3:nx-2,j]), (1-h_t) * v[3:nx-2,j-1] + h_t/(h_x^4) * B * u[:,j-1] + h_t * F[3:nx-2,j])
    #     fill!(@view(u[3:nx-2,j]), h_t * v[3:nx-2,j-1] + u[3:nx-2,j-1])
    #     fill!(@view(u[2,j]), (2/5) * (2 * u[3,j] - 0.5 * u[4,j]))
    #     fill!(@view(u[nx-1,j]), (2/5) * (2 * u[nx-2,j] - 0.5 * u[nx-3,j]))
    # end

    return 0
end
# error_matrix=double(u-u_exact)
# error_rate=double(abs(error_matrix ./ u_exact))
# surface(u)
@time Mainf()
GC.gc()