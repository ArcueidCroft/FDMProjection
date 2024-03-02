using DifferentialEquations
using LinearAlgebra
using SparseArrays
include("parameters.jl")
include("force.jl")

# * function of fdm-equations
# * the original equations:
# * ∂u/∂t = v
# * ∂v/∂t = f - v - (∂4)u/∂(t4), f = f(x, t)
# * du => ∂u/∂t or ∂v/∂t
# * u  => solutions
# * p  => parameters in fdm-equations
# * t  => continuous variable
function FinDiffEqua(du, u, p, t)
    # N = p # 离散节点数

    @inbounds for i = 1:N
        x = x_discrete[i]
        # * u[:, 1] = u
        # * u[:, 2] = v
        # * du[:, 1] = ∂u/∂t = v
        # * du[:, 2] = ∂v/∂t
        du[i, 1] = u[i, 2]
        du[i, 2] = Force(x, t)- u[i, 2] -hx^(-4)*(
            u[i+2, 1]-4*u[i+1, 1]+6*[i, 1]-4*u[i-1, 1]+u[i-2, 1]
            )
    end
end