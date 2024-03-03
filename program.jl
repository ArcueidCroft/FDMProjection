using DifferentialEquations
using Plots
include("./parameters.jl")
include("force.jl")
include("bound.jl")
include("equations.jl")

function Mainf()
    N, x_discrete, l, hx, timespan = GetParameters()
    u0 = BoundaryCondition(x_discrete)
    prob = ODEProblem(FinDiffEqua, u0, timespan, p)
    t = solve(prob)
    solution = t[:, 2, :]

    # * Plot image with Plots.jl
    Plots.heatmap(solution)
end

Mainf()