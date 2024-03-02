using DifferentialEquations
include("./parameters.jl")
include("force.jl")
include("bound.jl")
include("equations.jl")

function Mainf()
    u0 = BoundaryCondition(x_discrete)
    prob = ODEProblem(FinDiffEqua, u0, timespan, p)
end

Mainf()