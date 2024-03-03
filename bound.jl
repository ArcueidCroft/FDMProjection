include("exactsolution.jl")

function BoundaryCondition(discrete_vector)
    n = length(discrete_vector)
    u = zeros(n, 2)
    for i = 1:n
        x = discrete_vector[i]
        u[i, 1] = 0
        u[i, 2] = 0#TimeDerivative(x[i], 0)
    end
    return u
end

function PeriodicCondition(location, N)
    true_location = location
    if location == N+1
        true_location = 1
    elseif location == N+2
        true_location = 2
    elseif location == 0
        true_location = N
    elseif location == -1
        true_location = N-1
    end

    return true_location
end