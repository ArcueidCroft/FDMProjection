function BoundaryCondition(discrete_vector)
    n = length(discrete_vector)
    u = zeros(n, 2)
    for i = 1:n
        x = discrete_vector[i]
        u[i, 1] = 0.02
        u[i, 2] = 0.01
    end
    return u
end

function PeriodicCondition(a, N)

end