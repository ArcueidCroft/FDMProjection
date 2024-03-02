function BoundaryCondition(discrete_vector)
    n = length(discrete_vector)
    u = zeros(n, 2)
    for i = 1:n
        x = discrete_vector[i]
        u[i, 1] = 0
        u[i, 2] = 0
    end
    return u
end