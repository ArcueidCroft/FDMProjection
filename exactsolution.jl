using Flux
function ExactSolution(x, t)
    t^2*sin(pi*x/l)
end

function TimeDerivative(x, t)
    return gradient(ExactSolution, x, t)[2]
end