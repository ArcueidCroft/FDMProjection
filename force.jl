# include("parameters.jl")
function Force(x, t)
    p = sin(pi*x/l)
    return (pi/l)^4*t^2*p+2*t*p+2*p
end

function Force(x, t, l)
    p = sin(pi*x/l)
    return (pi/l)^4*t^2*p+2*t*p+2*p
end