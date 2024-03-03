"""
Order: N, x_discrete, l, hx, timespan
"""
function GetParameters()
    N = 32
    x_discrete = range(0, stop = 1, length = N)
    l = 2.93
    hx = step(x_discrete)
    timespan = (0, 10)
    return (N, x_discrete, l, hx, timespan)
end

