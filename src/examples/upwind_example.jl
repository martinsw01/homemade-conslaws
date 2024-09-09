include("../upwind.jl")
include("../viz.jl")

function main()
    x_L, x_R = 0, 1
    T = 1
    N = 50
    dx = (x_R - x_L)/N
    dt = 0.9*dx
    x = x_L+dx:dx:x_R
    t = dt:dt:T

    u0(x) = sin(2*pi*x)
    a(x, t) = 1

    U = upwind.upwind_scheme(x, t, a, u0)

    Viz.animate_solution(U,
                        (x, t) -> u0(x - t),
                        x, t)

end

main()