include("../upwind.jl")
include("../viz.jl")

function main()
    x_L, x_R = 0, 1
    T = 1
    N = 100
    dx = (x_R - x_L)/N
    dt = 0.9*dx
    x = x_L+dx:dx:x_R
    t = dt:dt:T
    M = length(t)
    
    u0(x) = sin(2*pi*x)

    A, f = upwind.upwind_scheme(N, M, dx, dt, u0)

    u = A\f


    Viz.animate_solution(reshape(u, N, M),
                         (x, t) -> u0(x - t),
                         x, t)
    # Viz.plot_matrix(A, N, M)

end

main()