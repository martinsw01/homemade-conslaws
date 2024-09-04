include("../upwind.jl")
include("../viz.jl")

function main()
    x_L, x_R = 0, 1
    T = 1
    N, M = 100, 110
    dx = (x_R - x_L)/N
    dt = T/M
    u0(x) = sin(2*pi*x)

    A, f = upwind.upwind_scheme(N, M, dx, dt, u0)

    u = A\f

    x = (1:N)/N
    t = (1:M)/M*T


    Viz.animate_solution(reshape(u, N, M), N, M)
    # Viz.plot_matrix(A, N, M)

end

main()