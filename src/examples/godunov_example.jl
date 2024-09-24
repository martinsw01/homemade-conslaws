include("../godunov.jl")
include("../viz.jl")

function main()
    x_L, x_R = -1, 1
    f(U) = U^2/2
    df(U) = U
    ω = 0
    N = 100
    dx = (x_R - x_L) / (N + 1)
    x_mid = x_L + 0:dx:x_R
    # U0 = (x_mid .> 0) * 2 .- 1
    U0 = sin.(4π*x_mid)
    # BC(t, U) = [1;; 0]
    # BC(t, U) = [-1;; 1]
    BC(t, U) = [U[end];; U[1]]
    # BC(t, U) = [U[end];; U[1]]
    dt = 0.1
    T = 1
    # U, t = Godunov.godunovmethod(f, df, ω, U0, BC, dx, dt, 0.1)

    # @time 1+1

    U, t = Godunov.godunovmethod(f, df, ω, U0, BC, dx, dt, T)
    
    Viz.animate_solution((U,),
                        "Approximate solution",
                        x_mid, t)
end

main()