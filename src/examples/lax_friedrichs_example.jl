include("../lax_friedrichs.jl")
include("../viz.jl")

function main()
    x_L, x_R = -1, 1
    f(U) = U^2/2
    df(U) = U

    N = 50
    dx = (x_R - x_L) / (N + 1)
    x_mid = x_L + 0:dx:x_R
    # U0 = sin.(4π*x_mid)
    U0 = x_mid .< 0
    # U0 = (x_mid .> 0) * 2 .- 1

    # BC(t, U) = [U[end];; U[1]]
    BC(t, U) = [1;; 0]
    # BC(t, U) = [-1;; 1]
    dt = dx
    T = 1
    U, t = LaxFriedrichs.lax_friedrichs_method(f, df, U0, BC, dx, dt, T)
    
    Viz.animate_solution(U,
                        (x, t) -> x < 0.5t,
                        x_mid, t)
    
    # Viz.animate_solution((U,),
    #                     "Approximate solution",
    #                     x_mid, t)
end

main()