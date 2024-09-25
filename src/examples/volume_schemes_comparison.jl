include("../finite_volumes.jl")
include("../viz.jl")

using LaTeXStrings

function main()
    x_L, x_R = -1, 1
    f(U) = U^2/2
    df(U) = U
    
    N = 50
    dx = (x_R - x_L) / (N + 1)
    dt = dx
    x_mid = x_L + 0:dx:x_R
    # U0 = x_mid .< 0
    # BC(t, U) = [1;; 0]
    U0 = sin.(4π * x_mid)
    BC(t, U) = [U[end]; U[1]]
    T = 1

    U_LxF, t_LxF = FiniteVolumes.lax_friedrichs_scheme(f, df, U0, BC, dx, dt, T)
    U_God, t_God = FiniteVolumes.godunov_scheme(f, df, 0, U0, BC, dx, dt, T)
    U_Rus, t_Rus = FiniteVolumes.rusanov_scheme(f, df, U0, BC, dx, dt, T)

    @assert length(t_LxF) == length(t_God) == length(t_Rus) "$(length(t_LxF)) ≠ $(length(t_God)) ≠ $(length(t_Rus))"

    Viz.animate_solution((U_LxF, U_God, U_Rus),#, x_mid' .< 0.5t),
                         [L"U_{LxF}" L"U_{God}" L"U_{Rus}"],# L"U_{exact}"],
                         x_mid, t_LxF, 3)

end

main()
