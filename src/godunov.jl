module Godunov

export godunovmethod

function godunovflux(f, U_l, U_r, ω)
    return max(f(max(U_l, ω)), f(min(U_r, ω)))
end

function cfl(df, U, dx)
    return dx/maximum(abs.(df.(U)))
end

"""
    godunovmethod(f, df, ω, U_0, BC, dx, dt, T)

Solves the scalar conservation law

``∂_t U + ∂_x f(U) = 0``

for a given flux function ``f`` with a unique minimizer ``ω`` using the Godunov method.

# Arguments
- `f`: Flux function with unique minimizer.
- `df`: Derivative of the flux.
- `ω`: The unique minimizer of the flux function.
- `U_0`: Initial condition array.
- `BC`: Boundary condition function. Takes time and the current state as arguments.
- `dx`: Spatial step size.
- `dt`: Maximal step size.
- `T`: Total simulation time.

# Returns
- Array representing the solution at each time step.
"""
function godunovmethod(f, df, ω, U_0, BC, dx, dt, T)
    U = copy(U_0)
    t = [0.]
    U_n = zeros(length(U[:, 1])+2)
    while t[end] < T
        U_n[2:end-1] = (U[:,end])
        U_n[[1 end]] = BC(t[end], U[:,end])
        dt_next = min(dt, cfl(df, U_n, dx))
        
        F = godunovflux.(f, U_n[1:end-1], U_n[2:end], ω)

        U_next = U_n[2:end-1] - dt_next * (F[2:end] - F[1:end-1]) / dx

        U = hcat(U, U_next)

        append!(t, t[end] + dt_next)
    end
    return U', t
end

end # module