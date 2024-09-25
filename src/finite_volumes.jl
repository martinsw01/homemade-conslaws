module FiniteVolumes

using ElasticArrays

export godunov_scheme, lax_friedrichs_scheme, rusanov_scheme

function cfl(df, U, dx)
    # println(maximum(abs.(df.(U))))
    return dx/maximum(abs.(df.(U)))
end


@views function step!(U_n, U, F, BC, numerical_flux, df, dx, dt, t)
    U_n[2:end-1] = (U[:,end])
    U_n[[1 end]] = BC(t[end], U[:,end])
    dt_next = min(dt, cfl(df, U_n, dx))

    println(dt_next)
    
    F[:] = numerical_flux.(U_n[1:end-1], U_n[2:end], dt_next)

    U_next = U_n[2:end-1] - dt_next * (F[2:end] - F[1:end-1]) / dx

    append!(U, U_next)
    append!(t, t[end] + dt_next)
end


"""
    finite_volumes(numerical_flux, df, U_0, BC, dx, dt, T)

Solves the scalar conservation law

``∂_t U + ∂_x f(U) = 0``

for a given flux function ``f`` using the finite volumes method.

# Arguments
- `numerical_flux`: Numerical flux function. Takes the cell averages left and right of the interface and the time step as arguments.
- `df`: Derivative of the flux.
- `U_0`: Initial condition array.
- `BC`: Boundary condition function. Takes time and the current state as arguments.
- `dx`: Spatial step size.
- `dt`: Maximal temporal step size.
- `T`: Total simulation time.

# Returns
- Matrix [time, space] representing the solution at each time step.
- Vector of time points.
"""
@views function finite_volumes(numerical_flux::Function, df, U_0, BC, dx, dt, T)
    U = ElasticMatrix(reshape(U_0 .+ 0., :, 1))
    t = ElasticVector([0.])
    U_n = zeros(length(U[:, 1])+2)
    F = zeros(length(U[:, 1])+1)
    while t[end] < T
        step!(U_n, U, F, BC, numerical_flux, df, dx, dt, t)
    end
    return U', t
end


"""
    godunov_scheme(f, df, ω, U_0, BC, dx, dt, T)

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
- Matrix [time, space] representing the solution at each time step.
- Vector of time points.
"""
function godunov_scheme(f, df, ω, U_0, BC, dx, dt, T)
    @assert df(ω) ≈ 0 "ω=$ω is not a minimizer of f"
    finite_volumes(df, U_0, BC, dx, dt, T) do U_l, U_r, dt
        return max(f(max(U_l, ω)), f(min(U_r, ω)))
    end
end


"""
    lax_friedrichs_scheme(f, df, U_0, BC, dx, dt, T)

Solves the scalar conservation law

``∂_t U + ∂_x f(U) = 0``

for a given flux function ``f`` using Lax-Friedrichs.

# Arguments
- `numerical_flux`: Numerical flux function. Takes the cell averages left and right of the interface and the time step as arguments.
- `df`: Derivative of the flux.
- `U_0`: Initial condition array.
- `BC`: Boundary condition function. Takes time and the current state as arguments.
- `dx`: Spatial step size.
- `dt`: Maximal temporal step size.
- `T`: Total simulation time.

# Returns
- Matrix [time, space] representing the solution at each time step.
- Vector of time points.
"""
function lax_friedrichs_scheme(f, df, U_0, BC, dx, dt, T)
    finite_volumes(df, U_0, BC, dx, dt, T) do U_l, U_r, dt
        return 0.5 * (f(U_l) + f(U_r)) - 0.5 * dx / dt * (U_r - U_l)
    end
end


"""
    rusanov_scheme(f, df, U_0, BC, dx, dt, T)

Solves the scalar conservation law

``∂_t U + ∂_x f(U) = 0``

for a given flux function ``f`` using Rusanov.

# Arguments
- `numerical_flux`: Numerical flux function. Takes the cell averages left and right of the interface and the time step as arguments.
- `df`: Derivative of the flux.
- `U_0`: Initial condition array.
- `BC`: Boundary condition function. Takes time and the current state as arguments.
- `dx`: Spatial step size.
- `dt`: Maximal temporal step size.
- `T`: Total simulation time.

# Returns
- Matrix [time, space] representing the solution at each time step.
- Vector of time points.
"""
function rusanov_scheme(f, df, U_0, BC, dx, dt, T)
    finite_volumes(df, U_0, BC, dx, dt, T) do U_l, U_r, dt
        return 0.5 * (f(U_l) + f(U_r)) - 0.5 * max(abs(df(U_l)), abs(df(U_r))) * (U_r - U_l)
    end
end

end # module

