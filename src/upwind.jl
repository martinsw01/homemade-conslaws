module upwind

using LinearAlgebra, BlockBandedMatrices

export upwind_scheme


"""
The N first entries determines the value of u after the first time step.
Zero otherwise.
"""
function load_vector(N, M, dx, dt, u0)
    u = zeros(N * M)
    nu = dt/dx

    x = (1:N) * dx

    u[1:N] = nu * u0.(x .- dx) + (1 - nu) * u0.(x)

    return u
end


"""
Creates the lower block of the upwind matrix. Determines the coefficients of the terms in the previous
time step ``u^{n-1}`` described in the equation (2.14) in numcl_notes. We have a periodic boundary condition, which
explains the entry in the top right corner.
"""
function create_lower(N, nu)
    L = zeros(N, N)
    L += Tridiagonal(fill(-nu, N-1),
                     fill(nu-1, N),
                     fill(0., N-1))
    L[1, N] = -nu
    return L
end

"""
Upper block of the upwind matrix. Determines the coefficients of the terms in the next time step
``u^{n+1}`` of equation (2.14).
This is an explicit scheme, so these are all zeros.
"""
function create_upper(N)
    return zeros(N, N)
end

"""
Diagonal block of the upwind matrix. Determines the coefficients of the terms in the current time step
``u^{n}`` of equation (2.14).
"""
function create_diag(N)
    D = zeros(N, N)
    return D + Diagonal(fill(1., N))
end

"""
Creates a block tridiagonal matrix of the above three blocks. Sets up the linear system described in
equation (2.14) for all time steps.
"""
function upwind_matrix(N, M, dx, dt)
    return BlockTridiagonal(fill(create_lower(N, dt/dx), M-1),
                            fill(create_diag(N), M),
                            fill(create_upper(N), M-1))
end

"""
Sets up the linear system described in equation (2.14) for all time steps. Assumes a periodic boundary
condition (see the equations below (2.8)). Here my solution differ slightly from the one in the notes.
As with the notes, I assign ``x_L = x_0 < x_1 < ... < x_N = x_R``. However, the notes says that
``u_0 = u_{N-1}`` and ``u_N = u_1``. This, I believe, is a typo or a mistake. I use instead
``u_{-1} = u_{N-1}`` and ``u_{0} = u_{N}``. 
As the solution is periodic, I solve for ``u_j`` for all ``j=1,...,N``, where ``x_0``and ``x_N`` are the
boundary points.
"""
function upwind_scheme(N::Integer, M::Integer, dx::Real, dt::Real, u0::Function)
    A = upwind_matrix(N, M, dx, dt)
    f = load_vector(N, M, dx, dt, u0)

    return A, f
end

end # module
