module central_difference

using LinearAlgebra, BlockBandedMatrices, Plots

export central_difference_scheme


function load_vector(N, M, dx, dt, u0)
    u = zeros(N * M)
    nu = dt/(2*dx)

    x = (0:N+1) * dx
    u0eval = u0.(x)

    u[1:N] = u0eval[2:end-1] + nu * (u0eval[1:end-2] - u0eval[3:end])

    return u
end


function create_lower(N, nu)
    L = zeros(N, N)
    L += Tridiagonal(fill(-nu, N-1), fill(-1., N), fill(nu, N-1))
    L[1, N] = -nu
    L[N, 1] = +nu
    return L
end

function create_upper(N)
    return zeros(N, N)
end

function create_diag(N)
    D = zeros(N, N)
    return D + Diagonal(fill(1., N))
end


function central_difference_matrix(N, M, dx, dt)
    return BlockTridiagonal(fill(create_lower(N, dt/(2*dx)), M-1),
                            fill(create_diag(N), M),
                            fill(create_upper(N), M-1))
end


function central_difference_scheme(N::Integer, M::Integer, dx::Real, dt::Real, u0::Function)
    A = central_difference_matrix(N, M, dx, dt)
    f = load_vector(N, M, dx, dt, u0)

    return A, f
end

end # module