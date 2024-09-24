module Rusanov

export rusanov_method

function rusanov_flux(U_l, U_r, f, df)
    return 0.5 * (f(U_l) + f(U_r)) - 0.5 * max(abs(df(U_l)), abs(df(U_r))) * (U_r - U_l)
end

function cfl(df, U, dx)
    return dx / maximum(abs.(df.(U)))
end


function rusanov_method(f, df, U_0, BC, dx, dt, T)
    U = copy(U_0)
    t = [0.]
    U_n = zeros(length(U[:, 1])+2)
    while t[end] < T
        U_n[2:end-1] = (U[:,end])
        U_n[[1 end]] = BC(t[end], U[:,end])
        dt_next = min(dt, cfl(df, U_n, dx))
        
        F = rusanov_flux.(U_n[1:end-1], U_n[2:end], f, df)

        U_next = U_n[2:end-1] - dt_next * (F[2:end] - F[1:end-1]) / dx

        U = hcat(U, U_next)

        append!(t, t[end] + dt_next)
    end
    return U', t
end

end # module