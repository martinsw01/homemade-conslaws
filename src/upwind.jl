module upwind

function upwind_scheme(x, t, a::Function, u0::Function)
    u = zeros(length(t) + 1, length(x))
    u[1,:] = u0.(x)

    dt = t[2] - t[1]
    dx = x[2] - x[1]
    nu = dt/dx

    for n in eachindex(t)
        a_eval = a.(x, t[n])
        u[n+1,:] = nu * a_eval .* circshift(u[n,:], 1) + (1 .- nu * a_eval) .* u[n,:]
    end

    return @view u[2:end, :]
end

end # module