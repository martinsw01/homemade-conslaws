module central_difference

export central_difference_scheme

function central_difference_scheme(x, t, u0::Function)
    u = zeros(length(t) + 1, length(x))
    u[1,:] = u0.(x)

    dt = t[2] - t[1]
    dx = x[2] - x[1]
    nu = dt/(2*dx)

    for n in 1:length(t)
        u[n+1,:] = u[n,:] - nu * (circshift(u[n,:], -1) - circshift(u[n,:], 1))
    end

    return @view u[2:end,:]
end



end # module