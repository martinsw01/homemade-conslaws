using StaticArrays

export ShallowWater1D

struct ShallowWater1D <: Equation
    g::Float64
end


function (swe::ShallowWater1D)(Q)
    h, hu = Q
    @SVector [hu, hu^2 / h + 0.5 * swe.g * h^2]
end

function compute_eigenvalues(swe::ShallowWater1D, Q)
    h, hu = Q
    u = hu / h
    @SVector [u - sqrt(swe.g * h), u + sqrt(swe.g * h)]
end

function compute_max_abs_eigenvalue(swe::ShallowWater1D, Q)
    maximum(abs.(compute_eigenvalues(swe, Q)))
end

function conserved_variables(::ShallowWater1D)
    :h, :hu
end