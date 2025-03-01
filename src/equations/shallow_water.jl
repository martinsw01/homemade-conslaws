using StaticArrays

export ShallowWater1D

"""
    ShallowWater1D(g)

Shallow water equations in 1D with gravity `g`. It is given by the flux function

```math
f(h, hu) = \\begin{pmatrix} hu \\\\ \\frac{hu^2}{h} + \\frac{1}{2} g h^2 \\end{pmatrix}
```

where ``h`` is the water height and ``hu`` is the momentum.
"""
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