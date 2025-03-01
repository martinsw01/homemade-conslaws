export GodunovFlux

"""
    GodunovFlux{Float <: AbstractFloat}(ω::Float)

An exact Riemann solver assuming the flux function has a unique minimum at ``ω``.

```math
    F^\\text{God}(U_L, U_R) = \\max\\Big\\{f\\big(\\max\\{U_L, \\omega\\}\\big), f\\big(\\min\\{U_L, \\omega\\}\\big)\\Big\\}
```
"""
struct GodunovFlux{Float <: AbstractFloat} <: NumericalFlux
    ω::Float
end

function (F::GodunovFlux)(eq::Equation, (left,), (right,), dx, dt)
    max(eq(max(left, F.ω)),
        eq(min(right, F.ω)))
end

stencil_size(::GodunovFlux) = 1, 1