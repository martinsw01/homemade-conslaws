export RusanovFlux

"""
    RusanovFlux()

An approximate Riemann solver that uses the average of the flux accross the interfaces and adjusting using the maximum eigenvalue of the flux function.

```math
    F^\\text{Rus}(U_L, U_R) = \\frac{1}{2} \\Big(F(U_L) + F(U_R)\\Big) - \\lambda (U_R - U_L)
```
"""
struct RusanovFlux <: NumericalFlux end

function (::RusanovFlux)(eq::Equation, (left,), (right,), dx, dt)
    λ = max(compute_max_abs_eigenvalue(eq, left), compute_max_abs_eigenvalue(eq, right))
    0.5(eq(left) + eq(right) - λ * (right - left))
end

stencil_size(::RusanovFlux) = 1, 1