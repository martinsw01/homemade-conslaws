export LaxFriedrichsFlux

"""
    LaxFriedrichsFlux()

Lax-Friedrichs numerical flux. It is a simple numerical flux that uses the average of the fluxes
across the interfaces and adjusts by the maximum allowed wave speed ``\\frac{\\Dx}{\\Dt}``.

```math
    F^\\text{LxF}(U_L, U_R) = \\frac{1}{2} \\Big(f(U_L) + f(U_R)\\Big) - \\frac{\\Dx}{2\\Dt} (U_R - U_L)
```
"""
struct LaxFriedrichsFlux <: NumericalFlux end

(::LaxFriedrichsFlux)(eq::Equation, (left,), (right,), dx, dt) = 0.5*(eq(left) + eq(right)) - 0.5 * dx/dt * (right - left)

stencil_size(::LaxFriedrichsFlux) = (1, 1)