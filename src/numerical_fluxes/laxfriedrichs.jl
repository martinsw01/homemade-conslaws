export LaxFriedrichsFlux

struct LaxFriedrichsFlux <: NumericalFlux end

(::LaxFriedrichsFlux)(eq::Equation, (left,), (right,), dx, dt) = 0.5*(eq(left) + eq(right)) - 0.5 * dx/dt * (right - left)

stencil_size(::LaxFriedrichsFlux) = (1, 1)