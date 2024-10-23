export LaxFriedrichsFlux

struct LaxFriedrichsFlux <: NumericalFlux end

(::LaxFriedrichsFlux)(eq::Equation, (left,), (right,), dx, dt) = 0.5*(eq(left) + eq(right)) - 0.5 * dx/dt * (right - left)

number_of_cells(::LaxFriedrichsFlux) = 1