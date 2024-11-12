export RusanovFlux

struct RusanovFlux <: NumericalFlux end

function (::RusanovFlux)(eq::Equation, (left,), (right,), dx, dt)
    λ = max(compute_max_abs_eigenvalue(eq, left), compute_max_abs_eigenvalue(eq, right))
    0.5(eq(left) + eq(right) - λ * (right - left))
end

stencil_size(::RusanovFlux) = 1, 1