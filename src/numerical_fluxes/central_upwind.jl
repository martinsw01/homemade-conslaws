export CentralUpwind

"""
    CentralUpwind()

Used for the shallow water equations. Well-balanced in the absence of dry land. Given by

```math
    F(Q_L, Q_R) = \\frac{a^+ F(Q_L) - a^- F(Q_R)}{a^+ - a^-} + \\frac{a^+ a^- (Q_R - Q_L)}{a^+ - a^-},
```
where ``a^+ = \\max(\\lambda(Q_L), \\lambda(Q_R), 0)`` and ``a^- = \\min(\\lambda(Q_L), \\lambda(Q_R), 0)`` are
the maximum and minimum eigenvalues of the Jacobian matrix of the flux function.
"""
struct CentralUpwind <: NumericalFlux end

function (::CentralUpwind)(eq::Equation, (left,), (right,), dx, dt)
    eig_left = compute_eigenvalues(eq, left)
    eig_right = compute_eigenvalues(eq, right)
    _zero = zero(eltype(left))
    aplus = max(eig_left..., eig_right..., _zero)
    aminus = min(eig_left..., eig_right..., _zero)

    (aplus .* eq(left) .- aminus .* eq(right) .+ aplus .* aminus .* (right .- left)) ./ (aplus .- aminus)
end

stencil_size(::CentralUpwind) = 1, 1