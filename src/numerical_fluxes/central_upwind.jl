export CentralUpwind

struct CentralUpwind <: NumericalFlux end

function (::CentralUpwind)(eq::Equation, (left,), (right,), dx, dt)
    eig_left = compute_eigenvalues(eq, left)
    eig_right = compute_eigenvalues(eq, right)
    _zero = zero(eltype(left))
    aplus = max(eig_left..., eig_right..., _zero)
    aminus = min(eig_left..., eig_right..., _zero)

    (aplus .* eq(left) .- aminus .* eq(right) .+ aplus .* aminus .* (right .- left)) ./ (aplus .- aminus)
end

number_of_cells(::CentralUpwind) = 1