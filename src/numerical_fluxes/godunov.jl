export GodunovFlux

struct GodunovFlux{Float <: AbstractFloat} <: NumericalFlux
    ω::Float
end

function (F::GodunovFlux)(eq::Equation, (left,), (right,), dx, dt)
    max(eq(max(left, F.ω)),
        eq(min(right, F.ω)))
end

stencil_size(::GodunovFlux) = 1, 1