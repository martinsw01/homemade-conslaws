import Base: convert, promote_rule
import Base: +, *, -, /, ^
import Base: ≈, isless

export DualNumber

struct DualNumber{d, T <: Number} <: Number
    real::T
    dual::SVector{d, T}

    function DualNumber(real::T, dual::AbstractVector{S}) where {T, S}
        d = length(dual)
        ST = promote_type(T, S)
        new{d, ST}(real, SVector{d, ST}(dual))
    end
    function DualNumber(real::T, dual::SVector{d, T}) where {d, T}
        new{d, T}(real, dual)
    end
end

function Base.Real(x::DualNumber)
    x.real
end

function promote_rule(::Type{DualNumber{d, S}}, ::Type{DualNumber{d, T}}) where {d, S, T}
    DualNumber{d, promote_type(S, T)}
end

function promote_rule(::Type{DualNumber{d, S}}, ::Type{T}) where {d, S, T<:Real}
    DualNumber{d, promote_type(T, S)}
end

function convert(::Type{DualNumber{d, S}}, x::DualNumber{d, T}) where {d, S, T}
    DualNumber(convert(S, x.real), SVector{d, S}(x.dual))
end

function convert(::Type{DualNumber{d, S}}, x::T) where {d, S, T<:Number}
    DualNumber(convert(S, x), zeros(SVector{d, S}))
end

function ≈(x::DualNumber{d, S}, y::DualNumber{d, T}; kwargs...) where {d, S, T}
    isapprox(x.real, y.real; kwargs...) && isapprox(x.dual, y.dual; kwargs...)
end
function ≈(x::DualNumber, y::T; kwargs...) where T <: Real
    isapprox(x.real, y; kwargs...) && all(isapprox.(x.dual, zero(eltype(x.dual)); kwargs...))
end
function ≈(x::T, y::DualNumber; kwargs...) where T <: Real
    y ≈ x
end 

function +(a::DualNumber{d, S}, b::DualNumber{d, T}) where {d, S, T}
    DualNumber(a.real + b.real, a.dual + b.dual)
end
function +(a::DualNumber{d, S}, b::T) where {d, S, T<:Real}
    DualNumber(a.real + b, a.dual)
end
function +(a::T, b::DualNumber{d, S}) where {d, S, T<:Real}
    b + a
end

function -(a::DualNumber{d, S}) where {d, S}
    DualNumber(-a.real, -a.dual)
end
function -(a::DualNumber{d, S}, b::DualNumber{d, T}) where {d, S, T}
    a + (-b)
end
function -(a::DualNumber, b::T) where T<:Real
    a + (-b)
end
function -(a::T, b::DualNumber) where T<:Real
    b + (-a)
end

function *(a::DualNumber{d, S}, b::DualNumber{d, T}) where {d, S, T}
    DualNumber(a.real * b.real, a.real .* b.dual + a.dual .* b.real)
end
function *(a::DualNumber, b::T) where T<:Real
    DualNumber(a.real * b, a.dual * b)
end
function *(a::T, b::DualNumber) where T<:Real
    b * a
end

function /(a::DualNumber{d, S}, b::DualNumber{d, T}) where {d, S, T}
    DualNumber(a.real / b.real, (a.dual * b.real - a.real * b.dual) / b.real^2)
end
function /(a::DualNumber, b::T) where T<:Real
    DualNumber(a.real / b, a.dual / b)
end
function /(a::T, b::DualNumber) where T<:Real
    DualNumber(a / b.real, -a * b.dual / b.real^2)
end

function ^(a::DualNumber{d, S}, b::DualNumber{d, T}) where {d, S, T}
    exp(b * log(a))
end
Base.:^(a::DualNumber, b::Real) = DualNumber(a.real^b, b * a.real^(b-1) * a.dual)
Base.:^(a::DualNumber, b::Integer) = DualNumber(a.real^b, b * a.real^(b-1) * a.dual)
Base.:^(x::Real, y::DualNumber) = x^y.real * DualNumber(1, log(x) * y.dual)


for (f, df) in [(:exp, :exp),
    (:log, :inv),
    (:sin, :cos),
    (:cos, x -> -sin(x)),
    (:sqrt, x -> 0.5 / sqrt(x)),
    (:cbrt, x -> 1 / (3 * cbrt(x)^2)),
    (:tan, x -> 1 + tan(x)^2),
    (:asin, x -> 1 / sqrt(1 - x^2)),
    (:acos, x -> -1 / sqrt(1 - x^2)),
    (:atan, x -> 1 / (1 + x^2)),
    (:sinh, :cosh),
    (:cosh, :sinh),
    (:tanh, x -> 1 - tanh(x)^2),
    (:asinh, x -> 1 / sqrt(x^2 + 1)),
    (:acosh, x -> 1 / sqrt(x^2 - 1)),
    (:atanh, x -> 1 / (1 - x^2)),
    (:inv, x -> -1 / x^2),
    (:abs, :sign)]
    @eval Base.$f(x::DualNumber) = DualNumber($f(x.real), $df(x.real) * x.dual)
end


function Base.isless(x::DualNumber{d, S}, y::DualNumber{d, T}) where {d, S, T}
    isless(x.real, y.real)
end
function Base.isless(x::DualNumber, y::Real)
    isless(x.real, y)
end
function Base.isless(x::Real, y::DualNumber)
    isless(x, y.real)
end