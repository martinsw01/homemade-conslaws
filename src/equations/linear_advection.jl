export ConstantLinearAdvection

struct ConstantLinearAdvection{RealType <: Real} <: Equation
    a::RealType
end

(eq::ConstantLinearAdvection)(u) = eq.a * u

compute_max_abs_speed(::ConstantLinearAdvection, U) = abs(a)