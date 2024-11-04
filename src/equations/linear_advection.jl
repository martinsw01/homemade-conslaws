export ConstantLinearAdvection

struct ConstantLinearAdvection{RealType <: Real} <: Equation
    a::RealType
end

(eq::ConstantLinearAdvection)(u) = eq.a * u

compute_eigenvalues(eq::ConstantLinearAdvection, U) = eq.a

compute_max_abs_eigenvalue(eq::ConstantLinearAdvection, U) = abs(eq.a)

conserved_variables(::ConstantLinearAdvection) = (:U,)