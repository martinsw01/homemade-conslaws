export BurgersEQ

struct BurgersEQ <: Equation end

(eq::BurgersEQ)(u) = 0.5 * u .^ 2


compute_eigenvalues(::BurgersEQ, U) = U

compute_max_abs_eigenvalue(::BurgersEQ, U) = abs(U)

conserved_variables(::BurgersEQ) = (:U,)