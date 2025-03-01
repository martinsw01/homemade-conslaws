export BurgersEQ

"""
    BurgersEQ

The standard Burgers equation ``u_t + \\qty(\\frac12 u^2)_x = 0`` used as a test case for hyperbolic PDE solvers.
"""
struct BurgersEQ <: Equation end

(eq::BurgersEQ)(u) = 0.5 * u .^ 2


compute_eigenvalues(::BurgersEQ, U) = U

compute_max_abs_eigenvalue(::BurgersEQ, U) = abs(U)