export BurgersEQ

struct BurgersEQ <: Equation end

(eq::BurgersEQ)(u) = 0.5 * u .^ 2

compute_max_abs_speed(::BurgersEQ, U) = abs(U)
