export Equation, compute_eigenvalues, compute_max_abs_eigenvalue


"""
    Equation

Describes the flux term ``f`` in a conservation law ``u_t + f(u)_x = 0``. Should be a callable object of the form `(::Equation)(U)`
computing the flux at `U`.
"""
abstract type Equation end

"""
    compute_eigenvalues(eq::Equation, U)

Compute the eigenvalues of the flux Jacobian ``f'(U)`` at `U` for the equation `eq`.
"""
function compute_eigenvalues end

"""
    compute_max_abs_eigenvalue(eq::Equation, U)

Compute the maximum absolute eigenvalue of the flux Jacobian ``f'(U)`` at `U` for the equation `eq`. Used to compute the time step
given the CFL condition.
"""
function compute_max_abs_eigenvalue end


include("burgers.jl")
include("shallow_water.jl")