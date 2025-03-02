export Equation, compute_eigenvalues, compute_max_abs_eigenvalue


"""
    Equation

Describes the flux term ``f`` in a conservation law ``u_t + f(u)_x = 0``. Should be a callable object of the form `(::Equation)(U)`
computing the flux at `U`.

# Example

```julia-console
julia> eq = BurgersEQ();
julia> eq.(-2:2)
5-element Vector{Float64}:
 2.0
 0.5
 0.0
 0.5
 2.0
```
"""
abstract type Equation end

"""
    compute_eigenvalues(eq::Equation, U)

Compute the eigenvalues of the flux Jacobian ``f'(U)`` at `U` for the equation `eq`.

# Example

```julia-console
julia> eq = BurgersEQ();
julia> compute_eigenvalues(eq, -2)
-2
julia> compute_eigenvalues(eq, 0)
0
julia> compute_eigenvalues(eq, 2)
2
```
"""
function compute_eigenvalues end

"""
    compute_max_abs_eigenvalue(eq::Equation, U)

Compute the maximum absolute eigenvalue of the flux Jacobian ``f'(U)`` at `U` for the equation `eq`. Used to compute the time step
given the CFL condition.

# Example
```julia-console
julia> eq = BurgersEQ();
julia> compute_max_abs_eigenvalue(eq, -2)
-2
julia> compute_max_abs_eigenvalue(eq, 0)
0
julia> compute_max_abs_eigenvalue(eq, 2)
2
```
"""
function compute_max_abs_eigenvalue end


include("burgers.jl")
include("shallow_water.jl")