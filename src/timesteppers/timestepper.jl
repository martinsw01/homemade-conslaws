export TimeStepper, integrate!


"""
    TimeStepper

Specifies how to advance the solution forward in time.
"""
abstract type TimeStepper end


"""
    integrate!(::Grid, ::TimeStepper, system, compute_max_dt)

Advance the solution forward in time, mutating the grid in place. Returns the time step size used.
"""

include("forwardeuler.jl")
include("rk2.jl")