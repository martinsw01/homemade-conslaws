export TimeStepper, integrate!


"""
    TimeStepper

Specifies how to advance the solution forward in time.
"""
abstract type TimeStepper end


"""
    integrate!(::Grid, ::TimeStepper, system, compute_max_dt)

Advance the solution forward in time, mutating the cell averages in place. Returns the time step size used.
"""
function integrate! end

include("forwardeuler.jl")
include("rk2.jl")