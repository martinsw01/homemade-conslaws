export NumericalFlux, stencil_size

"""
    NumericalFlux()
"""
abstract type NumericalFlux end

"""
    stencil_size(::NumericalFlux)

Gives a tuple of the number of left and right cells needed to compute the numerical flux.
"""
function stencil_size end

include("central_upwind.jl")
include("godunov.jl")
include("laxfriedrichs.jl")
include("rusanov.jl")