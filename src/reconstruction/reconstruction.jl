export Reconstruction, reconstruct

"""
    Reconstruction

Determines how to reconstruct the solution at the cell interfaces from the cell averages.
"""
abstract type Reconstruction end

"""
    reconstruct(::Reconstruction, ::Grid, cell_averages)

Returns the left and right reconstructed values at the cell interfaces, of the same shape as `cell_averages`.
"""
function reconstruct end

include("noreconstruction.jl")
include("linear_reconstruction.jl")