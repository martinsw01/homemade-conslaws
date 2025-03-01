export NoReconstruction

"""
    NoReconstruction

Zeroth order reconstruction, i.e. the cell averages are used as the left and right values at the cell interfaces.
"""
struct NoReconstruction <: Reconstruction end

function reconstruct(::NoReconstruction, grid::Grid, cell_averages=cells(grid))
    return cell_averages, cell_averages
end