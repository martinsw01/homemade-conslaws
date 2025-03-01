export Grid, cells, cell_centers, for_each_cell, for_each_interior_cell, for_each_boundary_cell, reduce_cells, create_buffer

"""
    Grid

A type representing a grid. Contains the cell averages and determines how to loop over the grid.
"""
abstract type Grid end


"""
    cells(::Grid)

Returns the cell averages of the grid.
"""
function cells end

"""
    cell_centers(::Grid)

Returns the cell centers of the grid.
"""
function cell_centers end

"""
    for_each_cell(f, ::Grid)

Calls `f` on the cells of the grid and each cell index.
"""
function for_each_cell end

"""
    for_each_interior_cell(f, ::Grid; p=1, q=p)

Calls `f` on the cells and the index of every interior cell its `p` and `q` left and right neighbours, specified by the
stencil size of the numerical flux.
"""
function for_each_interior_cell end

"""
    for_each_boundary_cell(f, ::Grid, cells)
    for_each_boundary_cell(f, ::Grid, left_reconstruction, right_reconstruction)

Calls `f` on the the stencils centered at the boundary cells and their indices.
"""
function for_each_boundary_cell end

"""
    reduce_cells(f, ::Grid, initial_value)

Reduces the cells of the grid using the function `f` and the initial value `initial_value`. Used to compute the maximum
time step given by the CFL condition.
"""
function reduce_cells end



include("uniform1d.jl")
include("uniform1dwalls.jl")