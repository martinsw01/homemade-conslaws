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

# Example
```julia
# Set all cells to zero
for_each_cell(grid) do cells, idx
    cells[idx] = zero(eltype(cells))
end
```
"""
function for_each_cell end

"""
    for_each_interior_cell(f, ::Grid; p=1, q=p)

Calls `f` on the cells and the index of every interior cell its `p` and `q` left and right neighbours, specified by the
stencil size of the numerical flux.

# Example

Compute the forward difference reconstruction of the interior cells:

```julia
left_reconstruction, right_reconstruction = ... # Preallocated buffers
for_each_interior_cell(grid; p=1, q=1) do cells, (left_idx, center_idx, right_idx)
    left, center, right = cells[left_idx], cells[center_idx], cells[right_idx]
    left_reconstruction[idx] = center - 0.5 * (right - center)
    right_reconstruction[idx] = center + 0.5 * (right - center)
end
````
"""
function for_each_interior_cell end

"""
    for_each_boundary_cell(f, ::Grid, cells)
    for_each_boundary_cell(f, ::Grid, left_reconstruction, right_reconstruction)

Calls `f` on the the stencils centered at the boundary cells and their indices.

# Example

Compute the forward difference reconstruction of the boundary cells 
```julia
left_reconstruction, right_reconstruction = ... # Preallocated buffers
for_each_boundary_cell(grid, cells(grid)) do (left, center, right), idx
    left_reconstruction[idx] = center - 0.5 * (right - center)
    right_reconstruction[idx] = center + 0.5 * (right - center)
end
```

Compute the temporal derivative of the cell averages at the boundary cells.
```julia
out = ... # Preallocated buffer
for_each_boundary_cell(grid, left_reconstruction, right_reconstruction) do (lleft, lcenter, lright), (rleft, rcenter, rright), idx
    left_flux = F(rcenter, lright, dx, dt)
    right_flux = F(rleft, lcenter, dx, dt)
    out[idx] = (left_flux - right_flux) / dx
end
```
"""
function for_each_boundary_cell end

"""
    reduce_cells(f, ::Grid, initial_value)

Reduces the cells of the grid using the function `f` and the initial value `initial_value`. Used to compute the maximum
time step given by the CFL condition.

# Example

```julia
# Compute the sum of all cells
reduce_cells(grid, 0.) do acc, cells, idx
    acc + cells[idx]
end
```
"""
function reduce_cells end


"""
    create_buffer(::Grid)

Preallocates a buffer of the same shape as the cells of the grid.
"""
function create_buffer end



include("uniform1d.jl")
include("uniform1dwalls.jl")