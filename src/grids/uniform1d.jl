using QuadGK, StaticArrays

export Grid1D, UniformGrid1D, for_each_cell, for_each_interior_cell, for_each_boundary_cell, get_neighbour_indices, get_dx, cells, separate_variables,
    cell_centers


abstract type Grid1D{BC <: BoundaryCondition} <: Grid end

"""
    UniformGrid1D{BC <: BoundaryCondition, Float <: AbstractFloat, Cell} <: Grid1D{BC}
"""
struct UniformGrid1D{BC <: BoundaryCondition, Float <: AbstractFloat, Cell} <: Grid1D{BC}
    dx::Float
    bc::BC
    cells::Vector{Cell}
    domain::Tuple{Float, Float}

    function UniformGrid1D(N, bc, u0, domain)
        x_L, x_R = domain
        dx = (x_R - x_L) / N
        return new{typeof(bc), typeof(dx), eltype(u0)}(
            dx, bc, u0[:], domain
        )
    end
    function UniformGrid1D(N, bc, u0::AbstractVector{<:AbstractVector{T}}, domain) where T
        x_L, x_R = domain
        dx = (x_R - x_L) / N
        d = length(u0[1])
        return new{typeof(bc), typeof(dx), SVector{d, T}}(
            dx, bc, SVector{d, T}.(u0), domain
        )
    end
end

function cell_centers(grid::UniformGrid1D)
    x_L, x_R = grid.domain
    x_L + 0.5grid.dx:grid.dx:x_R - 0.5grid.dx
end

function cell_centers(N, x_l, x_r)
    dx = (x_r - x_l) / N
    x_l + 0.5dx:dx:x_r - 0.5dx
end

function _reduce_cells(f, grid::Grid1D, indices, initial_value)
    reduce(indices, init=initial_value) do acc, i
        f(acc, cells(grid), i)
    end
end

reduce_cells(f, grid::Grid1D, initial_value) = _reduce_cells(f, grid, eachindex(cells(grid)), initial_value)

function create_buffer(grid::Grid1D)
    similar(cells(grid))
end


function _for_each_cell(f, grid::Grid1D, indices)
    for i in indices
        f(cells(grid), i)
    end
end

for_each_cell(f, grid::Grid1D) = _for_each_cell(f, grid, eachindex(cells(grid)))

function _for_each_interior_cell(f, grid::Grid1D, indices, p, q)
    for cell_idx in indices
        ileft = cell_idx - p
        iright = cell_idx + q
        f(cells(grid), ileft, cell_idx, iright)
    end
end

for_each_interior_cell(f, grid::UniformGrid1D; p=1, q=p) = _for_each_interior_cell(f, grid, eachindex(cells(grid))[p+1:end-q], p, q)

for_each_boundary_cell(f, grid::Grid1D) = for_each_boundary_cell(f, grid, cells(grid))

function for_each_boundary_cell(f, ::Grid1D{NeumannBC}, cells)
    f((cells[1], cells[1], cells[2]),
      1)
    f((cells[end-1], cells[end], cells[end]),
      lastindex(cells))
end

function for_each_boundary_cell(f, ::Grid1D{NeumannBC}, left_reconstruction, right_reconstruction)
    f((left_reconstruction[1], right_reconstruction[1], right_reconstruction[2]),
      (left_reconstruction[1], left_reconstruction[1], left_reconstruction[2]),
      1)
    f((left_reconstruction[end-1], left_reconstruction[end], right_reconstruction[end]),
      (right_reconstruction[end-1], right_reconstruction[end], right_reconstruction[end]),
      lastindex(left_reconstruction))
end

function for_each_boundary_cell(f, ::Grid1D{PeriodicBC}, left, right)
    f((left[end], left[1], left[2]), (right[end], right[1], right[2]), 1)
    f((left[end-1], left[end], left[1]), (right[end-1], right[end], right[1]), lastindex(left))
end

function for_each_boundary_cell(f, ::Grid1D{PeriodicBC}, cells)
    f((cells[end], cells[1], cells[2]), 1)
    f((cells[end-1], cells[end], cells[1]), lastindex(cells))
end

function create_wall_ghost_cell(Q)
    @SVector [Q[1], -Q[2]]
end

function for_each_boundary_cell(f, ::Grid1D{WallBC}, left, right)
    f((create_wall_ghost_cell(right[1]), left[1], left[2]),
      (create_wall_ghost_cell(left[1]), right[1], right[2]),
      1)
    f((left[end-1], left[end], create_wall_ghost_cell(right[end])),
      (right[end-1], right[end], create_wall_ghost_cell(left[end])),
      lastindex(left))
end
function for_each_boundary_cell(f, ::Grid1D{WallBC}, cells)
    f((create_wall_ghost_cell(cells[1]), cells[1], cells[2]),
      1)
    f((cells[end-1], cells[end], create_wall_ghost_cell(cells[end])),
      lastindex(cells))
end


function get_dx(grid::UniformGrid1D, cell_idx)
    grid.dx
end

function separate_variables(U::AbstractMatrix{SVector{N, T}}) where {N, T}
    stacked = stack(U)
    (stacked[i,:,:]' for i in 1:N)
end

function separate_variables(U::AbstractVector{SVector{N, T}}) where {N, T}
    [U[i][1] for i in eachindex(U)], [U[i][2] for i in eachindex(U)]
end


cells(grid::UniformGrid1D) = grid.cells
