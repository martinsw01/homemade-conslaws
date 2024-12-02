export UniformGrid1DWalls


struct UniformGrid1DWalls{Float <: AbstractFloat, Int <: Integer} <: Grid1D{WallsBC{Float}}
    dx::Float
    bc::WallsBC{Float}
    cells::Vector{SVector{2, Float}}
    domain::Tuple{Float, Float}

    left_of_walls_indices::Vector{Int}
    right_of_walls_indices::Vector{Int}
    cells_not_containing_walls_indices::Vector{Int}
    interior_cells_indices::Vector{Int}

    function UniformGrid1DWalls(N, bc::BoundaryCondition, u0::Function, domain)
        x_L, x_R = domain
        dx = (x_R - x_L) / (N+1)

        cells = [SVector{2, typeof(dx)}(u0(x)...) for x in _cell_centers(x_L, x_R, dx)]

        left_indices, right_indices = _compute_wall_neighbour_indices(bc, N, dx, x_L, x_R)
        cells_not_containing_walls_indices = _compute_inner_cells_indices(left_indices, right_indices, N)
        inner_cells_indices = _compute_inner_cell_indices(cells_not_containing_walls_indices, left_indices, right_indices)

        return new{typeof(dx), typeof(N)}(
            dx, bc, cells, domain,left_indices, right_indices, cells_not_containing_walls_indices, inner_cells_indices
        )
    end
end


function _compute_inner_cell_indices(cells_not_containing_walls_indices, left_of_walls_indices, right_of_walls_indices)
    [i for i in cells_not_containing_walls_indices if i ∉ left_of_walls_indices && i ∉ right_of_walls_indices]
end


function _compute_wall_neighbour_indices(bc::WallsBC, N::Int, dx, x_L, x_R) where Int
    left_of_walls_indices = Vector{Int}(undef, length(bc.walls)+1)
    right_of_walls_indices = Vector{Int}(undef, length(bc.walls)+1)

    right_of_walls_indices[1] = 1
    left_of_walls_indices[end] = N + 1

    for (i, (left, right)) in enumerate(bc.walls)
        left_of_walls_indices[i] = floor(Int, (left - x_L) / dx)
        right_of_walls_indices[i+1] = ceil(Int, (right - x_L) / dx) + 1
    end
    return left_of_walls_indices, right_of_walls_indices
end


@views function _compute_number_of_indices_covered_by_walls(left_of_walls_indices, right_of_walls_indices)
    reduce(zip(left_of_walls_indices[1:end-1], right_of_walls_indices[2:end]); init=zero(eltype(left_of_walls_indices))) do acc, (left, right)
        acc + right - left - 1
    end
end


function _compute_inner_cells_indices(left_of_walls_indices, right_of_walls_indices, N)
    number_of_indices = _compute_number_of_indices_covered_by_walls(left_of_walls_indices, right_of_walls_indices)
    indices = Vector{typeof(N)}(undef, N + 1 - number_of_indices)

    i=1

    for (end_wall_left, start_wall_right) in zip(right_of_walls_indices, left_of_walls_indices)
        number_of_indices_between_walls = start_wall_right-end_wall_left+1
        indices[i:i+number_of_indices_between_walls-1] .= end_wall_left:start_wall_right
        i += number_of_indices_between_walls
    end

    return indices
end


function _cell_centers(x_L, x_R, dx)
    x_L + 0.5dx:dx:x_R
end


cell_centers(grid::UniformGrid1DWalls) = _cell_centers(grid.domain..., grid.dx)

reduce_cells(f, grid::UniformGrid1DWalls, initial_value) = _reduce_cells(f, grid, grid.cells_not_containing_walls_indices, initial_value)

for_each_cell(f, grid::UniformGrid1DWalls) = _for_each_cell(f, grid, grid.cells_not_containing_walls_indices)

for_each_interior_cell(f, grid::UniformGrid1DWalls; p=1, q=p) = _for_each_interior_cell(f, grid, grid.interior_cells_indices, p, q)

function for_each_boundary_cell(f, grid::UniformGrid1DWalls, left, right)
    for i in grid.left_of_walls_indices
        f((left[i-1], left[i], create_wall_ghost_cell(right[i])),
          (right[i-1], right[i], create_wall_ghost_cell(left[i])),
          i)
    end
    for i in grid.right_of_walls_indices
        f((create_wall_ghost_cell(right[i]), left[i], left[i+1]),
          (create_wall_ghost_cell(left[i]), right[i], right[i+1]),
          i)
    end
end

function for_each_boundary_cell(f, grid::UniformGrid1DWalls, cells)
    for i in grid.left_of_walls_indices
        f((cells[i-1], cells[i], create_wall_ghost_cell(cells[i])),
          i)
    end
    for i in grid.right_of_walls_indices
        f((create_wall_ghost_cell(cells[i]), cells[i], cells[i+1]),
          i)
    end
end

cells(grid::UniformGrid1DWalls) = grid.cells

get_dx(grid::UniformGrid1DWalls, cell_idx) = grid.dx