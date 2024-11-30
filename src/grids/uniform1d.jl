using QuadGK, StaticArrays

export UniformGrid1D, for_each_cell, for_each_inner_cell, for_each_boundary_cell, get_neighbour_indices, get_dx, inner_cells, separate_variables,
    cell_centers


abstract type Grid1D{BC <: BoundaryCondition} <: Grid end

struct UniformGrid1D{BC <: BoundaryCondition, Float <: AbstractFloat, Cell} <: Grid1D{BC}
    dx::Float
    bc::BC
    cells::Vector{Cell}
    domain::Tuple{Float, Float}

    function UniformGrid1D(N, bc::BoundaryCondition, u0::Function, domain)
        x_L, x_R = domain
        dx = (x_R - x_L) / (N+1)
        x = x_L:dx:x_R

        averages = _calc_average.(u0, x[1:end-1], x[2:end])

        if eltype(averages) <: Number
            return new{typeof(bc), typeof(dx), eltype(averages)}(
                dx, bc, averages, domain
            )
        else
            conserved_variables::Int64 = length(averages[1])
            FloatType = eltype(averages[1])
            cells = [SVector{conserved_variables, FloatType}(average) for average in averages]
            return new{typeof(bc), typeof(dx), eltype(cells)}(
                dx, bc, cells, domain
            )
        end
    end
end

function cell_centers(grid::UniformGrid1D)
    x_L, x_R = grid.domain
    x_L + 0.5grid.dx:grid.dx:x_R - 0.5grid.dx
end

function _calc_average(f, a, b)
    quadgk(f, a, b)[1] / (b - a)
end

function reduce_cells(f, grid::Grid1D, initial_value)
    inner_cells_indices = eachindex(cells(grid))
    reduce(inner_cells_indices, init=initial_value) do acc, i
        f(acc, cells(grid), i)
    end
end


function create_buffer(grid::Grid1D)
    similar(cells(grid))
end


function for_each_cell(f, grid::Grid1D)
    for cell_idx in eachindex(cells(grid))
        f(cells(grid), cell_idx)
    end
end


function for_each_inner_cell(f, grid::Grid1D; p=1, q=p)
    for cell_idx in eachindex(cells(grid))[p+1:end-q]
        ileft = cell_idx - p
        iright = cell_idx + q
        f(cells(grid), ileft, cell_idx, iright)
    end
end


for_each_boundary_cell(f, grid::Grid1D) = for_each_boundary_cell(f, grid, (cells(grid),))

function for_each_boundary_cell(f, ::Grid1D{NeumannBC}, inputs)
    f(((cells[1], cells[1], cells[2]) for cells in inputs)..., 1)
    f(((cells[end-1], cells[end], cells[end]) for cells in inputs)..., lastindex(inputs[1]))
end

function for_each_boundary_cell(f, ::Grid1D{PeriodicBC}, inputs)
    f(((cells[end], cells[1], cells[2]) for cells in inputs)..., 1)
    f(((cells[end-1], cells[end], cells[1]) for cells in inputs)..., lastindex(inputs[1]))
end


function create_wall_ghost_cell(Q)
    @SVector [Q[1], -Q[2]]
end

function for_each_boundary_cell(f, ::Grid1D{WallBC}, inputs)
    f(((create_wall_ghost_cell(cells[1]), cells[1], cells[2]) for cells in inputs)..., 1)
    f(((cells[end-1], cells[end], create_wall_ghost_cell(cells[end])) for cells in inputs)..., lastindex(inputs[1]))
end


function get_dx(grid::UniformGrid1D, cell_idx)
    grid.dx
end

function separate_variables(U::AbstractMatrix{SVector{N, T}}) where {N, T}
    stacked = stack(U)
    (stacked[i,:,:]' for i in 1:N)
end


cells(grid::UniformGrid1D) = grid.cells

inner_cells(grid::UniformGrid1D, cells=cells(grid)) = cells
