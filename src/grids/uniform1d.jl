using QuadGK

export UniformGrid1D

export UniformGrid1D, for_each_cell

struct UniformGrid1D{BC <: BoundaryCondition, Float <: AbstractFloat} <: Grid
    dx::Float64
    bc::BC
    domain::Tuple{Float, Float}
    cells::Vector{Float}

    function UniformGrid1D(N, bc::BoundaryCondition, u0, domain)
        x_L, x_R = domain
        dx = (x_R - x_L) / (N+1)
        x = x_L:dx:x_R
        cells = _calc_average.(u0, x[1:end-1], x[2:end])

        new{typeof(bc), typeof(dx)}(dx, bc, domain, cells)
    end
end


function _calc_average(f, a, b)
    quadgk(f, a, b)[1] / (b - a)
end


function reduce_cells(f, grid::UniformGrid1D{BoundaryCondition, Float}, initial_value::Float) where {BoundaryCondition, Float}
    reduce(grid.cells, init=initial_value) do acc, cell
        f(acc, cell, grid.dx)
    end
end
    


function create_substep_buffers(grid::UniformGrid1D, num_substeps)
    [grid.cells[:] for _ in 1:num_substeps]
end
    


"""
    for_each_cell(f, grid)

Apply the function `f` to each cell in the grid, their width `dx` and the index of the cell.
"""
function for_each_cell end


function for_each_cell(f, grid::UniformGrid1D)
    for cell_idx in eachindex(grid.cells)
        f(grid.cells, grid.dx, cell_idx)
    end
end



function for_each_cell(f, grid::UniformGrid1D, compute_neighbouring_indices::Function)
    for cell_idx in eachindex(grid.cells)
        f(grid.cells, grid.dx, cell_idx, compute_neighbouring_indices(grid.cells, cell_idx))
    end
end


function for_each_cell(f, grid::UniformGrid1D{PeriodicBC, Float}, p::Integer) where {Float <: AbstractFloat}
    N_pluss_1 = length(grid.cells)
    for_each_cell(f, grid, 
                  (cells, cell_idx) -> periodic_indices(N_pluss_1, cell_idx, p),
                  )
end

periodic_indices(N, cell_idx, p) = DeferredRange(-p:p) do offset
    mod1(cell_idx + offset, N)
end


function for_each_cell(f, grid::UniformGrid1D{NeumannBC, Float}, p::Integer) where {Float <: AbstractFloat}
    N_pluss_1 = length(grid.cells)
    for_each_cell(f, grid,
                  (cells, cell_idx) -> neumann_indices(N_pluss_1, cell_idx, p)
                  )
end


neumann_indices(N, cell_idx, p) = DeferredRange(-p:p) do offset
    clamp(cell_idx + offset, 1, N)
end


"""
    DeferredRange(f, start, stop)

Similar to a `UnitRange(start, stop)`, but the values are computed on each access using the function `f`. Allows
for allocation free views of arrays similar to `UnitRange` but with deferred computation of the values.
"""
struct DeferredRange{T <: Integer, F} <: AbstractArray{T, 1}
    f::F
    start::T
    stop::T
    function DeferredRange(f, start::T, stop::T) where {T <: Integer}
        new{T, typeof(f)}(f, start, stop)
    end
    function DeferredRange(f, range::UnitRange)
        new{eltype(range), typeof(f)}(f, first(range), last(range))
    end
end

Base.size(r::DeferredRange) = (r.stop - r.start + one(r.start),)

function Base.getindex(r::DeferredRange, i::T) where {T <: Integer}
    r.f(r.start + i - one(T))
end
