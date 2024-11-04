using QuadGK, StaticArrays

export UniformGrid1D, for_each_cell, get_neighbour_indices, get_dx


abstract type Grid1D{BC <: BoundaryCondition} <: Grid end

struct UniformGrid1D{BC <: BoundaryCondition, Float <: AbstractFloat, Cell} <: Grid1D{BC}
    dx::Float
    bc::BC
    cells::Vector{Cell}

    function UniformGrid1D(N, bc::BoundaryCondition, u0, domain)
        x_L, x_R = domain
        dx = (x_R - x_L) / (N+1)
        x = x_L:dx:x_R

        averages = _calc_average.(u0, x[1:end-1], x[2:end])

        if eltype(averages) <: Number
            return new{typeof(bc), typeof(dx), eltype(averages)}(dx, bc, averages)
        else
            conserved_variables::Int64 = length(averages[1])
            FloatType = eltype(averages[1])
            cells = [SVector{conserved_variables, FloatType}(averages[j]) for j in 1:N+1]
            return new{typeof(bc), typeof(dx), eltype(cells)}(dx, bc, cells)
        end
    end
end


function _calc_average(f, a, b)
    quadgk(f, a, b)[1] / (b - a)
end

function reduce_cells(f, grid::Grid1D, initial_value)
    reduce(eachindex(cells(grid)), init=initial_value) do acc, i
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


function get_neighbour_indices(grid::Grid1D{PeriodicBC}, cell_idx, p, q=p)
    N = length(cells(grid))
    DeferredRange(-p:q) do offset
        mod1(cell_idx + offset, N)
    end
end

function get_neighbour_indices(grid::Grid1D{NeumannBC}, cell_idx, p, q=p)
    N = length(cells(grid))
    DeferredRange(-p:q) do offset
        clamp(cell_idx + offset, 1, N)
    end
end

function get_dx(grid::UniformGrid1D, cell_idx)
    grid.dx
end


function cells(grid::UniformGrid1D)
    grid.cells
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
    DeferredRange(f, range::UnitRange) = DeferredRange(f, first(range), last(range))
end

Base.size(r::DeferredRange) = (r.stop - r.start + one(r.start),)

function Base.getindex(r::DeferredRange, i::T) where {T <: Integer}
    r.f(r.start + i - one(T))
end
