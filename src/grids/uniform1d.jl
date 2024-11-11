using QuadGK, StaticArrays

export UniformGrid1D, for_each_cell, for_each_inner_cell, get_neighbour_indices, get_dx, inner_cells


abstract type Grid1D{BC <: BoundaryCondition} <: Grid end

struct UniformGrid1D{BC <: BoundaryCondition, Float <: AbstractFloat, Int <: Integer, Cell} <: Grid1D{BC}
    dx::Float
    bc::BC
    cells::Vector{Cell}
    gc::Int

    function UniformGrid1D(N, bc::BoundaryCondition, u0, domain, ghost_cells=1)
        x_L, x_R = domain
        dx = (x_R - x_L) / (N+1)
        x = -ghost_cells * dx + x_L:dx:x_R + ghost_cells*dx

        averages = _calc_average.(u0, x[1:end-1], x[2:end])

        if eltype(averages) <: Number
            return new{typeof(bc), typeof(dx), typeof(ghost_cells), eltype(averages)}(
                dx, bc, averages, ghost_cells
            )
        else
            conserved_variables::Int64 = length(averages[1])
            FloatType = eltype(averages[1])
            cells = [SVector{conserved_variables, FloatType}(average) for average in averages]
            return new{typeof(bc), typeof(dx), typeof(ghost_cells), eltype(cells)}(
                dx, bc, cells, ghost_cells
            )
        end
    end
end


function _calc_average(f, a, b)
    quadgk(f, a, b)[1] / (b - a)
end

function reduce_cells(f, grid::Grid1D, initial_value)
    inner_cells_indices = eachindex(cells(grid))[grid.gc+1:end-grid.gc]
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


function for_each_inner_cell(f, grid::Grid1D; gc=grid.gc, p=gc, q=p)
    for cell_idx in eachindex(cells(grid))[gc+1:end-gc]
        ileft = cell_idx - p
        iright = cell_idx + q
        f(cells(grid), ileft, cell_idx, iright)
    end
end

update_bc!(grid::Grid1D, eq::Equation) = update_bc!(cells(grid), grid, eq)

function update_bc!(out, grid::Grid1D{NeumannBC}, ::Equation)
    for i in 1:grid.gc
        out[i] = out[grid.gc+1]
    end
    for i in lastindex(out)-grid.gc+1:lastindex(out)
        out[i] = out[end-grid.gc]
    end
end

@views function update_bc!(out, grid::Grid1D{PeriodicBC}, ::Equation)
    out[1:grid.gc] .= out[end-2*grid.gc+1:end-grid.gc]
    out[end-grid.gc+1:end] .= out[grid.gc+1:2*grid.gc]
end

@views function update_bc!(out, grid::Grid1D{WallBC}, ::ShallowWater1D)
    for i in 1:grid.gc
        h, uh = out[grid.gc + i]
        out[grid.gc-i+1] = @SVector [h, -uh]

        h, uh = out[end-grid.gc-i+1]
        out[end-grid.gc+i] = @SVector [h, -uh]
    end
end


function get_dx(grid::UniformGrid1D, cell_idx)
    grid.dx
end


cells(grid::UniformGrid1D) = grid.cells

inner_cells(grid::UniformGrid1D, cells=cells(grid)) = cells[grid.gc+1:end-grid.gc]
