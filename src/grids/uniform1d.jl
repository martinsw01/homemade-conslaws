using QuadGK

export UniformGrid1D

export UniformGrid1D, for_each_cell

struct UniformGrid1D{BC <: BoundaryCondition} <: Grid
    dx::Float64
    bc::BC
    domain::Tuple{Float64, Float64}
    cells::Vector{Float64}

    function UniformGrid1D(N, bc::BC, u0, domain) where BC <: BoundaryCondition
        x_L, x_R = domain
        dx = (x_R - x_L) / (N+1)
        x = x_L:dx:x_R
        cells = _integrate.(u0, x[1:end-1], x[2:end]) ./ dx

        new{BC}(dx, bc, domain, cells)
    end
end


function _integrate(f, a, b)
    quadgk(f, a, b)[1]
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
    f.(Ref(grid.cells), grid.dx, eachindex(grid.cells))
end


function for_each_cell(f, grid::UniformGrid1D, compute_neighbouring_indices::Function)
    for cell_idx in eachindex(grid.cells)
        f(grid.cells, grid.dx, cell_idx, compute_neighbouring_indices(grid.cells, cell_idx))
    end
end


function for_each_cell(f, grid::UniformGrid1D{PeriodicBC}, p::Integer)
    N_pluss_1 = length(grid.cells)
    for_each_cell(f, grid, 
                  (cells, cell_idx) -> mod1.(cell_idx .+ (-p:p), N_pluss_1),
                #   cell_idx -> mod1.(cell_idx - p:cell_idx - 1, N),
                #   cell_idx -> mod1.(cell_idx + 1:cell_idx + p, N)
                  )
end


function for_each_cell(f, grid::UniformGrid1D{NeumannBC}, p::Integer)
    N_pluss_1 = length(grid.cells)
    for_each_cell(f, grid,
                  (cells, cell_idx) -> clamp.(cell_idx .+ (-p:p), 1, N_pluss_1)
                #   cell_idx -> max.(1, cell_idx - p:cell_idx - 1),
                #   cell_idx -> min.(N, cell_idx + 1:cell_idx + p)
                  )
end
