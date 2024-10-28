export NoReconstruction

struct NoReconstruction{Cells} <: Reconstruction
    reconstruction_buffer::Cells

    function NoReconstruction(grid::Grid)
        buffer = create_buffer(grid)
        new{typeof(buffer)}(buffer)
    end
end

function reconstruct(r::NoReconstruction, grid::Grid, cell_averages=cells(grid))
    for_each_cell(grid) do cells, cell_idx
        r.reconstruction_buffer[cell_idx] = cell_averages[cell_idx]
    end
    return r.reconstruction_buffer, r.reconstruction_buffer
end