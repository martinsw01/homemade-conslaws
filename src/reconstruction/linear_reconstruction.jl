export LinearReconstruction

struct LinearReconstruction{Cells} <: Reconstruction
    left_buffer::Cells
    right_buffer::Cells

    function LinearReconstruction(grid::Grid)
        left = create_buffer(grid)
        right = create_buffer(grid)
        new{typeof(left)}(left, right)
    end
end


function minmod(a, b)
    if a * b > 0
        return sign(a) * min(abs(a), abs(b))
    else
        return 0
    end
end


function reconstruct(reconstruction::LinearReconstruction, grid::Grid, cell_averages=cells(grid))
    for_each_inner_cell(grid, 1) do cells, ileft, imiddle, iright
        left = cell_averages[ileft]
        center = cell_averages[imiddle]
        right = cell_averages[iright]

        downwind = right - center
        upwind = center - left
        slope = minmod(downwind, upwind)

        reconstruction.left_buffer[imiddle] = center - 0.5 * slope
        reconstruction.right_buffer[imiddle] = center + 0.5 * slope
    end

    return reconstruction.left_buffer, reconstruction.right_buffer
end