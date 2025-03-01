export LinearReconstruction


"""
    LinearReconstruction{Cells}

Second order reconstruction using the minmod slope. Preallocates buffers for the left and right reconstructed values.
"""
struct LinearReconstruction{Cells} <: Reconstruction
    left_buffer::Cells
    right_buffer::Cells

    function LinearReconstruction(grid::Grid)
        left = create_buffer(grid)
        right = create_buffer(grid)
        new{typeof(left)}(left, right)
    end
end


function minmod(a, b, c)
    if sign(a) == sign(b) == sign(c)
        return sign(a) * min(abs(a), abs(b), abs(c))
    else
        return 0
    end
end


function _reconstruct_cell!(left_buffer, right_buffer, left, center, right, idx)
    downwind = right .- center
    upwind = center .- left
    central = 0.5 .* (downwind .+ upwind)

    slope = minmod.(downwind, central, upwind)

    left_buffer[idx] = center .- 0.5 .* slope
    right_buffer[idx] = center .+ 0.5 .* slope
end


function reconstruct(reconstruction::LinearReconstruction, grid::Grid, cell_averages=cells(grid))
    for_each_interior_cell(grid) do cells, ileft, imiddle, iright
        left = cell_averages[ileft]
        center = cell_averages[imiddle]
        right = cell_averages[iright]

        _reconstruct_cell!(reconstruction.left_buffer, reconstruction.right_buffer, left, center, right, imiddle)
    end

    for_each_boundary_cell(grid, cell_averages) do stencil, idx
        _reconstruct_cell!(reconstruction.left_buffer, reconstruction.right_buffer, stencil..., idx)
    end

    return reconstruction.left_buffer, reconstruction.right_buffer
end