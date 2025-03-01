export RK2

"""
    RK2{Cells}

A second order Runge-Kutta time-stepping scheme, a strong stability preserving version of the standard RK2 scheme,
applying two steps of the forward Euler scheme. Preallocates two buffers to store intermediate values.
"""
struct RK2{Cells} <: TimeStepper
    substep_buffers::Vector{Cells}

    function RK2(grid::Grid)
        substep_buffers = [create_buffer(grid), create_buffer(grid)]
        new{eltype(substep_buffers)}(substep_buffers)
    end
end

function integrate!(grid::Grid, timestepper::RK2, system, compute_max_dt)
    equation = system.eq
    F = system.numerical_flux
    reconstruction = system.reconstruction


    left, right = reconstruct(reconstruction, grid)

    dt = compute_max_dt(equation, grid, left, right)

    set_time_derivative!(timestepper.substep_buffers[1], grid, left, right, equation, F, dt)

    for_each_cell(grid) do cells, i
        timestepper.substep_buffers[1][i] = cells[i] + dt * timestepper.substep_buffers[1][i]
    end

    left, right = reconstruct(reconstruction, grid, timestepper.substep_buffers[1])

    set_time_derivative!(timestepper.substep_buffers[2], grid, left, right, equation, F, dt)

    for_each_cell(grid) do cells, i
        timestepper.substep_buffers[2][i] = timestepper.substep_buffers[1][i] + dt * timestepper.substep_buffers[2][i]
        cells[i] = 0.5 * (cells[i] + timestepper.substep_buffers[2][i])
    end

    return dt
end