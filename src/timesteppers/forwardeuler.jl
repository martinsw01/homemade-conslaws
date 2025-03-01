export ForwardEuler

"""
    ForwardEuler{Cells}

The standard first order TVD forward Euler time-stepping scheme. Preallocates a buffer to store the time derivatives.
"""
struct ForwardEuler{Cells} <: TimeStepper
    substep_buffer::Cells

    function ForwardEuler(grid::Grid)
        substep_buffer = create_buffer(grid)
        new{typeof(substep_buffer)}(substep_buffer)
    end
end


function integrate!(grid::Grid, timestepper::ForwardEuler, system, compute_max_dt)
    equation = system.eq
    reconstruction = system.reconstruction
    F = system.numerical_flux
    
    left, right = reconstruct(reconstruction, grid)

    dt = compute_max_dt(equation, grid, left, right)

    set_time_derivative!(timestepper.substep_buffer, grid, left, right, equation, F, dt)

    for_each_cell(grid) do cells, i
        cells[i] = cells[i] + dt * timestepper.substep_buffer[i]
    end

    return dt
end