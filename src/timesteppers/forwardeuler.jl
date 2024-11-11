export ForwardEuler

struct ForwardEuler{Cells} <: TimeStepper
    substep_buffer::Cells

    function ForwardEuler(grid::Grid)
        substep_buffer = create_buffer(grid)
        new{typeof(substep_buffer)}(substep_buffer)
    end
end


function integrate!(grid::Grid, system::ConservedSystem{E, R, NF, ForwardEuler{Cells}}, compute_max_dt) where {E, R, NF, Cells}
    equation = system.eq
    timestepper = system.timestepper
    reconstruction = system.reconstruction
    F = system.numerical_flux
    
    left, right = reconstruct(reconstruction, grid)

    dt = compute_max_dt(equation, grid, left, right)

    set_time_derivative!(timestepper.substep_buffer, grid, left, right, equation, F, dt)

    for_each_cell(grid) do cells, i
        cells[i] = cells[i] + dt * timestepper.substep_buffer[i]
    end

    update_bc!(grid, equation)

    return dt
end