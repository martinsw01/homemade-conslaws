export RK2

struct RK2 <: TimeStepper end

number_of_substeps(::RK2) = 2

function integrate!(add_time_derivative!, grid::Grid, substep_outputs, dt, ::RK2)
    for_each_cell(grid) do cells, dx, cell_idx
        substep_outputs[1][cell_idx] = cells[cell_idx]
        substep_outputs[2][cell_idx] = zero(eltype(cells))
    end

    add_time_derivative!(substep_outputs[2], dt)

    for_each_cell(grid) do cells, dx, cell_idx
        cells[cell_idx] = cells[cell_idx] + 0.5dt * substep_outputs[2][cell_idx]
        substep_outputs[2][cell_idx] = zero(eltype(cells))
    end

    add_time_derivative!(substep_outputs[2], dt)

    for_each_cell(grid) do cells, dx, cell_idx
        cells[cell_idx] = substep_outputs[1][cell_idx] + dt*substep_outputs[2][cell_idx]
    end
    
end