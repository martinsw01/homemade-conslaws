export ForwardEuler

struct ForwardEuler <: TimeStepper end

number_of_substeps(::ForwardEuler) = 1


function integrate!(add_time_derivative!, grid::Grid, substep_outputs, dt, ::ForwardEuler)
    add_time_derivative!(substep_outputs[1], dt)

    for_each_cell(grid) do cells, dx, cell_idx
        cells[cell_idx] = cells[cell_idx] + dt * substep_outputs[1][cell_idx]
    end
end
    