export Simulator, simulate!

struct Simulator{SystemType <: System, GridType <: Grid, TimeStepperType <: TimeStepper, Float <: AbstractFloat}
    system::SystemType
    grid::GridType
    timestepper::TimeStepperType
    substep_buffers::Vector
    t::Ref{Float}
    dt::Float

    function Simulator(system::System, grid::Grid, timestepper::TimeStepper, t::Float=zero(Float)) where Float <: AbstractFloat
        substep_buffers = create_substep_buffers(grid, number_of_substeps(timestepper))
        new{typeof(system), typeof(grid), typeof(timestepper), Float}(system, grid, timestepper, substep_buffers, Ref(t), zero(Float))
    end
end


# function integrate!(output, grid::Grid, dt)
#     for_each_cell(grid) do cells, dx, cell_idx
#         output[cell_idx] = dt * output[cell_idx] + cells[cell_idx]
#     end
# end


function add_time_derivative!(output, eq::Equation, F::NumericalFlux, grid::Grid, dt)
    p = number_of_cells(F)
    for_each_cell(grid, p) do cells, dx, cell_idx, neighbour_cell_idx
        left_minus = cells[neighbour_cell_idx[1:p]]
        right_minus = cells[neighbour_cell_idx[end-p:end-1]]
        left_plus = cells[neighbour_cell_idx[2:p+1]]
        right_plus = cells[neighbour_cell_idx[end-p+1:end]]

        flux_minus = F(eq, left_minus, right_minus, dx, dt)
        flux_plus = F(eq, left_plus, right_plus, dx, dt)

        output[cell_idx] += (flux_minus - flux_plus) / dx
    end
end

function set_zero!(output, grid::Grid)
    for_each_cell(grid) do cells, dx, cell_idx
        output[cell_idx] = zero(eltype(output))
    end
end

# function set_new_state!(grid::Grid, new_states)
#     for_each_cell(grid) do cells, dx, cell_idx
#         cells[cell_idx] = new_states[cell_idx]
#     end
# end

function calc_max_dt(eq::Equation, grid::Grid, max_dt)
    for_each_cell(grid) do cells, dx, cell_idx
        speed = compute_max_abs_speed(eq, cells[cell_idx])
        if !iszero(speed)
            max_dt = min(dx/speed, max_dt)
        end
    end
    return max_dt
end


function perform_step!(simulator::Simulator, max_dt)
    dt = calc_max_dt(simulator.system.eq, simulator.grid, max_dt)

    set_zero!.(simulator.substep_buffers, Ref(simulator.grid))

    integrate!(simulator.grid, simulator.substep_buffers, dt, simulator.timestepper) do substep_output, dt
        add_time_derivative!(substep_output,
                             simulator.system.eq,
                             simulator.system.numerical_flux,
                             simulator.system.grid,
                             dt)
    end

    simulator.t[] += dt
end


function simulate!(simulator::Simulator, T, max_dt, callbacks)
    # for callback in callbacks
    #     callback(simulator)
    # end
    while simulator.t[] < T
        perform_step!(simulator, max_dt)

        for callback in callbacks
            callback(simulator)
        end
    end
end


