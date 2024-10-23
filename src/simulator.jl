export Simulator, simulate!

struct Simulator{SystemType <: System, GridType <: Grid, TimeStepperType <: TimeStepper, Float <: AbstractFloat, Cells}
    system::SystemType
    grid::GridType
    timestepper::TimeStepperType
    substep_buffers::Vector{Cells}
    t::Base.RefValue{Float}
    dt::Float

    function Simulator(system::System, grid::Grid, timestepper::TimeStepper, t::Float=zero(Float)) where Float <: AbstractFloat
        substep_buffers = create_substep_buffers(grid, number_of_substeps(timestepper))
        new{typeof(system), typeof(grid), typeof(timestepper), Float, eltype(substep_buffers)}(system, grid, timestepper, substep_buffers, Base.RefValue(t), zero(Float))
    end
end



@views function add_time_derivative!(output, eq::Equation, F::NumericalFlux, grid::Grid, dt)
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


function calc_max_dt(eq::Equation, grid::Grid, max_dt)
    reduce_cells(grid, max_dt) do acc_max_dt, cell, dx
        speed = compute_max_abs_speed(eq, cell)
        if !iszero(speed)
            min(dx/speed, acc_max_dt)
        else
            acc_max_dt
        end
    end
end


function perform_step!(simulator::Simulator, max_dt)
    dt = calc_max_dt(simulator.system.eq, simulator.grid, max_dt)

    integrate!(simulator.grid, simulator.substep_buffers, dt, simulator.timestepper) do substep_output, dt
        add_time_derivative!(substep_output,
                             simulator.system.eq,
                             simulator.system.numerical_flux,
                             simulator.system.grid,
                             dt)
    end

    simulator.t[] += dt
end


function simulate!(simulator::Simulator, T, max_dt, callbacks::Vector{Callback}=[]) where Callback
    while simulator.t[] < T
        perform_step!(simulator, max_dt)

        for callback in callbacks
            callback(simulator)
        end
    end
end


