export Simulator, simulate!

struct Simulator{
        SystemType <: System,
        GridType <: Grid,
        Float <: AbstractFloat
    }
    system::SystemType
    grid::GridType
    t::Base.RefValue{Float}
    function Simulator(
            system::System,
            grid::Grid,
            t::Float=zero(Float)
        ) where Float <: AbstractFloat

        new{typeof(system), typeof(grid), Float}(
            system,
            grid,
            Base.RefValue(t)
        )
    end
end


@views function set_time_derivative!(out, grid::Grid, left, right, eq::Equation, F::NumericalFlux, dt)
    p = number_of_cells(F)
    for_each_inner_cell(grid, p) do cells, ileft, imiddle, iright
        left_minus = right[ileft:imiddle-1]
        right_minus = left[imiddle:iright-1]
        left_plus = right[ileft+1:imiddle]
        right_plus = left[imiddle+1:iright]

        dx = get_dx(grid, imiddle)

        flux_minus = F(eq, left_minus, right_minus, dx, dt)
        flux_plus = F(eq, left_plus, right_plus, dx, dt)

        out[imiddle] = (flux_minus - flux_plus) / dx
    end
end


function calc_max_dt(eq::Equation, grid::Grid, left, right, max_dt)
    reduce_cells(grid, max_dt) do acc_max_dt, cells, idx
        dx = get_dx(grid, idx)
        speed_left = compute_max_abs_eigenvalue(eq, left[idx])
        speed_right = compute_max_abs_eigenvalue(eq, right[idx])
        max_abs_speed = max(speed_left, speed_right)
        if !iszero(max_abs_speed)
            min(dx/max_abs_speed, acc_max_dt)
        else
            acc_max_dt
        end
    end
end


function perform_step!(simulator::Simulator, max_dt)
    dt = integrate!(simulator.grid, simulator.system,
               (eq, grid, left, right) -> calc_max_dt(eq, grid, left, right, max_dt))

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
