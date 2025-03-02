using ElasticArrays

export Simulator, simulate!, simulate_and_aggregate!

"""
    Simulator{GridType <: Grid, Float <: AbstractFloat}

A type representing a simulation of a conservation law on a grid. The simulation is described by a `ConservedSystem` and a `GridType`.
"""
struct Simulator{GridType <: Grid, Float <: AbstractFloat}
    system::ConservedSystem
    grid::GridType
    t::Base.RefValue{Float}
    function Simulator(
            system::ConservedSystem,
            grid::Grid,
            t::Float=zero(Float)
        ) where Float <: AbstractFloat

        new{typeof(grid), Float}(
            system,
            grid,
            Base.RefValue(t)
        )
    end
end


@views function set_time_derivative!(out, grid::Grid, left, right, eq::Equation, F::NumericalFlux, dt)
    p, q = stencil_size(F)
    for_each_interior_cell(grid; p=p, q=q) do cells, ileft, imiddle, iright
        left_minus = (right[ileft],)
        right_minus = (left[imiddle],)
        left_plus = (right[imiddle],)
        right_plus = (left[iright],)

        dx = get_dx(grid, imiddle)

        flux_minus = F(eq, left_minus, right_minus, dx, dt)
        flux_plus = F(eq, left_plus, right_plus, dx, dt)

        out[imiddle] = (flux_minus - flux_plus) / dx
    end

    for_each_boundary_cell(grid, left, right) do (LL, ML, RL), (LR, MR, RR), idx
        dx = get_dx(grid, idx)

        flux_minus = F(eq, (LR,), (ML,), dx, dt)
        flux_plus = F(eq, (MR,), (RL,), dx, dt)
        
        out[idx] = (flux_minus - flux_plus) / dx
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
    dt = integrate!(simulator.grid, simulator.system.timestepper, simulator.system,
               (eq, grid, left, right) -> calc_max_dt(eq, grid, left, right, max_dt))

    simulator.t[] += Real(dt) # May be a dual number
end


"""
    simulate!(simulator::Simulator, T, max_dt, callbacks::Vector{Callback}=[]) where Callback

Simulate the system described by `simulator` until time `T` using a time step size of at most `max_dt`. The `callbacks` are called at each time step.
The cell averages are mutated in place.

# Example
```julia
eq = BurgersEQ()
bc = NeumannBC()

N = 20
x_L, x_R = 0.0, 1.0
x = cell_centers(N, x_L, x_R)
u0(x) = x .< 0.5 ? 1.0 : 0.0
U0 = u0.(x)
grid = UniformGrid1D(N, bc, U0, (x_L, x_R))

F = LaxFriedrichsFlux()
reconstruction = NoReconstruction()
timestepper = ForwardEuler(grid)

system = ConservedSystem(eq, reconstruction, F, timestepper)
simulator = Simulator(system, grid, 0.)

simulate!(simulator, 1., 0.1)

U = cells(grid)
```
"""
function simulate!(simulator::Simulator, T, max_dt, callbacks::Vector{Callback}=[]) where Callback
    while simulator.t[] < T
        perform_step!(simulator, max_dt)

        for callback in callbacks
            callback(simulator)
        end
    end
end


"""
    simulate_and_aggregate!(simulator::Simulator, T, max_dt, callbacks::Vector{Callback}=[]) where Callback

The same as [`simulate!`](@ref) but keeps track of and returns the cell averages and time at each time step.

# Example
```julia
eq = BurgersEQ()
bc = NeumannBC()

N = 20
x_L, x_R = 0.0, 1.0
x = cell_centers(N, x_L, x_R)
u0(x) = x .< 0.5 ? 1.0 : 0.0
U0 = u0.(x)
grid = UniformGrid1D(N, bc, U0, (x_L, x_R))

F = LaxFriedrichsFlux()
reconstruction = NoReconstruction()
timestepper = ForwardEuler(grid)

system = ConservedSystem(eq, reconstruction, F, timestepper)
simulator = Simulator(system, grid, 0.)

U, t = simulate_and_aggregate!(simulator, 1., 0.1)

u(x, t) = u0(x-0.5t)
Viz.animate_solution(U', u, grid, t, 2)
```
"""
function simulate_and_aggregate!(simulator::Simulator, T, max_dt, callbacks::Vector{Callback}=[]) where Callback
    U = ElasticMatrix(reshape(cells(simulator.grid), :, 1))
    t = ElasticVector([simulator.t[]])

    function collect_state(simulator)
        append!(U, copy(cells(simulator.grid)))
        append!(t, simulator.t[])
    end

    simulate!(simulator, T, max_dt, [collect_state, callbacks...])

    return U, t
end