# using Revise
using StaticArrays
using ElasticArrays
using SinSWE
using homemade_conslaws

function run_SinSWE_simulation(T, N, u0)
    backend = make_cpu_backend()

    grid = SinSWE.CartesianGrid(N; gc=2, boundary=SinSWE.WallBC())
    
    equation = SinSWE.ShallowWaterEquations1DPure(1., 1.)
    reconstruction = SinSWE.LinearReconstruction(1.)
    numericalflux = SinSWE.CentralUpwind(equation)
    timestepper = SinSWE.RungeKutta2()
    conserved_system = SinSWE.ConservedSystem(backend, reconstruction, numericalflux, equation, grid, [])
    simulator = SinSWE.Simulator(backend, conserved_system, timestepper, grid, cfl=1.) # Should be cfl=0.2, but it is not yet implemented in homemade_conslaws
    
    x = SinSWE.cell_centers(grid)
    initial = u0.(x)
    SinSWE.set_current_state!(simulator, initial)

    H = ElasticMatrix(reshape([Q[1] for Q in initial], :, 1))
    UH = ElasticMatrix(reshape([Q[2] for Q in initial], :, 1))
    t = ElasticVector([0.])

    function collect_state(t_j, simulator)
        UH_next = [Q[2] for Q in SinSWE.current_interior_state(simulator)]
        H_next = [Q[1] for Q in SinSWE.current_interior_state(simulator)]

        append!(UH, UH_next)
        append!(H, H_next)
        append!(t, t_j)
    end

    
    SinSWE.simulate_to_time(simulator, T; callback=collect_state)
    
    return H', UH', t, SinSWE.cell_centers(grid)
end

function run_homemade_conslaws_simulation(T, N, u0)
    eq = ShallowWater1D(1.)
    bc = WallBC()
    grid = UniformGrid1D(N, bc, u0.(cell_centers(N, 0, 1)), (0, 1))
    x = cell_centers(grid)
    grid.cells .= u0.(x)

    dt = grid.dx

    F = CentralUpwind()

    reconstruction = LinearReconstruction(grid)
    timestepper = RK2(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    Q, t = simulate_and_aggregate!(simulator, T, dt)
    H, UH = separate_variables(Q)

    return H, UH, t, x
end

u0(x) = @SVector[exp.(-(x - 0.5)^2 / 0.001) .+ 1.5, 0.0 .* x]

N = 50
T = 2.
H_SinSWE, UH_SinSwe, t_SinSWE, x_SinSWE = run_SinSWE_simulation(T, N, u0)
H, UH, t, x = run_homemade_conslaws_simulation(T, N-1, u0)

max_abs_difference = maximum(abs.(H_SinSWE .- H)[1:end-1,:])
@show max_abs_difference

homemade_conslaws.Viz.animate_solution(
    H_SinSWE, "H_SinSWE", x_SinSWE, t_SinSWE, T)
