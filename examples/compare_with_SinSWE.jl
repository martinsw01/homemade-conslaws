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
    simulator = SinSWE.Simulator(backend, conserved_system, timestepper, grid, cfl=1.)
    
    x = SinSWE.cell_centers(grid)
    initial = u0.(x)
    SinSWE.set_current_state!(simulator, initial)

    H = ElasticMatrix(reshape([Q[1] for Q in initial], :, 1))
    UH = ElasticMatrix(reshape([Q[2] for Q in initial], :, 1))
    t = ElasticVector([0.])

    function collect_state(t_j, simulator)
        append!(UH, [Q[2] for Q in SinSWE.current_interior_state(simulator)])
        append!(H, [Q[1] for Q in SinSWE.current_interior_state(simulator)])
        append!(t, t_j)

    end

    
    SinSWE.simulate_to_time(simulator, T; callback=collect_state)
    
    return H', UH', t, SinSWE.cell_centers(grid)
end

function run_homemade_conslaws_simulation(T, N, u0)
    eq = ShallowWater1D(1.)
    bc = WallBC()
    grid = UniformGrid1D(N, bc, u0, (0, 1), 2)
    x = -1.5grid.dx:grid.dx:1+1.5grid.dx
    grid.cells .= u0.(x)

    dt = grid.dx

    F = CentralUpwind()

    reconstruction = LinearReconstruction(grid)
    timestepper = RK2(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    H = ElasticMatrix(reshape([Q[1] for Q in inner_cells(grid)], :, 1))
    UH = ElasticMatrix(reshape([Q[2] for Q in inner_cells(grid)], :, 1))
    t = ElasticVector([0.])

    function collect_state(simulator)
        append!(UH, [Q[2] for Q in inner_cells(grid)])
        append!(H, [Q[1] for Q in inner_cells(grid)])
        append!(t, simulator.t[])
    end

    simulate!(simulator, T, dt, [collect_state])

    return H', UH', t, x[2:end-1]
end

u0(x) = @SVector[exp.(-(x - 0.5)^2 / 0.001) .+ 1.5, 0.0 .* x]

N = 50
T = 2
H_SinSWE, UH_SinSwe, t_SinSWE, x_SinSWE = run_SinSWE_simulation(T, N, u0);
H, UH, t, x = run_homemade_conslaws_simulation(T, N-1, u0);

homemade_conslaws.Viz.animate_solution(
    (H_SinSWE, H), ["H_SinSWE" "H"], x_SinSWE, t, 4)
