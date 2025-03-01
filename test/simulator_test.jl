using Test
using homemade_conslaws

@testset "Test riemann problem single step" begin
    bc = NeumannBC()
    u0 = [1., 0.5, 0.]
    N = 3
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))
    dt = grid.dx
    T = dt
    eq = BurgersEQ()
    F = LaxFriedrichsFlux()
    reconstruction = NoReconstruction()
    timestepper = ForwardEuler(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    simulate!(simulator, T, dt)

    expected_cells = [15/16, 3/4, 5/16]

    @test cells(grid) ≈ expected_cells
end


@testset "Test cfl condition" begin
    bc = NeumannBC()
    u0 = [1., 0.5, 0.]
    N = 3
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))
    max_dt = 2*grid.dx
    T = max_dt
    eq = BurgersEQ()
    F = LaxFriedrichsFlux()
    reconstruction = NoReconstruction()
    timestepper = ForwardEuler(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    simulate!(simulator, T, max_dt)

    expected_cells = [297 / 320, 5 / 6, 629 / 960]

    @test cells(grid) ≈ expected_cells
end

@testset "Test riemann problem 2 steps" begin
    bc = NeumannBC()
    u0 = [1., 0.5, 0.]
    N = 3
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))
    max_dt = grid.dx
    T = 2*max_dt
    eq = BurgersEQ()
    F = LaxFriedrichsFlux()
    reconstruction = NoReconstruction()
    timestepper = ForwardEuler(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    simulate!(simulator, T, max_dt)

    expected_cells = [945 / 1024, 105 / 128, 663 / 1024]

    @test cells(grid) ≈ expected_cells
end

@testset "Test accumulation of states" begin
    bc = NeumannBC()
    u0 = [1., 0.5, 0.]
    N = 3
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))
    max_dt = grid.dx
    T = 2*max_dt
    eq = BurgersEQ()
    F = LaxFriedrichsFlux()
    reconstruction = NoReconstruction()
    timestepper = ForwardEuler(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    U, t = simulate_and_aggregate!(simulator, T, max_dt)

    expected_cells = [1. 0.5 0.;
                      15/16 3/4 5/16;
                      945/1024 105/128 663/1024]'
    expected_t = [0., max_dt, T]

    @test U ≈ expected_cells
    @test t ≈ expected_t
end