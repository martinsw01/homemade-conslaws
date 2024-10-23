using Test
using homemade_conslaws

@testset "Test riemann problem single step" begin
    bc = NeumannBC()
    u0(x) = x < 0 ? 1 : 0
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))
    dt = grid.dx
    T = dt
    eq = BurgersEQ()
    F = LaxFriedrichsFlux()
    system = ConservedSystem(eq, grid, F)
    timestepper = ForwardEuler()
    simulator = Simulator(system, grid, timestepper, 0.)

    @test grid.cells ≈ [1., 0.5, 0.] #rtol=1e-16

    grid.cells[:] = [1., 0.5, 0.]

    simulate!(simulator, T, dt)

    expected_cells = [15/16, 3/4, 5/16]

    @test grid.cells ≈ expected_cells
end


@testset "Test cfl condition" begin
    bc = NeumannBC()
    u0(x) = x < 0 ? 1 : 0
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))
    grid.cells[:] = [1., 0.5, 0.]
    max_dt = 2*grid.dx
    T = max_dt
    eq = BurgersEQ()
    F = LaxFriedrichsFlux()
    system = ConservedSystem(eq, grid, F)
    timestepper = ForwardEuler()
    simulator = Simulator(system, grid, timestepper, 0.)

    simulate!(simulator, T, max_dt)

    expected_cells = [297 / 320, 5 / 6, 629 / 960]

    @test grid.cells ≈ expected_cells
end

@testset "Test riemann problem 2 steps" begin
    bc = NeumannBC()
    u0(x) = x < 0 ? 1 : 0
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))
    max_dt = grid.dx
    T = 2*max_dt
    eq = BurgersEQ()
    F = LaxFriedrichsFlux()
    system = ConservedSystem(eq, grid, F)
    timestepper = ForwardEuler()
    simulator = Simulator(system, grid, timestepper, 0.)

    grid.cells[:] = [1., 0.5, 0.]

    simulate!(simulator, T, max_dt)

    expected_cells = [945 / 1024 , 105 / 128, 663 / 1024]

    @test grid.cells ≈ expected_cells
end