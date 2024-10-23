using homemade_conslaws
using Test

@testset "Test preserves 0" begin
    bc = NeumannBC()
    u0(x) = 0
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))
    dt = grid.dx
    T = dt
    eq = BurgersEQ()
    F = LaxFriedrichsFlux()
    system = ConservedSystem(eq, grid, F)
    timestepper = RK2()
    simulator = Simulator(system, grid, timestepper, 0.)

    simulate!(simulator, T, dt)

    @test grid.cells ≈ [0., 0., 0.] atol=1e-32
end

@testset "Test one step" begin
    bc = NeumannBC()
    u0(x) = x < 0
    N = 3
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))
    grid.cells[:] = [1., 1., 0., 0.]
    dt = grid.dx
    T = dt
    eq = BurgersEQ()
    F = LaxFriedrichsFlux()
    system = ConservedSystem(eq, grid, F)
    timestepper = RK2()
    simulator = Simulator(system, grid, timestepper, 0.)

    simulate!(simulator, T, dt)

    @test grid.cells ≈ [255/256, 263/256, 65/256, 57/256]
end