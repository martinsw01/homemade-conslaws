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
    reconstruction = NoReconstruction(grid)
    timestepper = RK2(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)


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
    reconstruction = NoReconstruction(grid)
    timestepper = RK2(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    simulate!(simulator, T, dt)


    grid2 = UniformGrid1D(N, bc, u0, (x_L, x_R))
    grid2.cells[:] = [1., 1., 0., 0.]
    reconstruction2 = NoReconstruction(grid2)
    timestepper2 = ForwardEuler(grid2)
    system2 = ConservedSystem(eq, reconstruction2, F, timestepper2)
    simulator2 = Simulator(system2, grid2, 0.)

    simulate!(simulator2, 2*T, dt)

    @test grid.cells ≈ 0.5 .* ([1., 1., 0., 0.] .+ grid2.cells)

    @test grid.cells ≈ [127/128, 127/128, 33/128, 33/128]
end