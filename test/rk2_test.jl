using homemade_conslaws
using Test

@testset "Test preserves 0" begin
    bc = NeumannBC()
    u0 = [0., 0., 0.]
    N = 3
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))
    dt = grid.dx
    T = dt
    eq = BurgersEQ()
    F = LaxFriedrichsFlux()
    reconstruction = NoReconstruction()
    timestepper = RK2(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)


    simulate!(simulator, T, dt)

    @test cells(grid) ≈ [0., 0., 0.] atol=1e-32
end

@testset "Test one step" begin
    bc = NeumannBC()
    u0 = [1., 1., 0., 0.]
    N = 4
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))
    dt = grid.dx
    T = dt
    eq = BurgersEQ()
    F = LaxFriedrichsFlux()
    reconstruction = NoReconstruction()
    timestepper = RK2(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    simulate!(simulator, T, dt)


    grid2 = UniformGrid1D(N, bc, u0, (x_L, x_R))
    reconstruction2 = NoReconstruction()
    timestepper2 = ForwardEuler(grid2)
    system2 = ConservedSystem(eq, reconstruction2, F, timestepper2)
    simulator2 = Simulator(system2, grid2, 0.)

    simulate!(simulator2, 2*T, dt)

    @test cells(grid) ≈ 0.5 .* ([1., 1., 0., 0.] .+ cells(grid2))

    @test cells(grid) ≈ [127/128, 127/128, 33/128, 33/128]
end