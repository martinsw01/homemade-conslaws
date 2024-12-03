# using Revise
using homemade_conslaws
using StaticArrays
using ElasticArrays

function main()
    eq = ShallowWater1D(1.)   # Gravity is 1
    h0(x) = 1. + 0.3exp(-10x^2)
    q0(x) = @SVector [h0(x), 0.]
    N = 100
    bc = WallBC()
    grid = UniformGrid1D(N, bc, q0, (-1, 1))
    dt = grid.dx
    T = 2.

    F = CentralUpwind()

    reconstruction = LinearReconstruction(grid)
    timestepper = RK2(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    Q, t = simulate_and_aggregate!(simulator, T, dt)
    H, UH = separate_variables(Q)

    homemade_conslaws.Viz.animate_water(
        H, UH,
        grid,
        t
    )
end

main()