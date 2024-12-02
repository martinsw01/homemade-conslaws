
# using Revise
using homemade_conslaws

function main()
    N = 50

    eq = ShallowWater1D(1.)
    h0(x) = 1. + 0.3exp(-10x^2)
    q0(x) = [h0(x), 0.]
    bc = WallsBC([[-0.5,-0.1]])
    grid = UniformGrid1DWalls(N, bc, q0, (-1, 1))
    dt = grid.dx
    T = 2.

    F = CentralUpwind()

    reconstruction = LinearReconstruction(grid)
    timestepper = RK2(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    Q, t = simulate_and_aggregate!(simulator, T, dt)
    H, UH = separate_variables(Q)

    homemade_conslaws.Viz.animate_solution(
        (H,),
        "h",
        cell_centers(grid),
        t
    )
end

main()