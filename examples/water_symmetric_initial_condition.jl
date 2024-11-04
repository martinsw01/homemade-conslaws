using homemade_conslaws
using StaticArrays
using ElasticArrays

function main()
    eq = ShallowWater1D(1.)   # Gravity is 1
    h0(x) = 1. + 0.3exp(-10x^2)
    u0(x) = @SVector [h0(x), 0.] # Constaint height of 1, no momentum
    N = 100
    bc = PeriodicBC()            # Periodic boundary condition actes as a wall for symmetric initial value
    grid = UniformGrid1D(N, bc, u0, (-1, 1))
    dt = grid.dx
    T = 2.

    F = CentralUpwind()

    reconstruction = NoReconstruction(grid)
    timestepper = ForwardEuler(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    U = ElasticMatrix(reshape([Q[1] for Q in inner_cells(grid)], :, 1))
    t = ElasticVector([0.])

    function add_solution(simulator)
        append!(U, [Q[1] for Q in inner_cells(grid)])
        append!(t, simulator.t[])
    end

    simulate!(simulator, T, dt, [add_solution])

    homemade_conslaws.Viz.animate_solution(
        (U',),
        "h",
        -1:grid.dx:1-grid.dx,
        t,
        3
    )
end

main()