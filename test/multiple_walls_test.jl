using homemade_conslaws
using Test

@testset "Test no interior walls" begin
    bc = WallsBC()
    u0(x) = [x, -x]
    N = 3
    x_L, x_R = -1, 1

    grid = UniformGrid1DWalls(N, bc, u0, (x_L, x_R))

    expected_cell_centers = [-0.75, -0.25, 0.25, 0.75]
    expected_dx = 1/2

    @test grid.dx ≈ expected_dx
    @test cell_centers(grid) ≈ expected_cell_centers

    expected_left_of_walls_indices = [4]
    expected_right_of_walls_indices = [1]

    @test grid.left_of_walls_indices == expected_left_of_walls_indices
    @test grid.right_of_walls_indices == expected_right_of_walls_indices


    expected_inner_cells_indices = 1:4
    @test grid.cells_not_containing_walls_indices == expected_inner_cells_indices
end


@testset "Test two interior walls respectively starting and stopping at interface" begin
    bc = WallsBC([[-0.5, -0.4], [0.4, 0.5]])
    u0(x) = [x, -x]
    N = 7
    x_L, x_R = -1, 1

    grid = UniformGrid1DWalls(N, bc, u0, (x_L, x_R))

    expected_left_of_walls_indices = [2, 5, 8]
    expected_right_of_walls_indices = [1, 4, 7]

    @test grid.left_of_walls_indices == expected_left_of_walls_indices
    @test grid.right_of_walls_indices == expected_right_of_walls_indices

    @test grid.cells_not_containing_walls_indices == [1, 2, 4, 5, 7, 8]
end

@testset "Test one interior wall entirely contained in cell" begin
    bc = WallsBC([[0.1, 0.2]])
    u0(x) = [x, -x]
    N = 3
    x_L, x_R = -1, 1

    grid = UniformGrid1DWalls(N, bc, u0, (x_L, x_R))

    expected_left_of_walls_indices = [2, 4]
    expected_right_of_walls_indices = [1, 4]

    @test grid.left_of_walls_indices == expected_left_of_walls_indices
    @test grid.right_of_walls_indices == expected_right_of_walls_indices

    expected_inner_cells_indices = [1, 2, 4]
    @test grid.cells_not_containing_walls_indices == expected_inner_cells_indices
end

@testset "Test wall covering multiple cells" begin
    bc = WallsBC([[-0.4, 0.4]])
    u0(x) = [x, -x]
    N = 7
    x_L, x_R = -1, 1
    grid = UniformGrid1DWalls(N, bc, u0, (x_L, x_R))

    expected_left_of_walls_indices = [2, 8]
    expected_right_of_walls_indices = [1, 7]

    @test grid.left_of_walls_indices == expected_left_of_walls_indices
    @test grid.right_of_walls_indices == expected_right_of_walls_indices

    expected_inner_cells_indices = [1, 2, 7, 8]
    @test grid.cells_not_containing_walls_indices == expected_inner_cells_indices
end

@testset "Test infinitesimal wall" begin
    bc = WallsBC([[0., 0.]])
    u0(x) = [x, -x]
    N = 3
    x_L, x_R = -1, 1
    grid = UniformGrid1DWalls(N, bc, u0, (x_L, x_R))

    @test grid.left_of_walls_indices == [2, 4]
    @test grid.right_of_walls_indices == [1, 3]
    @test grid.cells_not_containing_walls_indices == 1:4
end


function simulate_simple_wall_bc(N, T, F, q0, x_L, x_R, eq)
    bc = WallBC()
    grid = UniformGrid1D(N, bc, q0, (x_L, x_R))
    cells(grid)[:] = q0.(cell_centers(grid))
    dt = grid.dx

    reconstruction = NoReconstruction(grid)
    timestepper = ForwardEuler(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    Q, t = simulate_and_aggregate!(simulator, T, dt)
    H, UH = separate_variables(Q)
    return H, UH, cell_centers(grid), t
end


function simulate_new_walls_bc(N, T, F, q0_left, x_L, x_mid, x_R, eq)
    q0(x) = x <= x_mid ? q0_left(x) : [1., 0.]

    bc = WallsBC([[x_mid, x_mid]])
    grid = UniformGrid1DWalls(N, bc, q0, (x_L, x_R))
    dt = grid.dx

    reconstruction = NoReconstruction(grid)
    timestepper = ForwardEuler(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    Q, t = simulate_and_aggregate!(simulator, T, dt)
    H, UH = separate_variables(Q)
    return H, UH, cell_centers(grid), t
end


@testset "Compare with simple wall boundary condition" begin
    """
    Place a wall at x=0 and compare with the the solution for WallBC on [0,1].
    """
    x_L, x_mid, x_R = -1., 0., 1.
    N = 4
    T = 50.
    q0_left(x) = [1. + 0.1x, 0.5x]

    eq = ShallowWater1D(1.)
    F = CentralUpwind()

    H_simple_left, UH_simple_left, x_simple_left, t_simple_left = simulate_simple_wall_bc(N, T, F, q0_left, x_L, x_mid, eq)

    H, UH, x, t = simulate_new_walls_bc(2N+1, T, F, q0_left, x_L, x_mid, x_R, eq)

    @test H[:,1:N+1] ≈ H_simple_left
    @test UH[:,1:N+1] ≈ UH_simple_left
    @test x[1:N+1] ≈ x_simple_left
    @test t ≈ t_simple_left
end

nothing