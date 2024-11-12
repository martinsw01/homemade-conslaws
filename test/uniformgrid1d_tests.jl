using homemade_conslaws
using Test

@testset "Test cell vales at shock" begin
    bc = PeriodicBC()
    u0(x) = x < 0
    N = 3
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    @test inner_cells(grid) ≈ [1, 1, 0, 0]
end

@testset "Test grid periodic BC" begin
    bc = PeriodicBC()
    u0(x) = x < 0 ? 1 : 0
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    expected_cells = [1, 0.5, 0]
    expected_dx = (x_R - x_L) / (N + 1)

    @test inner_cells(grid) ≈ expected_cells
    @test grid.dx ≈ expected_dx
    for_each_inner_cell(grid) do cells, _, cell_idx, _
        @test cells[cell_idx] ≈ expected_cells[cell_idx-1]
        @test get_dx(grid, cell_idx) ≈ expected_dx
    end

    homemade_conslaws.update_bc!(grid, BurgersEQ())

    @test grid.cells ≈ [0, 1, 0.5, 0, 1]
end

@testset "Test update_bc two ghost cells, PeriodicBC" begin
    bc = PeriodicBC()
    u0(x) = x
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R), ghost_cells=2)

    homemade_conslaws.update_bc!(grid, BurgersEQ())

    @test grid.cells ≈ [0, 2/3, -2/3, 0, 2/3, -2/3, 0]
end

@testset "Test update_bc two ghost cells, NeumannBC" begin
    bc = NeumannBC()
    u0(x) = x
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R); ghost_cells=2)

    homemade_conslaws.update_bc!(grid, BurgersEQ())

    @test grid.cells ≈ [-2/3, -2/3, -2/3, 0, 2/3, 2/3, 2/3]
end

@testset "Test grid neumann BC" begin
    bc = NeumannBC()
    u0(x) = x
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    homemade_conslaws.update_bc!(grid, BurgersEQ())

    @test grid.cells ≈ [-2/3, -2/3, 0., 2/3, 2/3]
end

@testset "Test grid Wall BC" begin
    bc = WallBC()
    u0(x) = [x, -x]
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R); ghost_cells=2)

    homemade_conslaws.update_bc!(grid, ShallowWater1D(1.))

    @test grid.cells ≈ [[0, 0], [-2/3, -2/3], [-2/3, 2/3], [0, 0], [2/3, -2/3], [2/3, 2/3], [0, 0]]
end

