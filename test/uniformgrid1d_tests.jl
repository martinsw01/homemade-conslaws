using homemade_conslaws
using Test

@testset "Test cell vales at shock" begin
    bc = PeriodicBC()
    u0(x) = x < 0
    N = 3
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    @test grid.cells ≈ [1, 1, 0, 0]
end

@testset "Test grid periodic BC" begin
    bc = PeriodicBC()
    u0(x) = x < 0 ? 1 : 0
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    expected_cells = [1, 0.5, 0]
    expected_dx = (x_R - x_L) / (N + 1)

    @test grid.cells ≈ expected_cells
    @test grid.dx ≈ expected_dx
    for_each_cell(grid) do cells, cell_idx
        @test cells[cell_idx] ≈ expected_cells[cell_idx]
        @test get_dx(grid, cell_idx) ≈ expected_dx
    end
    for_each_cell(grid) do cells, center_idx1
        left_idx, center_idx2, right_idx = get_neighbour_indices(grid, center_idx1, 1)

        @test mod1(left_idx + 1, N+1) == center_idx1 == center_idx2 == mod1(right_idx - 1, N+1)
        @test cells[center_idx1] ≈ expected_cells[center_idx1]
        @test get_dx(grid, center_idx1) ≈ expected_dx
    end

    for_each_cell(grid) do cells, center_idx1
        left_idx, _, center_idx2, _, right_idx = get_neighbour_indices(grid, center_idx1, 2)

        @test mod1(left_idx + 2, N+1) == center_idx1 == center_idx2 == mod1(right_idx - 2, N+1)
    end
end

@testset "Test grid neumann BC" begin
    bc = NeumannBC()
    u0(x) = x
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    expected_neighbour_indices = [
        [1, 1, 1, 2, 3],
        [1, 1, 2, 3, 3],
        [1, 2, 3, 3, 3]
    ]
    for_each_cell(grid) do cells, cell_idx
        @test get_neighbour_indices(grid, cell_idx, 2) == expected_neighbour_indices[cell_idx]
    end
end

