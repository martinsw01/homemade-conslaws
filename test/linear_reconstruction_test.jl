using Test
using homemade_conslaws

@testset "Test shock reconstruction" begin
    bc = NeumannBC()
    u0(x) = x < 0 ? 1 : 0
    N = 3
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    @test inner_cells(grid) ≈ [1., 1., 0., 0.]

    reconstruction = LinearReconstruction(grid)

    left, right = homemade_conslaws.reconstruct(reconstruction, grid)

    @test inner_cells(grid, left) ≈ [1., 1., 0., 0.] ≈ inner_cells(grid, right)
end

@testset "Test smooth reconstruction" begin
    bc = NeumannBC()
    u0 = identity
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    reconstruction = LinearReconstruction(grid)

    left, right = homemade_conslaws.reconstruct(reconstruction, grid, [1., 1., 0.5, 0., 0.])

    @test inner_cells(grid, left) ≈ [1., 0.75, 0.]
    @test inner_cells(grid, right) ≈ [1., 0.25, 0.]

    left, right = homemade_conslaws.reconstruct(reconstruction, grid)

    @test inner_cells(grid, left) ≈ [-1., -1/3, 1/3]
    @test inner_cells(grid, right) ≈ [-1/3, 1/3, 1.]

end

@testset "Test system reconstruction" begin
    bc = NeumannBC()
    u0(x) = [x, -x]
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R), 2)

    reconstruction = LinearReconstruction(grid)

    left, right = homemade_conslaws.reconstruct(reconstruction, grid)

    @test inner_cells(grid, left) ≈ [[-1., 1.], [-1/3, 1/3], [1/3, -1/3]]
    @test inner_cells(grid, right) ≈ [[-1/3, 1/3], [1/3, -1/3], [1., -1.]]
end