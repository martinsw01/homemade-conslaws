using Test
using homemade_conslaws

@testset "Test shock reconstruction" begin
    bc = NeumannBC()
    u0(x) = x < 0 ? 1 : 0
    N = 3
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    @test grid.cells ≈ [1., 1., 0., 0.]

    reconstruction = LinearReconstruction(grid)

    left, right = homemade_conslaws.reconstruct(reconstruction, grid)

    @test left ≈ [1., 1., 0., 0.] ≈ right
end

@testset "Test smooth reconstruction" begin
    bc = NeumannBC()
    u0 = identity
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    reconstruction = LinearReconstruction(grid)

    left, right = homemade_conslaws.reconstruct(reconstruction, grid, [1., 0.5, 0.])

    @test left ≈ [1., 0.75, 0.]
    @test right ≈ [1., 0.25, 0.]
end