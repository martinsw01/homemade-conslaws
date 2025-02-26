using Test
using homemade_conslaws

@testset "Test shock reconstruction" begin
    bc = NeumannBC()
    u0 = [1., 1., 0., 0.]
    N = 3
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    @test cells(grid) ≈ [1., 1., 0., 0.]

    reconstruction = LinearReconstruction(grid)

    left, right = homemade_conslaws.reconstruct(reconstruction, grid)

    @test left ≈ [1., 1., 0., 0.] ≈ right
end

@testset "Test smooth reconstruction" begin
    bc = NeumannBC()
    u0 = [-2/3, 0, 2/3]
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    reconstruction = LinearReconstruction(grid)

    left, right = homemade_conslaws.reconstruct(reconstruction, grid, [1., 0.5, 0.])

    @test left ≈ [1., 0.75, 0.]
    @test right ≈ [1., 0.25, 0.]

    left, right = homemade_conslaws.reconstruct(reconstruction, grid)

    @test left ≈ [-2/3, -1/3, 2/3]
    @test right ≈ [-2/3, 1/3, 2/3]

end

@testset "Test system reconstruction" begin
    bc = NeumannBC()
    u0 = [[x, -x] for x in [-2/3, 0, 2/3]]
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    reconstruction = LinearReconstruction(grid)

    left, right = homemade_conslaws.reconstruct(reconstruction, grid)

    @test left ≈ [[-2/3, 2/3], [-1/3, 1/3], [2/3, -2/3]]
    @test right ≈ [[-2/3, 2/3], [1/3, -1/3], [2/3, -2/3]]
end