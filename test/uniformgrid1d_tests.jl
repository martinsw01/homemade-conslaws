using homemade_conslaws
using Test

@testset "Test grid periodic BC" begin
    bc = PeriodicBC()
    u0 = [-0.5, 0.5]
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    expected_stencils = [(0.5, -0.5, 0.5),
                         (-0.5, 0.5, -0.5)]

    for_each_boundary_cell(grid) do stencil, idx
        @test all(stencil .≈ expected_stencils[idx])
    end
end


@testset "Test grid neumann BC" begin
    bc = NeumannBC()
    u0 = [-0.5, 0.5]
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    expected_stencils = [(-0.5, -0.5, 0.5),
                         (-0.5, 0.5, 0.5)]

    for_each_boundary_cell(grid) do stencil, idx
        @test all(stencil .≈ expected_stencils[idx])
    end
end


@testset "Test grid Wall BC" begin
    bc = WallBC()
    u0 = [[-0.5, 0.5], [0.5, -0.5]]
    N = 2
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    expected_stencils = [([-0.5, -0.5], [-0.5, 0.5], [0.5, -0.5]),
                         ([-0.5, 0.5], [0.5, -0.5], [0.5, 0.5])]

    for_each_boundary_cell(grid) do stencil, idx
        @test all(stencil .≈ expected_stencils[idx])
    end
end

