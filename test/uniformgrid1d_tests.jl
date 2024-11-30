using homemade_conslaws
using Test

@testset "Test cell vales at shock" begin
    bc = PeriodicBC()
    u0(x) = x < 0
    N = 3
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    @test cells(grid) ≈ [1, 1, 0, 0]
end

@testset "Test grid periodic BC" begin
    bc = PeriodicBC()
    u0(x) = x
    N = 1
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
    u0(x) = x
    N = 1
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
    u0(x) = [x, -x]
    N = 1
    x_L, x_R = -1, 1
    grid = UniformGrid1D(N, bc, u0, (x_L, x_R))

    expected_stencils = [([-0.5, -0.5], [-0.5, 0.5], [0.5, -0.5]),
                         ([-0.5, 0.5], [0.5, -0.5], [0.5, 0.5])]

    for_each_boundary_cell(grid) do stencil, idx
        @test all(stencil .≈ expected_stencils[idx])
    end
end

