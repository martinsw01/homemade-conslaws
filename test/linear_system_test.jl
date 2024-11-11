using Test
using homemade_conslaws

struct ConstantLinearAdvection <: Equation
    a::Matrix{Float64}
    eigvals::Vector{Float64}
    ConstantLinearAdvection(a = [1. 0; 0. -1.], eigvals = [1., -1.]) = new(a, eigvals)
end

(eq::ConstantLinearAdvection)(u) = eq.a * u

homemade_conslaws.compute_eigenvalues(eq::ConstantLinearAdvection, U) = eq.eigvals

homemade_conslaws.compute_max_abs_eigenvalue(eq::ConstantLinearAdvection, U) = maximum(abs.(eq.eigvals))

homemade_conslaws.conserved_variables(::ConstantLinearAdvection) = (:U,)


@testset "Test cell values" begin
    eq = ConstantLinearAdvection()
    u0(x) = [x, -x]
    N = 2
    bc = NeumannBC()
    grid = UniformGrid1D(N, bc, u0, (-1, 1))

    @test inner_cells(grid) ≈ [[-2/3, 2/3], [0., 0.], [2/3, -2/3]]

    expected_cells = [[-2/3, 2/3], [0., 0.], [2/3, -2/3]]

    for_each_inner_cell(grid) do cells, _, cell_idx, _
        @test cells[cell_idx] ≈ expected_cells[cell_idx-1] atol=1e-15
    end

    cell_sum = homemade_conslaws.reduce_cells(grid, 0.) do acc, cells, cell_idx
        acc + sum(cells[cell_idx])
    end

    @test cell_sum ≈ 0.
end

@testset "Test flux computation" begin
    eq = ConstantLinearAdvection()
    u0(x) = [x > 0, -x > 0]
    N = 3
    bc = NeumannBC()
    grid = UniformGrid1D(N, bc, u0, (-1, 1))
    dx = grid.dx
    dt = dx
    T = dt

    F = CentralUpwind()    
    reconstruction = NoReconstruction(grid)
    timestepper = ForwardEuler(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    simulate!(simulator, T, dt)

    expected_cells = [[0., 1.], [0., 0.], [0., 0.], [1., 0.]]

    @test inner_cells(grid) ≈ expected_cells
end