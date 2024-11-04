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

    @test grid.cells ≈ [[-2/3, 2/3], [0., 0.], [2/3, -2/3]]

    expected_cells = [[-2/3, 2/3], [0., 0.], [2/3, -2/3]]

    for_each_cell(grid) do cells, cell_idx
        @test cells[cell_idx] ≈ expected_cells[cell_idx] atol=1e-15
    end

    cell_sum = homemade_conslaws.reduce_cells(grid, 0.) do acc, cells, cell_idx
        acc + sum(cells[cell_idx])
    end

    @test cell_sum ≈ 0.

    expected_neighbour_values = [
        [[-2/3, 2/3], [-2/3, 2/3], [0., 0.]],
        [[-2/3, 2/3], [0., 0.], [2/3, -2/3]],
        [[0., 0.], [2/3, -2/3], [2/3, -2/3]]
    ]

    for_each_cell(grid) do cells, cell_idx
        neighbour_cells = cells[get_neighbour_indices(grid, cell_idx, 1)]
        @test expected_neighbour_values[cell_idx] ≈ neighbour_cells
    end
end

@testset "Test flux computation" begin
    eq = ConstantLinearAdvection()
    u0(x) = [x > 0, -x > 0]
    N = 3
    bc = NeumannBC()
    grid = UniformGrid1D(N, bc, u0, (-1, 1))

    F = CentralUpwind()

    reconstruction = NoReconstruction(grid)

    left, right = homemade_conslaws.reconstruct(reconstruction, grid)

    dx = grid.dx
    dt = dx
    T = dt

    expected_flux = [
        [0., -1.],
        [0., -1.],
        [0., 0.],
        [1., 0.],
        [1., 0.]
    ]   

    p = homemade_conslaws.number_of_cells(F)
    for_each_cell(grid) do cells, i
        neighbour_cell_idx = get_neighbour_indices(grid, i, p)
        right_minus = left[neighbour_cell_idx[p+1:end-1]]
        left_minus = right[neighbour_cell_idx[1:p]]
        right_plus = left[neighbour_cell_idx[p+2:end]]
        left_plus = right[neighbour_cell_idx[p+1:end-1]]

        dx = get_dx(grid, i)

        flux_minus = F(eq, left_minus, right_minus, dx, dt)
        flux_plus = F(eq, left_plus, right_plus, dx, dt)

        @test flux_minus ≈ expected_flux[i]
        @test flux_plus ≈ expected_flux[i+1]
    end
    
    reconstruction = NoReconstruction(grid)
    timestepper = ForwardEuler(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    simulate!(simulator, T, dt)

    expected_cells = [[0., 1.], [0., 0.], [0., 0.], [1., 0.]]

    @test grid.cells ≈ expected_cells
end
    

nothing