using SinSWE, homemade_conslaws, StaticArrays, ElasticArrays, LaTeXStrings

function compute_steady_state_solution(T, N, H0, B_data)
    grid = SinSWE.CartesianGrid(N; gc=2, boundary=SinSWE.WallBC(), extent=[0.0  1.], )

    backend = make_cpu_backend()
    B = SinSWE.BottomTopography1D(B_data, backend, grid)
    eq = SinSWE.ShallowWaterEquations1D(B)
    rec = SinSWE.LinearReconstruction(1.)
    flux = SinSWE.CentralUpwind(eq)
    bst = SinSWE.SourceTermBottom()
    conserved_system = SinSWE.ConservedSystem(backend, rec, flux, eq, grid, bst)

    timestepper = SinSWE.ForwardEulerStepper()
    simulator = SinSWE.Simulator(backend, conserved_system, timestepper, grid, cfl=0.2)

    x = SinSWE.cell_centers(grid)
    xf = SinSWE.cell_faces(grid)
    u0(x) = @SVector[H0, 0.0]
    initial = u0.(x)

    initial .= [(@SVector [max(initial[i][1], SinSWE.B_cell(B, i+2)), 0.]) for i in eachindex(initial)]

    SinSWE.set_current_state!(simulator, initial)

    SinSWE.simulate_to_time(simulator, T)

    H = SinSWE.current_interior_state(simulator).h
    HU = SinSWE.current_interior_state(simulator).hu
    return H, HU, xf
end


function compute_max_error(H, HU, lakeatrest_height, wall_indices)
    H_difference = copy(H) .- lakeatrest_height
    HU_difference = copy(HU)

    H_difference[wall_indices] .= 0.
    HU_difference[wall_indices] .= 0.

    H_max_error = maximum(abs.(H_difference))
    HU_max_error = maximum(abs.(HU_difference))
    return H_max_error, HU_max_error
end


function compute_square_error(H, HU, lakeatrest_height, wall_indices)
    H_difference = copy(H) .- lakeatrest_height
    HU_difference = copy(HU)

    H_difference[wall_indices] .= 0.
    HU_difference[wall_indices] .= 0.

    square_error_H = sqrt(sum(abs2, H_difference)/length(H))
    square_error_HU = sqrt(sum(abs2, HU_difference)/length(HU))
    return square_error_H, square_error_HU
end


function compute_abs_error(H, HU, lakeatrest_height, wall_indices)
    H_difference = copy(H) .- lakeatrest_height
    HU_difference = copy(HU)

    H_difference[wall_indices] .= 0.
    HU_difference[wall_indices] .= 0.

    abs_error_H = sum(abs, H_difference)/length(H)
    abs_error_HU = sum(abs, HU_difference)/length(HU)
    return abs_error_H, abs_error_HU
end


function compute_errors(grid_sizes::AbstractArray, compute_error, wall_height)
    errors_H = Vector{Float64}(undef, length(grid_sizes))
    errors_HU = Vector{Float64}(undef, length(grid_sizes))

    for (i, N) in enumerate(grid_sizes)
        @show i, N
        T = 10.
        B = zeros(N+5)

        wall_interface_indices = (N+5) ÷ 2 .+ (1:1)

        H0 = 1.
        B[wall_interface_indices] .= wall_height
        H, HU, xf = compute_steady_state_solution(T, N, H0, B)

        H_error, HU_error = compute_error(H, HU, H0, wall_interface_indices .- 2)
        errors_H[i] = H_error
        errors_HU[i] = HU_error
    end

    return errors_H, errors_HU
end

function compute_errors(N, compute_error, wall_heights::AbstractArray)
    errors_H = Vector{Float64}(undef, length(wall_heights))
    errors_HU = Vector{Float64}(undef, length(wall_heights))

    for (i, wall_height) in enumerate(wall_heights)
        T = 10.
        B = zeros(N+5)

        wall_interface_indices = (N+5) ÷ 2 .+ (1:1)

        H0 = 1.
        B[wall_interface_indices] .= wall_height
        H, HU, xf = compute_steady_state_solution(T, N, H0, B)

        H_error, HU_error = compute_error(H, HU, H0, wall_interface_indices .- 2)
        errors_H[i] = H_error
        errors_HU[i] = HU_error
    end

    return errors_H, errors_HU
end


function main()
    # N = 10
    # N = [10, 15, 20, 25]
    N = floor.(Int64, exp.(2.5:0.5:6.))
    # N = floor.(Int64, exp.(2.:0.5:6.))
    # wall_heights = 1.1:0.5:3#floor.(Int64, exp.(0.2:0.3:2.3))
    wall_heights = exp.(0.1:0.15:2.3)
    wall_height = 3.
    # errors_H, errors_HU = compute_errors(N, compute_max_error, wall_heights)
    # errors_H, errors_HU = compute_errors(N, compute_max_error, wall_heights)
    errors_H, errors_HU = compute_errors(N, compute_max_error, wall_height)
    # errors_H, errors_HU = compute_errors(N, compute_max_error)
    # homemade_conslaws.Viz.plot_error_wall(wall_heights, errors_H, errors_HU, L"L^\infty"*"-error", "max_error_walls")
    # homemade_conslaws.Viz.plot_error_wall(wall_heights, errors_H, errors_HU, L"L^\infty"*"-error", "max_error_walls")
    homemade_conslaws.Viz.plot_error(N, errors_H, errors_HU, L"L^\infty"*"-error", "max_error")
    # @show errors
    # @show errors_H
    # @show errors_HU
    # nothing
end

main()