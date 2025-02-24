using SinSWE, homemade_conslaws, StaticArrays, ElasticArrays

function simulate_with_bottom_topography(T, N, H0, B_data)
    grid = SinSWE.CartesianGrid(N; gc=2, boundary=SinSWE.WallBC(), extent=[0.0  1.], )

    backend = make_cpu_backend()
    B = SinSWE.BottomTopography1D(B_data, backend, grid)
    eq = SinSWE.ShallowWaterEquations1D(B)
    rec = SinSWE.LinearReconstruction(1.) # TODO: set to 1.
    flux = SinSWE.CentralUpwind(eq)
    bst = SinSWE.SourceTermBottom()
    conserved_system = SinSWE.ConservedSystem(backend, rec, flux, eq, grid, bst)

    timestepper = SinSWE.ForwardEulerStepper()
    simulator = SinSWE.Simulator(backend, conserved_system, timestepper, grid, cfl=0.2)

    x = SinSWE.cell_centers(grid)
    xf = SinSWE.cell_faces(grid)
    u0(x) = @SVector[H0(x), 0.0]
    initial = u0.(x)
    # initial[6] = @SVector[1.0, -1.5]
    # initial[7] = @SVector[1.0, 1.5]
    # initial .= [(@SVector [h, 0.0]) for h in [1.0, 1.0, 0.9925263921686729, 0.9723037772896017, 0.943940820816304, 0.9121840966367951, 0.8815082446828364, 0.8567190670896799, 0.8400837887632049, 0.8303343963654161]]
    initial .= [(@SVector [max(initial[i][1], SinSWE.B_cell(B, i+2)), 0.]) for i in eachindex(initial)]

    SinSWE.set_current_state!(simulator, initial)

    SinSWE.reconstruct!(backend, rec, simulator.system.left_buffer, simulator.system.right_buffer, SinSWE.current_state(simulator), grid, eq, SinSWE.XDIR)

    H_left = ElasticMatrix(reshape(simulator.system.left_buffer.h, :, 1))
    H_right = ElasticMatrix(reshape(simulator.system.right_buffer.h, :, 1))
    UH_left = ElasticMatrix(reshape(simulator.system.left_buffer.hu, :, 1))
    UH_right = ElasticMatrix(reshape(simulator.system.right_buffer.hu, :, 1))
    t = ElasticVector([0.])

    function collect_state(t_j, simulator)
        append!(UH_left, simulator.system.left_buffer.hu)
        append!(UH_right, simulator.system.right_buffer.hu)
        append!(H_left, simulator.system.left_buffer.h)
        append!(H_right, simulator.system.right_buffer.h)
        append!(t, t_j)
    end

    SinSWE.simulate_to_time(simulator, T, callback=collect_state)

    return H_left'[:,3:end-2], H_right'[:,3:end-2], UH_left'[:,3:end-2], UH_right'[:,3:end-2], t, xf
end


function compare_wall_height(wall_height)
    H0(x) = 1.# + 0.1x

    N = 10
    T = 10
    B = zeros(N+5)

    wall_interface_indices = (N+5) ÷ 2 .+ (-1:1)

    B[wall_interface_indices] .= wall_height
    H_wall_top_left, H_wall_top_right, UH_wall_top_left, UH_wall_top_right, t_wall_top, x_wall_top = simulate_with_bottom_topography(T, N, H0, B)

    average_left_momentum = 0.5 * (UH_wall_top_left[end,3] + UH_wall_top_right[end,3])
    # average_left_hight = 0.5 * (H_wall_top_left[end,3] + wall_height) #H_wall_top_right[end,3])
    # @show average_left_momentum
    # @show average_left_hight
    # left_momentum_difference = UH_wall_top_left[end,3] - UH_wall_top_right[end,3]
    # left_height_difference = H_wall_top_left[end,3] - H_wall_top_right[end,3]
    # @show left_momentum_difference
    # @show left_height_difference
    diff(t_wall_top[end-10:end])[3]
end

function main()
    # wall_heights = 1.:0.1:3.
    # wall_heights = 0.9:0.03:10.
    wall_heights = 0.9:0.3:10.
    dt = compare_wall_height.(wall_heights)
    dt_ref = compare_wall_height(0.)
    @show dt_ref
    dt, wall_heights
    homemade_conslaws.Viz.plot_timesteps(wall_heights, dt, dt_ref)
end

main()