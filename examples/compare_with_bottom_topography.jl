using SinSWE, homemade_conslaws, StaticArrays, ElasticArrays

function simulate_with_bottom_topography(T, N, H0, B_data)
    grid = SinSWE.CartesianGrid(N; gc=2, boundary=SinSWE.WallBC(), extent=[0.  1.], )

    # x0 = 0.5
    # B_data = [x < x0 ? 2. : 0. for x in SinSWE.cell_faces(grid, interior=false)]

    # B_data = zero(SinSWE.cell_faces(grid, interior=false))
    # # B_data[length(B_data) ÷ 2] = 2.
    # B_data[1] = 3.

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
    # H_left = ElasticMatrix(reshape([Q[1] for Q in simulator.system.left_buffer.h], :, 1))
    # H_right = ElasticMatrix(reshape([Q[1] for Q in simulator.system.right_buffer], :, 1))
    # UH_left = ElasticMatrix(reshape([Q[2] for Q in simulator.system.left_buffer], :, 1))
    # UH_right = ElasticMatrix(reshape([Q[2] for Q in simulator.system.right_buffer], :, 1))
    t = ElasticVector([0.])

    function collect_state(t_j, simulator)
        # H_left_next = [Q[1] for Q in simulator.system.left_buffer]
        # H_right_next = [Q[1] for Q in simulator.system.right_buffer]
        # UH_left_next = [Q[2] for Q in simulator.system.left_buffer]
        # UH_right_next = [Q[2] for Q in simulator.system.right_buffer]

        append!(UH_left, simulator.system.left_buffer.hu)
        append!(UH_right, simulator.system.right_buffer.hu)
        append!(H_left, simulator.system.left_buffer.h)
        append!(H_right, simulator.system.right_buffer.h)
        append!(t, t_j)
    end

    SinSWE.simulate_to_time(simulator, T, callback=collect_state)

    return H_left'[:,3:end-2], H_right'[:,3:end-2], UH_left'[:,3:end-2], UH_right'[:,3:end-2], t, xf
end

# function simulate_with_wall_bc(T, N, H0)
#     grid = SinSWE.CartesianGrid(N; gc=2, boundary=SinSWE.WallBC(), extent=[0.0  1.], )

#     backend = make_cpu_backend()
#     eq = SinSWE.ShallowWaterEquations1D(1., 1.)
#     rec = SinSWE.LinearReconstruction(1.)
#     flux = SinSWE.CentralUpwind(eq)
#     timestepper = SinSWE.ForwardEulerStepper()
#     conserved_system = SinSWE.ConservedSystem(backend, rec, flux, eq, grid, [])
#     simulator = SinSWE.Simulator(backend, conserved_system, timestepper, grid, cfl=0.2)

#     x = SinSWE.cell_centers(grid)
#     u0(x) = @SVector[H0(x), 0.0]
#     initial = u0.(x)
#     # initial .= [1.0, 1.0, 0.9925263921686729, 0.9723037772896017, 0.9438200176627859, 0.9113101113849261, 0.8786416537317597, 0.8485718552364452, 0.8226615579783141, 0.802614148189579]
#     SinSWE.set_current_state!(simulator, initial)

#     H = ElasticMatrix(reshape([Q[1] for Q in SinSWE.current_state(simulator)], :, 1))
#     UH = ElasticMatrix(reshape([Q[2] for Q in SinSWE.current_state(simulator)], :, 1))
#     t = ElasticVector([0.])

#     function collect_state(t_j, simulator)
#         UH_next = [Q[2] for Q in SinSWE.current_state(simulator)]
#         H_next = [Q[1] for Q in SinSWE.current_state(simulator)]

#         append!(UH, UH_next)
#         append!(H, H_next)
#         append!(t, t_j)
#     end 

#     SinSWE.simulate_to_time(simulator, T, callback=collect_state)
#     return H', UH', t, x, reconstruct
# end


function compute_steady_state_max_momentum(H0, N)
    B = zeros(N+5)
    T = 5.
    wall_height = 1.5
    wall_interface_indices = ((N÷10*3):(N÷10*4)) .+ 3
    B[wall_interface_indices] .= wall_height
    H_right, H_right, UH_left, UH_right, t, x = simulate_with_bottom_topography(T, N, H0, B)

    # initial_momentum = max(maximum(abs.(UH_left[20,:])), maximum(abs.(UH_right[20,:])))
    # final_momentum = max(maximum(abs.(UH_left[end,:])), maximum(abs.(UH_right[end,:])))

    final_momentum = (sum(UH_left[end,:] .+ UH_right[end,:])) / 2

    # return [initial_momentum, final_momentum]
    return final_momentum
end


function compare_wall_height(wall_height)
    H0(x) = 1.# + 0.1x

    N = 10
    T = 0.1
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

function main3()
    # wall_heights = 1.:0.1:3.
    wall_heights = 1.1:0.1:3.
    dt = compare_wall_height.(wall_heights)
    dt_ref = compare_wall_height(0.)
    @show dt_ref
    dt, wall_heights
    homemade_conslaws.Viz.plot_timesteps(wall_heights, dt, dt_ref)
end


function compare_steady_state_max_momenta()
    H0(x) = 1.
    N = 10:30:100
    UH = compute_steady_state_max_momentum.(H0, N)
    homemade_conslaws.Viz.plot_momenta(UH, N)
end

main2() = compare_steady_state_max_momenta()





function main(wall_height)
    H0(x) = 1.# + 0.1x

    N = 10
    T = 10.
    B = zeros(N+5)
    # H_const_top, UH_const_top, t_const_top, x_const_top, _ = simulate_with_bottom_topography(T, N, H0, B)

    # wall_height = 3.
    # wall_interface_indices = (N+5) ÷ 2 .+ (1:1)
    wall_interface_indices = (N+5) ÷ 2 .+ (-1:1)
    # # wall_interface_indices = 6:8
    # wall_interface_indices = 32:52
    # @show wall_interface_indices
    B[wall_interface_indices] .= wall_height
    H_wall_top_left, H_wall_top_right, UH_wall_top_left, UH_wall_top_right, t_wall_top, x_wall_top = simulate_with_bottom_topography(T, N, H0, B)
    # @show x_wall_top
    # @show x_wall_top[wall_interface_indices .- 2]

    # H_wall_top_left[:, N ÷ 2 .+ (-1:1)] .+= wall_height
    # H_wall_top_right[:, N ÷ 2 .+ (-2:0)] .+= wall_height

    # @show wall_interface_indices[2:end] .- 3

    H_wall_top_left[:, wall_interface_indices .- 2] .+= wall_height
    H_wall_top_right[:, wall_interface_indices .- 3] .+= wall_height
    # walls = [x_wall_top[N ÷ 2 .+ [-1, 1]]]
    

    # H_wall_bc, UH_wall_bc, t_wall_bc, x_wall_bc = simulate_with_wall_bc(T, N, H0)

    # println(diff(t_wall_top[1:10])[3])
    # println(diff(t_wall_top[end-10:end])[3])

    println(maximum(abs.(UH_wall_top_left[end,:])))
    println(maximum(abs.(UH_wall_top_right[end,:])))

    # homemade_conslaws.Viz.plot_momentum(UH_wall_top_left[end,:], UH_wall_top_right[end,:], x_wall_top, B[3:end-2])
    homemade_conslaws.Viz.plot_water(H_wall_top_left[end,:], H_wall_top_right[end,:], UH_wall_top_left[end,:], UH_wall_top_right[end,:], x_wall_top, B[3:end-2], "water_wall_height_$(Int64(10wall_height))")
    # homemade_conslaws.Viz.animate_water_rec2(H_wall_top_left, H_wall_top_right, UH_wall_top_left, UH_wall_top_right, x_wall_top, t_wall_top, B[3:end-2], 3)
    # homemade_conslaws.Viz.animate_water_rec(H_wall_top_left, H_wall_top_right, UH_wall_top_left, UH_wall_top_right, x_wall_top, t_wall_top, walls, 10)
    # homemade_conslaws.Viz.animate_solutions((H_const_top, H_wall_top), ["H constant topography" "H wall topography"], x_const_top, t_const_top, 2T)
    # homemade_conslaws.Viz.animate_water(H_wall_top, UH_wall_top, x_wall_top, t_wall_top, [wall], 2T)
end

for height in [0.9, 1.1, 1.3, 1.7, 2.0]
    main(height)
end

# main(10.)