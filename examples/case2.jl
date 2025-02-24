using SinSWE, homemade_conslaws, StaticArrays, ElasticArrays

function simulate_with_bottom_topography(T, N, H0, B_data)
    grid = SinSWE.CartesianGrid(N; gc=2, boundary=SinSWE.WallBC(), extent=[0.0  5.], )

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
    initial .= [(@SVector [max(initial[i][1], SinSWE.B_cell(B, i+2)), 0.]) for i in eachindex(initial)]
    initial[22:24] .= [(@SVector [h, uh]) for (h, uh) in [(0.9991266339695727, -0.0008150944836122437), (0.960169091665025, -0.23886602368659338), (1.2497810640985383, -1.01483619696961)]]

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


function simulate_with_wall_bc2(T, N, H0)
    grid = SinSWE.CartesianGrid(N-5; gc=2, boundary=SinSWE.WallBC(), extent=[0.0  4.], )
    B_data = zeros(N-5+5)
    backend = make_cpu_backend()
    B = SinSWE.BottomTopography1D(B_data, backend, grid)
    eq = SinSWE.ShallowWaterEquations1D(B)
    rec = SinSWE.LinearReconstruction(1.) # TODO: set to 1.
    flux = SinSWE.CentralUpwind(eq)
    bst = SinSWE.SourceTermBottom()
    conserved_system = SinSWE.ConservedSystem(backend, rec, flux, eq, grid, bst)

    timestepper = SinSWE.ForwardEulerStepper()
    simulator = SinSWE.Simulator(backend, conserved_system, timestepper, grid, cfl=0.01)

    x = SinSWE.cell_centers(grid)
    xf = SinSWE.cell_faces(grid)
    u0(x) = @SVector[H0(x), 0.0]
    initial = u0.(x)

    SinSWE.set_current_state!(simulator, initial)

    SinSWE.simulate_to_time(simulator, T)

    H_left = simulator.system.left_buffer.h[3:end-2]
    H_right = simulator.system.right_buffer.h[3:end-2]
    UH_left = simulator.system.left_buffer.hu[3:end-2]
    UH_right = simulator.system.right_buffer.hu[3:end-2]

    return H_left, H_right, UH_left, UH_right, 1, xf
end

function simulate_with_wall_bc(T, N, H0)
    eq = ShallowWater1D(9.81)

    q0(x) = [H0(x), 0.]
    # bc = WallsBC([[4.,4. +(1/(2N))]])
    bc = WallsBC([[4., 4.5]])
    grid = UniformGrid1DWalls(N, bc, q0, (0., 5.))
    dt = grid.dx

    F = CentralUpwind()

    reconstruction = LinearReconstruction(grid)
    timestepper = RK2(grid)
    system = ConservedSystem(eq, reconstruction, F, timestepper)
    simulator = Simulator(system, grid, 0.)

    Q, t = simulate_and_aggregate!(simulator, T, dt)

    Q_left, Q_right = homemade_conslaws.reconstruct(reconstruction, grid, Q[:,end])

    H_left, UH_left = separate_variables(Q_left)
    H_right, UH_right = separate_variables(Q_right)
    return H_left, H_right, UH_left, UH_right, t, cell_centers(grid)
end

H0(x) = @. 1. + 0.5 * exp(-2. * ((x - 2.))^2)
H0_tall(x) = @. 1. + 1.5 * exp(-2. * ((x - 2.))^2)


function compute_cell_faces(N, x_L, x_R)
    grid = SinSWE.CartesianGrid(N; gc=2, boundary=SinSWE.WallBC(), extent=[x_L  x_R], )
    return SinSWE.cell_faces(grid)
end


function main_wall_bc()
    H_left, H_right, UH_left, UH_right, t, _ = simulate_with_wall_bc(T, N-1, H0_tall)
    xf = compute_cell_faces(N, 0., 5.)

    # homemade_conslaws.Viz.plot_water_walls_bc(H_left, H_right, UH_left, UH_right, xf, 4., "case_2_wave_wall_bc")
    println(extrema(diff(t)))
    # println(minimum(diff(t)))
    homemade_conslaws.Viz.plot_time_step(t[2:end], diff(t), "time-step_case_2_wall_bc")
end
function main_wall_bc2(h0, name)
    H_left, H_right, UH_left, UH_right, t, _ = simulate_with_wall_bc2(T, N-1, h0)
    xf = compute_cell_faces(N, 0., 5.)

    A = zeros(6)
    H_left = [H_left; A]
    H_right = [H_right; A]
    UH_left = [UH_left; A]
    UH_right = [UH_right; A]

    println(maximum(abs.(UH_left)), " ", maximum(abs.(UH_right)))

    homemade_conslaws.Viz.plot_water_walls_bc(H_left, H_right, UH_left, UH_right, xf, 4., name)

    # homemade_conslaws.Viz.plot_time_step(t[2:end], diff(t), "time-step_case_2_wall_bc")
end


function main(wall_height, h0, name)
    B = zeros(N+5)

    wall_start_index = 4N ÷ 5 + 3
    @show wall_start_index

    B[wall_start_index:end] .= wall_height

    H_wall_top_left, H_wall_top_right, UH_wall_top_left, UH_wall_top_right, t_wall_top, x_wall_top = simulate_with_bottom_topography(T, N, h0, B)
    # H_wall_top_left, H_wall_top_right, UH_wall_top_left, UH_wall_top_right, t_wall_top, x_wall_top = simulate_with_bottom_topography(T, N, H0_tall, B)

    H_wall_top_left[:, wall_start_index-2:end] .+= wall_height
    H_wall_top_right[:, wall_start_index-3:end] .+= wall_height

    H = 0.5 .* (H_wall_top_left .+ H_wall_top_right)[end,:]
    UH = 0.5 .* (UH_wall_top_left .+ UH_wall_top_right)[end,:]

    println(maximum(abs.(UH_wall_top_left[end,:])), " ", maximum(abs.(UH_wall_top_right[end,:])))

    # homemade_conslaws.Viz.plot_water(H_wall_top_left[end,:], H_wall_top_right[end,:], UH_wall_top_left[end,:], UH_wall_top_right[end,:], x_wall_top, B[3:end-2], "$name$(Int64(10wall_height))")
    # homemade_conslaws.Viz.plot_water(H_wall_top_left[end,:], H_wall_top_right[end,:], UH_wall_top_left[end,:], UH_wall_top_right[end,:], x_wall_top, B[3:end-2], "case_2_min_dt2_$(Int64(10wall_height))")
    
    @show argmin(diff(t_wall_top))
    @show t_wall_top[argmin(diff(t_wall_top))]
    @show minimum(diff(t_wall_top))
    # @show maximum(diff(t_wall_top))
    homemade_conslaws.Viz.plot_time_step(t_wall_top[2:end-1], diff(t_wall_top)[1:end-1], "time-step_case_2_corrected$(Int64(10wall_height))")
    # homemade_conslaws.Viz.plot_time(t_wall_top, "case_2_submerged_time$(Int64(10wall_height))")

    # # homemade_conslaws.Viz.plot_momentum(UH_wall_top_left[end,:], UH_wall_top_right[end,:], x_wall_top, B[3:end-2])
    # homemade_conslaws.Viz.animate_water_rec2(H_wall_top_left, H_wall_top_right, UH_wall_top_left, UH_wall_top_right, x_wall_top, t_wall_top, B[3:end-2], 7)
    # # homemade_conslaws.Viz.animate_water_rec(H_wall_top_left, H_wall_top_right, UH_wall_top_left, UH_wall_top_right, x_wall_top, t_wall_top, walls, 10)
    # # homemade_conslaws.Viz.animate_solutions((H_const_top, H_wall_top), ["H constant topography" "H wall topography"], x_const_top, t_const_top, 2T)
    # # homemade_conslaws.Viz.animate_water(H_wall_top, UH_wall_top, x_wall_top, t_wall_top, [wall], 2T)


    # H0(x) = @. 1. + 0.5 * exp(-(2. * (x - 2.))^2)

    # N = 100
    # T = 3.
    # B = zeros(N+5)
    # B[83:end] .= wall_height
    # wall_interface_indices = findall(b -> b > 0.5, B)




    # # B[wall_interface_indices] .= wall_height
    # H_wall_top_left, H_wall_top_right, UH_wall_top_left, UH_wall_top_right, t_wall_top, x_wall_top = simulate_with_bottom_topography(T, N, H0, B)

    # H_wall_top_left[:, 81:end] .+= wall_height
    # H_wall_top_right[:, 80:end] .+= wall_height

    # # # homemade_conslaws.Viz.plot_momentum(UH_wall_top_left[end,:], UH_wall_top_right[end,:], x_wall_top, B[3:end-2])
    # # homemade_conslaws.Viz.plot_water(H_wall_top_left[end,:], H_wall_top_right[end,:], UH_wall_top_left[end,:], UH_wall_top_right[end,:], x_wall_top, B[3:end-2], "water_wall_height_$(Int64(10wall_height))")
    # homemade_conslaws.Viz.animate_water_rec2(H_wall_top_left, H_wall_top_right, UH_wall_top_left, UH_wall_top_right, x_wall_top, t_wall_top, B[3:end-2], 3)
    # # # homemade_conslaws.Viz.animate_water_rec(H_wall_top_left, H_wall_top_right, UH_wall_top_left, UH_wall_top_right, x_wall_top, t_wall_top, walls, 10)
    # # # homemade_conslaws.Viz.animate_solutions((H_const_top, H_wall_top), ["H constant topography" "H wall topography"], x_const_top, t_const_top, 2T)
    # # # homemade_conslaws.Viz.animate_water(H_wall_top, UH_wall_top, x_wall_top, t_wall_top, [wall], 2T)
end

N = 30
# T = 0.7313047188009953
T = 1.
# main(1.5, H0_tall, "case_2_corrected")
# main_wall_bc2(H0_tall, "case_2_tall_wave_wall_bc")
main_wall_bc()