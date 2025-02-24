include("shallow_water_recipe.jl")
include("wall_recipe.jl")
include("bathometry_recipe.jl")

using ..homemade_conslaws: Grid1D, WallsBC, cell_centers, walls
using StaticArrays, LaTeXStrings

export animate_water


# Plots.scalefontsizes(1.5)
default(;fontfamily="Computer Modern", linewidth=2)
color_cycle = [Plots.RGB([30, 136, 229]./255...),
               Plots.RGB([255, 193, 7]./255...),
               Plots.RGB([216, 27, 96]./255...),]


function combine_matrices(U, V)
    M, N = size(U)
    combined = Matrix{eltype(U)}(undef, M, 2N)
    combined[:, 1:2:end] .= U
    combined[:, 2:2:end] .= V
    return combined
end


function animate_water_rec(H_left, H_right, UH_left, UH_right, cell_faces, t, walls, animation_duration)
    x = vec(stack([[cell_faces[i], cell_faces[i+1]] for i in eachindex(cell_faces)[1:end-1]]))
    H_rec = combine_matrices(H_left, H_right)
    UH_rec = combine_matrices(UH_left, UH_right)

    animate_water(H_rec, UH_rec, x, t, walls, animation_duration)
end

function animate_water_rec2(H_left, H_right, UH_left, UH_right, cell_faces, t, B, animation_duration)
    x = vec(stack([[cell_faces[i], cell_faces[i+1]] for i in eachindex(cell_faces)[1:end-1]]))
    H = combine_matrices(H_left, H_right)
    UH = combine_matrices(UH_left, UH_right)

    H_lim = calc_ylim(extrema(H)..., 0.1)
    UH_lim = calc_ylim(extrema(abs.(UH))..., 0.1)

    # l = [H[i] + H[i+1] for i in eachindex(H[end,1:2:end])] .* 0.5
    # k = [UH[i] + UH[i+1] for i in eachindex(UH[end,1:2:end])] .* 0.5
    # @show H[end,:]a
    # @show size(l)

    fps, skip_frames = calc_fps(t, animation_duration)

    anim = @animate for n in eachindex(t)[1:skip_frames:end]
        wateranim(x, H[n, :], UH[n, :], "H", H_lim, UH_lim)
        bathometryplot!(cell_faces, B)

        timeanim!(t[n], t[end], x[1], x[end], H_lim)
    end
    # for n in eachindex(t)[1:skip_frames:end]
    #     wateranim(x, H[n, :], UH[n, :], "H", H_lim, UH_lim)
    #     bathometryplot!(cell_faces, B)

    #     timeanim!(t[n], t[end], x[1], x[end], H_lim)
    # end

    gif(anim, fps=fps)
end

function animate_water(H, UH, x, t, walls, animation_duration)
    H_lim = calc_ylim(extrema(H)..., 0.1)
    UH_lim = calc_ylim(extrema(UH)..., 0.1)

    fps, skip_frames = calc_fps(t, animation_duration)

    anim = @animate for n in eachindex(t)[1:skip_frames:end]
        wateranim(x, H[n, :], UH[n, :], "H", H_lim, UH_lim)
        for wall in walls
            wallplot!(wall..., H_lim...)
        end
        timeanim!(t[n], t[end], x[1], x[end], H_lim)
    end

    gif(anim, fps=fps)
end

flatten = vec ∘ stack

function plot_water(H_left, H_right, UH_left, HU_right, cell_faces, B, fig_name=nothing)
    x = flatten([[cell_faces[i], cell_faces[i+1]] for i in eachindex(cell_faces)[1:end-1]])
    H = flatten([[H_left[i], H_right[i]] for i in eachindex(H_left)])
    # UH = flatten([[UH_left[i], HU_right[i]] for i in eachindex(UH_left)])
    UH = flatten([[0.5*(HU_right[i] + UH_left[i]), HU_right[i]] for i in eachindex(UH_left)])

    @assert size(H) == size(UH) == size(x) "$(size(H)) == $(size(UH)) == $(size(x))"

    println("$(size(H)) == $(size(UH)) == $(size(x))")

    wateranim(x, H, abs.(UH), "H", (-0.25, 2.5*1.1), (0., 2.4))
    # wateranim(x, H, abs.(UH), "H", (0., maximum(H)*1.1), (0., maximum(abs.(UH))))
    # wateranim(x, H, UH, "H", calc_ylim(extrema(H)..., 0.1), calc_ylim(extrema(UH)..., 0.1))
    bathometryplot!(cell_faces, B)
    if !isnothing(fig_name)
        savefig("$fig_name.pdf")
    end
end

function plot_water_walls_bc(H_left, H_right, UH_left, HU_right, cell_faces, wall_start, fig_name=nothing)
    x = flatten([[cell_faces[i], cell_faces[i+1]] for i in eachindex(cell_faces)[1:end-1]])
    H = flatten([[H_left[i], H_right[i]] for i in eachindex(H_left)])
    # UH = flatten([[UH_left[i], HU_right[i]] for i in eachindex(UH_left)])
    UH = flatten([[0.5*(HU_right[i] + UH_left[i]), 0.] for i in eachindex(UH_left)])

    @assert size(H) == size(UH) == size(x) "$(size(H)) == $(size(UH)) == $(size(x))"

    # wateranim(x[1:2:end], H_left, abs.(HU_right), "H", (-0.25, 2.5*1.1), (0., 1.4))
    wateranim(x, H, abs.(UH), "H", (-0.25, 2.5*1.1), (0., 1.4))
    # wateranim(x, H, abs.(UH), "H", (0., maximum(H)*1.1), (0., maximum(abs.(UH))))
    # wateranim(x, H, UH, "H", calc_ylim(extrema(H)..., 0.1), calc_ylim(extrema(UH)..., 0.1))
    bathometryplot!([cell_faces[1], wall_start, wall_start, cell_faces[end]], [0., 0., 10., 10.])
    if !isnothing(fig_name)
        savefig("$fig_name.pdf")
    end
end


function plot_time_step(t, dt, fig_name=nothing)
    plot(#t, dt,
        dt,
        # xlabel="Time " * L"(\mathrm{s})",
        xlabel="Iteration " * L"n",
        ylabel="Time-step size " * L"(\mathrm{s})",
        label=L"\Delta t",
        # dotted
        linestyle=:dash,
        marker=:circle,
        color=color_cycle[1],
        ylims=(0.5minimum(dt), maximum(dt)*1.1),
        yaxis=:log)
    if !isnothing(fig_name)
        savefig("$fig_name.pdf")
    end
end

function plot_time(t, fig_name=nothing)
    plot(t,
        xlabel="Iteration",
        ylabel="Time " * L"(\mathrm{s})",
        label=L"t",
        color=color_cycle[1],
    )
    if !isnothing(fig_name)
        savefig("$fig_name.pdf")
    end
end


function plot_error(grid_sizes, H_error, HU_error, error_type, name=nothing)
    plot(grid_sizes, H_error, 
        xlabel="Grid resolution " * L"N",
        ylabel=error_type,
        linestyle=:dash,
        marker=:circle,
        markercolor=color_cycle[1],
        label=L"h"*" error "*L"(\mathrm{m})",
        legend=:left,
        linecolor=color_cycle[1],
        xaxis=:log, yaxis=:log)
    plot!(grid_sizes, HU_error,
        label=L"hu"*" error "*L"(\mathrm{m}^2\mathrm{s}^{-1})",
        linestyle=:dash,
        marker=:circle,
        linecolor=color_cycle[2],
        markercolor=color_cycle[2])
    if !isnothing(name)
        savefig("$name.pdf")
    end
end

function plot_error_wall(wall_heights, H_error, HU_error, error_type, name=nothing)
    plot(wall_heights, H_error, 
        xlabel="Wall heights " * L"h_W\, (\mathrm{m})",
        ylabel=error_type,
        linestyle=:dash,
        marker=:circle,
        markercolor=color_cycle[1],
        label=L"h"*" error "*L"(\mathrm{m})",
        legend=:right,
        linecolor=color_cycle[1],
        xaxis=:log, yaxis=:log
    )
    plot!(wall_heights, HU_error,
        label=L"hu"*" error "*L"(\mathrm{m}^2\mathrm{s}^{-1})",
        linestyle=:dash,
        marker=:circle,
        linecolor=color_cycle[2],
        markercolor=color_cycle[2])
    if !isnothing(name)
        savefig("$name.pdf")
    end
end


function plot_timesteps(wall_heights, dt, dt_ref)
    plot(wall_heights, dt, 
        ylim=(0.5minimum(dt), maximum(dt)*1.1),
        xlabel="Wall height " * L"(\mathrm{m})", 
        ylabel="Time-step size " * L"(\mathrm{s})", 
        label=L"\Delta t_{TOP}",
        linecolor=color_cycle[1])
    plot!(wall_heights[[1,end]], [dt_ref, dt_ref], label=L"\Delta t_{BC}",
        linecolor=color_cycle[2])
    savefig("timestep-wall_height.pdf")
end


function plot_momenta(UH, N)
    plot(N, UH, label="Max absolute momentum")
end


function plot_momentum(UH_left, UH_right, cell_faces, B)
    x = flatten([[cell_faces[i], cell_faces[i+1]] for i in eachindex(cell_faces)[1:end-1]])
    UH = flatten([[UH_left[i], UH_right[i]] for i in eachindex(UH_left)])

    wateranim(x, UH, abs.(UH), "UH", calc_ylim(extrema(UH)..., 0.1), (0., maximum(abs.(UH))))
end

function animate_water(H, UH, grid::Grid1D, t, animation_duration=t[end])
    animate_water(H, UH, cell_centers(grid), t, [], animation_duration)
end

function animate_water(H, UH, grid::Grid1D{WallsBC{Float}}, t, animation_duration=t[end]) where Float
    animate_water(H, UH, cell_centers(grid), t, walls(grid), animation_duration)
end
