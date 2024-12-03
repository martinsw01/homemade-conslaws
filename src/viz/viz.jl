module Viz

export plot_solution, animate_solution, animate_solutions

using Plots, LaTeXStrings
using ..homemade_conslaws: Grid1D
include("solution_recipe.jl")
include("time_recipe.jl")
include("water_animation.jl")

animate_solutions(solutions, names, grid::Grid1D, t, T=t[end]) = animate_solutions(solutions, names, cell_centers(grid), t, T)

animate_solution(solution, name, grid::Grid1D, t, animation_duration=t[end]) = animate_solution(solution, name, cell_centers(grid), t, animation_duration)

animate_solution(approximation, exact::Function, grid::Grid1D, t, animation_duration=t[end]) = animate_solution(approximation, exact, cell_centers(grid), t, animation_duration)


function plot_solution(u, x, t)
    heatmap(x, t, u, color = :greys, clim=(-1.3,1.3),
            xlabel="x", ylabel="t", show=true)
end


function calc_ylim(min, max, padding)
    return min - padding * (max - min),
           max + padding * (max - min)
end

function calc_fps(frames, duration)
    fps_real_time = length(frames) / duration

    if fps_real_time > 50 # Max supported fps is 50 on most brwosers
        skip_frames = Int64(fps_real_time ÷ 50)
        return 50, skip_frames
    else
        return fps_real_time, 1
    end
end


function _animate_solutions(update_animation, solutions, names, x, t, anim_duration=t[end])
    ylim = calc_ylim(extrema(solutions[end])..., 0.1)

    fps, skip_frames = calc_fps(t, anim_duration)

    anim = @animate for n in eachindex(t)[1:skip_frames:end]
        update_animation([u[n, :] for u in solutions], names, x, t, n, ylim)
    end

    gif(anim, fps=fps)
end


function animate_solutions(solutions, names, x, t, anim_duration=t[end])
    _animate_solutions(solutions, names, x, t, anim_duration) do solutions, names, x, t, n, ylim
        solutionanim(x, solutions, names, ylim)
        timeanim!(t[n], t[end], x[1], x[end], ylim)
    end
end

function animate_solution(u_approx, u_exact::Function, x, t, animation_duration=t[end])
    animate_solutions((u_approx, u_exact.(x', t)),
                      [L"U_\operatorname{approx}" L"U_\operatorname{exact}"],
                      x, t, animation_duration)
end

function animate_solution(u_approx, name::Union{AbstractString, Symbol}, x, t, animation_duration=t[end])
    animate_solutions((u_approx,),
                      name,
                      x, t, animation_duration)
end

animate_solution(u_approx, x, t, animation_duration=t[end]) = animate_solution(u_approx, :none, x, t, animation_duration)

end