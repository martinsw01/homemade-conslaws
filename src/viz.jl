module Viz

using Plots, LaTeXStrings

export plot_solution, animate_solution, plot_matrix

function plot_solution(u, x, t)
    heatmap(x, t, u, color = :greys, clim=(-1.3,1.3),
            xlabel="x", ylabel="t", show=true)
end


@userplot SolutionAnim
@recipe function f(anim::SolutionAnim)
    x, solutions, names, ylim = anim.args
    
    ylim --> ylim
    legend --> :topleft
    # legendcolumns --> length(names)
    label --> names

    [x], [solutions]
end

@userplot TimeAnim
@recipe function f(tp::TimeAnim)
    t_n, T, x_L, x_R, (_, ymax) = tp.args
    ymax

    label --> ["" "\$t=$(round(t_n, digits=2))\$"]
    linewidth --> 10
    color --> [:lightgrey :grey]

    [[x_L, x_R], [x_L, (x_R - x_L) * t_n / T + x_L]],
    [[ymax, ymax], [ymax, ymax]]
end

function calc_ylim(min, max, padding)
    return min - padding * (max - min),
           max + padding * (max - min)
end

function calc_fps(t)
    min_time = 1
    fps_real_time = length(t) / max(t[end], min_time)

    return min(50, fps_real_time) # Max supported fps is 50 on most brwosers
end

function animate_solution(solutions, names, x, t)
    ylim = calc_ylim(extrema(solutions[end])..., 0.1)

    anim = @animate for (n, t_n) in enumerate(t)
        solutionanim(x,
                     [u[n, :] for u in solutions],
                     names, ylim)
        timeanim!(t_n, t[end], x[1], x[end], ylim)
    
    end

    gif(anim, fps=calc_fps(t))
end

function animate_solution(u_approx, u_exact::Function, x, t)
    animate_solution((u_approx, u_exact.(x', t)),
                     [L"U_\operatorname{approx}" L"U_\operatorname{exact}"],
                     x, t)
end

function animate_solution(u_approx, x, t)
    animate_solution((u_approx,),
                     (:none,),
                     x, t)
end

function plot_matrix(A, N, M; plot_kargs...)
    B = A[end:-1:1, 1:end]

    hm = heatmap(B, plot_kargs...)
    vline!(N+0.5:N:N*M, color=:black, label="")
    hline!(N+0.5:N:N*M, color=:black, label="")
    return hm
end

end