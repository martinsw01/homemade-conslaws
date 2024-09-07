module Viz

using Plots, LaTeXStrings

export plot_solution, animate_solution, plot_matrix

function plot_solution(u, x, t)
    heatmap(x, t, u, color = :greys, clim=(-1.3,1.3),
            xlabel="x", ylabel="t", show=true)
end


@userplot SolutionAnim
@recipe function f(anim::SolutionAnim)
    x, t_n, solutions, names, ylim = anim.args
    
    title --> "t=$(round(t_n, digits=2))"
    ylim --> ylim
    legend --> :topleft
    # legendcolumns --> length(names)
    label --> names

    [x], [solutions]
end

@userplot CirclePlot
@recipe function f(cp::CirclePlot)
    x, y, i = cp.args
    n = length(x)
    inds = circshift(1:n, 1 - i)
    linewidth --> range(0, 10, length = n)
    seriesalpha --> range(0, 1, length = n)
    aspect_ratio --> 1
    label --> false
    x[inds], y[inds]
end

function calc_ylim(min, max, padding)
    return min - padding * (max - min),
           max + padding * (max - min)
end

function animate_solution(solutions, names, x, t)
    ylim = calc_ylim(extrema(solutions[end])..., 0.1)

    fps = length(t) * min(0.5, 1 / t[end])

    anim = @animate for (n, t_n) in enumerate(t)
        solutionanim(x, t_n,
                     [u[:, n] for u in solutions],
                     names, ylim)
    
    end

    gif(anim, fps=fps)
end

function animate_solution(u_approx, u_exact::Function, x, t)
    animate_solution((u_approx, u_exact.(x, t')),
                     [L"U_\text{approx}" L"U_\text{exact}"],
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