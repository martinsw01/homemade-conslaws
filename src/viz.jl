module Viz

using Plots

export plot_solution, animate_solution, plot_matrix

function plot_solution(u, x, t)
    heatmap(x, t, u, color = :greys, clim=(-1.3,1.3),
            xlabel="x", ylabel="t", show=true)
end

function animate_solution(u, N, M, name="animation.gif")
    x = (1:N)/N
    t = (1:M)/M

    anim = @animate for n in 1:M
        plot(x, u[1:N, n],
             ylim=(-1.3, 1.3),
             label="")
    end

    return gif(anim, name, fps=15)
end

function plot_matrix(A, N, M)
    B = A[end:-1:1, 1:end]

    hm = heatmap(B, show=true)
    vline!(N+0.5:N:N*M, color=:black, label="")
    hline!(N+0.5:N:N*M, color=:black, label="")
    return hm
end

end