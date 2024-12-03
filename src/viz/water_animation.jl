include("shallow_water_recipe.jl")
include("wall_recipe.jl")

using ..homemade_conslaws: Grid1D, WallsBC, cell_centers, walls

export animate_water


function _animate_water(H, UH, x, t, walls, animation_duration)
    H_lim = calc_ylim(extrema(H)..., 0.1)
    U_lim = calc_ylim(extrema(abs.(UH ./ H))..., 0.1)

    fps, skip_frames = calc_fps(t, animation_duration)

    anim = @animate for n in eachindex(t)[1:skip_frames:end]
        wateranim(x, H[n, :], UH[n, :], "H", H_lim, U_lim)
        for wall in walls
            wallplot!(wall..., H_lim...)
        end
        timeanim!(t[n], t[end], x[1], x[end], H_lim)
    end

    gif(anim, fps=fps)
end


function animate_water(H, UH, grid::Grid1D, t, animation_duration=t[end])
    _animate_water(H, UH, cell_centers(grid), t, [], animation_duration)
end

function animate_water(H, UH, grid::Grid1D{WallsBC{Float}}, t, animation_duration=t[end]) where Float
    _animate_water(H, UH, cell_centers(grid), t, walls(grid), animation_duration)
end
