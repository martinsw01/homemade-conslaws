@userplot WaterAnim
@recipe function f(anim::WaterAnim)
    x, H, UH, name, H_ylim, U_ylim = anim.args
    
    ylim --> H_ylim
    legend --> :topleft
    # legendcolumns --> length(names)
    label --> name

    seriestype --> :path
    fillrange --> H_ylim[1]
    fill_z --> abs.(UH)./H

    # Define the color gradient
    c --> :blues
    clims --> U_ylim

    [x], [H]
end