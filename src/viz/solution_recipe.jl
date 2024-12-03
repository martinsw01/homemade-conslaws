@userplot SolutionAnim
@recipe function f(anim::SolutionAnim)
    x, solutions, names, ylim = anim.args
    
    ylim --> ylim
    legend --> :topleft
    # legendcolumns --> length(names)
    label --> names

    [x], [solutions]
end