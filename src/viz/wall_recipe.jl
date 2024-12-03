@userplot WallPlot
@recipe function f(wp::WallPlot)
    start, stop, y_min, y_max = wp.args
    label --> :none
    seriestype --> :shape
    # fillalpha --> 0.5
    fillcolor --> :black
    [[start, start, stop, stop, start]], [[y_min, y_max, y_max, y_min, y_min]]
end