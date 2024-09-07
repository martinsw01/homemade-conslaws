module TempJmd

using Weave, Plots

export weave_all_files, temporarily_weave_all_files

function Weave.WeavePlots.add_plots_figure(report::Weave.Report, plot::Plots.AnimatedGif, ext)
    chunk = report.cur_chunk
    full_name, rel_name = Weave.get_figname(report, chunk, ext = ext)

    # An `AnimatedGif` has been saved somewhere temporarily, so make a copy to `full_name`.
    cp(plot.filename, full_name; force = true)
    push!(report.figures, rel_name)
    report.fignum += 1
    return full_name
end

function Base.display(report::Weave.Report, m::MIME"text/plain", plot::Plots.AnimatedGif)
    Weave.WeavePlots.add_plots_figure(report, plot, ".gif")
end

"""
Finds paths to all files in `dir` with a given `suffix` relative to `dir`.
"""
function find_files(dir, suffix = ".jmd")
    return [joinpath(root, f)
            for (root, dirs, files) in walkdir(dir)
                for f in files
                    if endswith(f, suffix)]
end


function weave_all_files(dir)
    files = find_files(dir)

    generated_dir = joinpath(dir, "generated")

    for f in files
        path_relative_to_dir = relpath(dirname(f), dir)
        out_path = joinpath(generated_dir,path_relative_to_dir)
        cache_path = joinpath("cache", path_relative_to_dir)
        mkpath(out_path)
        weave(f, doctype="github", out_path=out_path, cache=:user, cache_path=cache_path)
    end

    return relpath(generated_dir, dir)
end

end # module