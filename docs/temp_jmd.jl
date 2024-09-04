module TempJmd

using Weave

export weave_all_files

function find_files(dir, suffix = ".jmd")
    files = []
    for f in readdir(dir)
        if isfile(joinpath(dir, f)) && endswith(f, suffix)
            push!(files, joinpath(dir, f))
        elseif isdir(joinpath(dir, f))
            files = vcat(files, find_files(joinpath(dir, f)))
        end
    end
    return files
end

function weave_all_files(dir)
    files = find_files(dir)

    generated_dir = dir*"/generated"

    for f in files
        out_path = generated_dir*dirname(f)[length(dir)+1:end]
        mkpath(out_path)
        weave(f, doctype="github", out_path=out_path)
    end

    return generated_dir
end

end # module