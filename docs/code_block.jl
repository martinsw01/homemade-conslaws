using Highlights, Highlights.Lexers
import Base: replace

function _convert_to_pygments_classes(html)
    substitutes = [
        "hljl-oB" => "o",
        "hljl-nfB" => "mf",
        "class='hljl-" => "class='",
    ]
    reduce(substitutes; init=html) do html, (pattern, replacement)
        replace(html, pattern => replacement)
    end
end

_wrap_in_div(html) = "<div class='highlight'><pre><span></span><code>$(html[20:end-7])</code></pre></div>"

Base.replace(f::Function, A, r::Regex) = replace(A, r => f)

function _highlight_code(code)
    io = IOBuffer()
    highlight(io, MIME("text/html"), code, Lexers.JuliaLexer)
    String(take!(io))
end

function _hljl_to_pygments(html)
    html |> _convert_to_pygments_classes |> _wrap_in_div
end

function _highlight_code_blocks(content::String)
    code_block_regex = r"```julia\n(.*?)\n```"s

    replace(content, code_block_regex) do match
        code = match[10:end-4]
        code |> _highlight_code |> _hljl_to_pygments
    end
end

function _process_markdown_files!(f!, folder)
    for (root, dirs, files) in walkdir(folder)
        for file in files
            if endswith(file, ".md")
                f!(file, root)
            end
        end
    end
end

format_codeblocks_in_folder!(folder::String) = _process_markdown_files!(folder) do mdfile, root
    filepath = joinpath(root, mdfile)
    content = read(filepath, String)
    highlighted_content = _highlight_code_blocks(content)
    open(filepath, "w") do io
        write(io, highlighted_content)
    end
end