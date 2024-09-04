module MathEngine

using Documenter

export create_mathengine

macros() = Dict(
    :R => raw"\mathbb{R}",
    :bm => [raw"{\boldsymbol{#1}}", 1],
    :x => raw"\bm{\text{x}}",
    :s => raw"\sigma",
    :div => [raw"\operatorname{div}\left(#1\right)", 1],
    )


function create_mathengine()
    return Documenter.MathJax3(Dict(
        :loader => Dict("load" => ["[tex]/physics"]),
        :tex => Dict(
            :inlineMath => [["\$","\$"], ["\\(","\\)"]],
            :tags       => "ams",
            :packages   => ["base", "ams", "autoload", "configmacros", "physics"],
            :macros     => macros(),
        ),
    ))
end

end