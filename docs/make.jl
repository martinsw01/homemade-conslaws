using Documenter

include("mathengine.jl")

makedocs(
    sitename="Homemade Conslaws";
    format = Documenter.HTML(; mathengine=MathEngone.create_mathengine()),
    pages = [
        "Conservation Laws" => "index.md",
    ]
)