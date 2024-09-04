using Documenter

include("mathengine.jl")
include("temp_jmd.jl")

generated_dir = TempJmd.weave_all_files("docs/src")

makedocs(
    sitename="Homemade Conslaws";
    format = Documenter.HTML(; mathengine=MathEngine.create_mathengine()),
    pages = [
        "Conservation Laws" => "index.md",
        "Theory" => [
            "Linear transport equations" => "theory/linear_transport_eqs.md",
        ]
    ]
)

deploydocs(repo="github.com/martinsw01/homemade-conslaws.jl.git")

rm(generated_dir; recursive=true)