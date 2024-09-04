using Documenter

include("mathengine.jl")
include("temp_jmd.jl")

generated_dir = TempJmd.weave_all_files("docs/src")

makedocs(
    sitename="Homemade Conslaws";
    format = Documenter.HTML(; mathengine=MathEngone.create_mathengine()),
    pages = [
        "Conservation Laws" => "index.md",
    ]
)

rm(generated_dir; recursive=true)