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
            "Linear transport equations" => "$generated_dir/theory/linear_transport_eqs.md",
            "Scalar conservation laws" => "theory/scalar_cons_laws.md",
        ]
    ]
)

if "deploy" in ARGS
    deploydocs(repo="gihub.com/martinsw01/homemade-conslaws.jl.git")
end