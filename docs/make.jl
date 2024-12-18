using Documenter, DocumenterMarkdown

makedocs(
    sitename="Homemade Conslaws";
    # draft = true, # for livereload. Disables code execution in @example blocks
    format = Markdown(),
    pages = [
        "Conservation Laws" => "index.md",
        "Theory" => [
            "Linear transport equations" => "theory/linear_transport_eqs.md",
            "Scalar conservation laws" => "theory/scalar_cons_laws.md",
            "Finite volume schemes" => "theory/finite_volume_schemes.md",
            "Second order schemes" => "theory/2nd_order_schemes.md",
            "Nonlinear hyperbolic systems" => "theory/nonlinear_hyperbolic_systems.md"
        ]
    ]
)

if "deploy" in ARGS
    deploydocs(
        repo="gihub.com/martinsw01/homemade-conslaws.jl.git",
        deps=Deps.pip("mkdocs", "mkdocs-material"),
        make=() -> run(`mkdocs build`),
        target="site"
    )
end