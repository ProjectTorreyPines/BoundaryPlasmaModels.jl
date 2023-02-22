using BoundaryPlasmaModels
using Documenter

DocMeta.setdocmeta!(BoundaryPlasmaModels, :DocTestSetup, :(using BoundaryPlasmaModels); recursive=true)

makedocs(;
    modules=[BoundaryPlasmaModels],
    authors="Jerome Guterl",
    repo="https://github.com/jguterl/BoundaryPlasmaModels.jl/blob/{commit}{path}#{line}",
    sitename="BoundaryPlasmaModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jguterl.github.io/BoundaryPlasmaModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jguterl/BoundaryPlasmaModels.jl",
    devbranch="main",
)
