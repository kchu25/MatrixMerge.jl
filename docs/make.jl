using MatrixMerge
using Documenter

DocMeta.setdocmeta!(MatrixMerge, :DocTestSetup, :(using MatrixMerge); recursive=true)

makedocs(;
    modules=[MatrixMerge],
    authors="Shane Kuei Hsien Chu (skchu@wustl.edu)",
    repo="https://github.com/kchu25/MatrixMerge.jl/blob/{commit}{path}#{line}",
    sitename="MatrixMerge.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kchu25.github.io/MatrixMerge.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/MatrixMerge.jl",
    devbranch="main",
)
