using Dmk
using Documenter

DocMeta.setdocmeta!(Dmk, :DocTestSetup, :(using Dmk); recursive=true)

makedocs(;
    modules=[Dmk],
    authors="Timo Betcke <timo.betcke@gmail.com> and contributors",
    sitename="Dmk.jl",
    format=Documenter.HTML(;
        canonical="https://tbetcke.github.io/Dmk.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tbetcke/Dmk.jl",
    devbranch="main",
)
