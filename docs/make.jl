using DmspElfinConjunction
using Documenter

list_pages(dir) = map(f -> "$dir/$f", readdir(joinpath(@__DIR__, "src", dir)))

makedocs(;
    modules = [DmspElfinConjunction],
    sitename = "DmspElfinConjunction.jl",
    format = Documenter.HTML(;
        canonical = "https://beforerr.github.io/DmspElfinConjunction.jl",
    ),
    pages = [
        "Home" => "index.md",
        "demo.md"
    ],
)

deploydocs(;
    repo = "github.com/beforerr/DmspElfinConjunction.jl",
    push_preview = true
)
