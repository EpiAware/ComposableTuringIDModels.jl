using EpiAwarePrototype
using Documenter

DocMeta.setdocmeta!(
    EpiAwarePrototype, :DocTestSetup, :(using EpiAwarePrototype, Distributions);
    recursive = true)

makedocs(;
    modules = [EpiAwarePrototype],
    authors = "Sam Abbott",
    sitename = "EpiAwarePrototype.jl",
    format = Documenter.HTML(;
        canonical = "https://epiawareprototype.epiaware.org",
        edit_link = "main",
        assets = String[]),
    pages = [
        "Home" => "index.md",
        "Composable design" => "design.md",
        "API reference" => "api.md"
    ])

deploydocs(;
    repo = "github.com/EpiAware/EpiAwarePrototype.jl",
    devbranch = "main")
