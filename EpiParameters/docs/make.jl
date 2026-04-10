using Documenter
using EpiParameters

makedocs(;
    sitename = "EpiParameters.jl",
    authors = "Sebastian Funk and contributors",
    modules = [EpiParameters],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)
