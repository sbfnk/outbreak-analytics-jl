using Documenter
using Epidemics

makedocs(;
    sitename = "Epidemics.jl",
    authors = "Sebastian Funk and contributors",
    modules = [Epidemics],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)
