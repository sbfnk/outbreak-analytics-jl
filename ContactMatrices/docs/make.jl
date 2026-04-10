using Documenter
using ContactMatrices

makedocs(;
    sitename = "ContactMatrices.jl",
    authors = "Sebastian Funk and contributors",
    modules = [ContactMatrices],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)
