using Documenter
using FinalSize

makedocs(;
    sitename = "FinalSize.jl",
    authors = "Sebastian Funk and contributors",
    modules = [FinalSize],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)
