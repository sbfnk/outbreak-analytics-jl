using Documenter
using SocialMixr

makedocs(;
    sitename = "SocialMixr.jl",
    authors = "Sebastian Funk and contributors",
    modules = [SocialMixr],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)
