using Documenter
using CFR

makedocs(;
    sitename = "CFR.jl",
    authors = "Sebastian Funk and contributors",
    modules = [CFR],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)
