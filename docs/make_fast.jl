## Build docs excluding EpiNow2 tutorials (which require MCMC and take ~30min).
## Use `make.jl` for the full build.

using Documenter
using DocumenterVitepress

makedocs(;
    sitename = "Outbreak Analytics in Julia",
    authors = "Sebastian Funk and contributors",
    remotes = nothing,
    warnonly = true,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/sbfnk/outbreak-analytics-jl",
        devbranch = "main",
        devurl = "dev",
    ),
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Early tasks" => [
                "Read case data" => "tutorials/read-case-data.md",
                "Clean case data" => "tutorials/clean-case-data.md",
                "Validate case data" => "tutorials/validate-case-data.md",
                "Aggregate and visualise" => "tutorials/aggregate-visualise.md",
            ],
            "Middle tasks" => [
                "Epidemiological parameters" => "tutorials/epidemiological-parameters.md",
                "Outbreak severity" => "tutorials/outbreak-severity.md",
                "Account for superspreading" => "tutorials/superspreading.md",
                "Transmission chains" => "tutorials/transmission-chains.md",
            ],
            "Late tasks" => [
                "Contact matrices" => "tutorials/contact-matrices.md",
                "Simulating transmission" => "tutorials/simulating-transmission.md",
                "Choosing an appropriate model" => "tutorials/choosing-model.md",
                "Modelling disease burden" => "tutorials/disease-burden.md",
            ],
        ],
    ],
)
