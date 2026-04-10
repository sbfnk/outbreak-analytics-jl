### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4c5d
begin
    import Pkg
    Pkg.develop(path=joinpath(@__DIR__, "..", "EpiParameters"))
    using EpiParameters
    using Distributions
    using CairoMakie
end

# ╔═╡ 1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c6d
md"""
# Access Epidemiological Delay Distributions

Infectious diseases follow an **infection cycle** with distinct time periods:
exposure, presymptomatic infection, symptomatic infection, and recovery (or death).
The durations of these periods — and the delays between observable events like
symptom onset, hospitalisation, and death — are critical for informing disease
prevention and control measures. For example, knowing how long someone is infectious
before showing symptoms determines whether contact tracing can outpace transmission.

Early in an epidemic, these parameters are uncertain and evolving. Manually searching
the literature for the latest estimates is slow, error-prone, and difficult to
reproduce. Projects like **EpiParameters.jl** address this by building programmatic
access to curated catalogues of parameters synthesised from published studies, so
that analysts can query and use them directly in their pipelines.

## Learning objectives

- Understand why epidemiological delay distributions matter for outbreak analysis
- Query the epidemiological parameter database using `epiparameter()`
- Extract distribution objects for use in analysis
- Visualise and compare delay distributions across diseases
"""

# ╔═╡ 2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6c7d
md"""
## Generation time vs serial interval

Two commonly confused delay distributions:

- **Generation time**: time from infection of a primary case to infection of a secondary case (infection → infection)
- **Serial interval**: time from symptom onset of a primary case to symptom onset of a secondary case (symptoms → symptoms)

The generation time, jointly with the reproduction number *R*, determines how fast
an epidemic grows. A larger *R* and/or a shorter generation time both lead to faster
epidemic growth. In practice, we often use the **serial interval** as a proxy for
the generation time, because symptom onset dates are much easier to observe than
the moment of infection itself.

However, the serial interval can be **negative** for diseases with pre-symptomatic
transmission (e.g. COVID-19) — if a secondary case develops symptoms before the
primary case does — while the generation time is always positive.

Importantly, not all transmission pairs have the same delay. If we measured the
serial interval for many case pairs, we would see considerable variation. This is
why we need **full distributions** rather than just mean values: the shape and
spread of the distribution capture the range of plausible delays and are essential
for accurate modelling.
"""

# ╔═╡ 3a4b5c6d-7e8f-9a0b-1c2d-3e4f5a6b7c8d
md"""
## Exploring the database

The `epiparameter()` function is the main entry point for querying the database.
Called without arguments, it returns all entries; with keyword arguments, it acts
as a filter. Not all entries in the database have **fitted distributions** — some
only record summary statistics (e.g. mean and standard deviation) reported in the
original paper, without a parameterised probability distribution. When we need a
distribution object for modelling, we filter for entries where `distribution` is
not `nothing`.

Let's start by seeing what diseases and parameter types are available:
"""

# ╔═╡ 4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8c9d
list_diseases()

# ╔═╡ 5a6b7c8d-9e0f-1a2b-3c4d-5e6f7a8b9c0d
list_parameters()

# ╔═╡ 6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0c1d
md"""
## Querying for a specific disease

Let's find all COVID-19 incubation periods. The database may return multiple
entries from different studies. We then filter for those that have a fitted
distribution, since only these can be used directly as `Distributions.jl` objects
in downstream analysis:
"""

# ╔═╡ 7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1c2d
covid_ip = epiparameter(disease="COVID-19", epi_name="incubation period")

# ╔═╡ 8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2c3d
covid_ip_fitted = filter(p -> !isnothing(p.distribution), covid_ip)

# ╔═╡ 9a0b1c2d-3e4f-5a6b-7c8d-9e0f1a2b3c4d
md"""
## Visualising a distribution

The distribution is a standard `Distributions.jl` object — we can plot its
probability density function:
"""

# ╔═╡ aa1b2c3d-4e5f-6a7b-8c9d-ae1f2a3b4c5d
best = covid_ip_fitted[1]

# ╔═╡ ba2b3c4d-5e6f-7a8b-9c0d-be2f3a4b5c6d
(mean=mean(best.distribution), std=std(best.distribution), median=median(best.distribution))

# ╔═╡ ca3b4c5d-6e7f-8a9b-0c1d-ce3f4a5b6c7d
let d = best.distribution
    x = range(0, quantile(d, 0.99); length=200)
    fig = Figure(size=(600, 350))
    ax = Axis(fig[1, 1];
              xlabel="Days", ylabel="Density",
              title="COVID-19 incubation period ($(typeof(d).name.name))")
    lines!(ax, x, pdf.(d, x); linewidth=2.5)
    vlines!(ax, [mean(d)]; color=:red, linestyle=:dash, linewidth=1.5,
            label="Mean = $(round(mean(d); digits=1))d")
    axislegend(ax; position=:rt)
    fig
end

# ╔═╡ da4b5c6d-7e8f-9a0b-1c2d-de4f5a6b7c8d
md"""
## Comparing distributions across diseases

Let's compare incubation periods for different diseases:
"""

# ╔═╡ ea5b6c7d-8e9f-0a1b-2c3d-ee5f6a7b8c9d
let
    diseases = ["COVID-19", "Influenza", "SARS", "Ebola Virus Disease", "Mpox"]
    fig = Figure(size=(700, 400))
    ax = Axis(fig[1, 1];
              xlabel="Days", ylabel="Density",
              title="Incubation periods by disease")
    x = range(0, 25; length=300)

    for disease in diseases
        results = epiparameter(disease=disease, epi_name="incubation period")
        fitted = filter(p -> !isnothing(p.distribution), results)
        isempty(fitted) && continue
        d = fitted[1].distribution
        lines!(ax, x, pdf.(d, x); linewidth=2,
               label="$(disease) (μ=$(round(mean(d); digits=1))d)")
    end

    axislegend(ax; position=:rt)
    fig
end

# ╔═╡ ab12cd34-ef56-7890-ab12-cd34ef567890
md"""
The comparison above highlights how different diseases can have markedly different
incubation periods. COVID-19 and influenza tend to have relatively short incubation
periods (a few days), while diseases like Ebola and Mpox have longer ones. These
differences have direct implications for surveillance and control: shorter incubation
periods mean less time to identify and isolate cases before they become symptomatic
and potentially infectious.
"""

# ╔═╡ fa6b7c8d-9e0f-1a2b-3c4d-fe6f7a8b9c0d
md"""
## Citation and metadata

Always check the source when selecting a parameter. The reliability of an estimate
depends on factors like **sample size**, the **region** where data were collected,
and the **inference method** used. A parameter estimated from a handful of cases in
a single hospital may not generalise to a different setting or epidemic phase.
"""

# ╔═╡ 0b7b8c9d-0e1f-2a3b-4c5d-0f7f8a9b0c1d
best.citation["year"], best.citation["doi"]

# ╔═╡ 1b8b9c0d-1e2f-3a4b-5c6d-1f8f9a0b1c2d
best.metadata

# ╔═╡ 2b9b0c1d-2e3f-4a5b-6c7d-2f9f0a1b2c3d
md"""
## Onset-to-death distributions

These are critical for CFR estimation. Let's visualise the Ebola onset-to-death
distribution:
"""

# ╔═╡ 3ba01c2d-3e4f-5a6b-7c8d-3fa01a2b3c4d
let
    results = epiparameter(disease="Ebola", epi_name="onset to death")
    fitted = filter(p -> !isnothing(p.distribution), results)
    d = fitted[1].distribution
    x = range(0, quantile(d, 0.99); length=200)

    fig = Figure(size=(600, 350))
    ax = Axis(fig[1, 1];
              xlabel="Days from onset", ylabel="Density",
              title="Ebola onset-to-death delay ($(typeof(d).name.name))")
    band!(ax, x, zeros(length(x)), pdf.(d, x); color=(:steelblue, 0.3))
    lines!(ax, x, pdf.(d, x); linewidth=2.5, color=:steelblue)
    vlines!(ax, [mean(d)]; color=:red, linestyle=:dash, linewidth=1.5,
            label="Mean = $(round(mean(d); digits=1))d")
    axislegend(ax; position=:rt)
    fig
end

# ╔═╡ 4ba12c3d-4e5f-6a7b-8c9d-4fa12a3b4c5d
md"""
## Challenge: Ebola parameters

1. How many Ebola parameters are in the database? How many have fitted distributions?
2. Find the serial interval for Ebola — what is the mean?
3. Plot the onset-to-death distribution alongside the incubation period.

**Bonus question**: Look up the serial interval for both COVID-19 and SARS.
If both diseases had a similar reproduction number *R*, which would be harder to
control and why? *(Hint: think about how the serial interval relates to the speed
of epidemic growth and the time available for interventions.)*
"""

# ╔═╡ 5ba23c4d-5e6f-7a8b-9c0d-5fa23a4b5c6d
# Try it yourself!
ebola_all = epiparameter(disease="Ebola")

# ╔═╡ 6ba34c5d-6e7f-8a9b-0c1d-6fa34a5b6c7d
md"""
## Key points

- **EpiParameters.jl** provides access to a curated database of epidemiological
  delay distributions from the literature
- Use `epiparameter(disease=..., epi_name=...)` to query by disease and parameter type
- Entries with fitted parameters return `Distributions.jl` objects that can be
  plotted and used directly in analysis pipelines
- Use `list_diseases()` and `list_parameters()` to explore what's available
- Always check the citation and metadata (sample size, region, inference method)
  when selecting a parameter for your analysis
"""

# ╔═╡ Cell order:
# ╟─1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c6d
# ╠═0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4c5d
# ╟─2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6c7d
# ╟─3a4b5c6d-7e8f-9a0b-1c2d-3e4f5a6b7c8d
# ╠═4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8c9d
# ╠═5a6b7c8d-9e0f-1a2b-3c4d-5e6f7a8b9c0d
# ╟─6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0c1d
# ╠═7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1c2d
# ╠═8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2c3d
# ╟─9a0b1c2d-3e4f-5a6b-7c8d-9e0f1a2b3c4d
# ╠═aa1b2c3d-4e5f-6a7b-8c9d-ae1f2a3b4c5d
# ╠═ba2b3c4d-5e6f-7a8b-9c0d-be2f3a4b5c6d
# ╠═ca3b4c5d-6e7f-8a9b-0c1d-ce3f4a5b6c7d
# ╟─da4b5c6d-7e8f-9a0b-1c2d-de4f5a6b7c8d
# ╠═ea5b6c7d-8e9f-0a1b-2c3d-ee5f6a7b8c9d
# ╟─ab12cd34-ef56-7890-ab12-cd34ef567890
# ╟─fa6b7c8d-9e0f-1a2b-3c4d-fe6f7a8b9c0d
# ╠═0b7b8c9d-0e1f-2a3b-4c5d-0f7f8a9b0c1d
# ╠═1b8b9c0d-1e2f-3a4b-5c6d-1f8f9a0b1c2d
# ╟─2b9b0c1d-2e3f-4a5b-6c7d-2f9f0a1b2c3d
# ╠═3ba01c2d-3e4f-5a6b-7c8d-3fa01a2b3c4d
# ╟─4ba12c3d-4e5f-6a7b-8c9d-4fa12a3b4c5d
# ╠═5ba23c4d-5e6f-7a8b-9c0d-5fa23a4b5c6d
# ╟─6ba34c5d-6e7f-8a9b-0c1d-6fa34a5b6c7d
