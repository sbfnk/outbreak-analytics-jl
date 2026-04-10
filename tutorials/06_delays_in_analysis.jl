### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4d02
begin
    import Pkg
    Pkg.develop(path=joinpath(@__DIR__, "..", "EpiParameters"))
    Pkg.develop(path=expanduser("~/code/EpiNow2.jl"))
    using EpiNow2
    using EpiParameters
    using Distributions
    using DataFrames
    using Dates
    using CairoMakie
end

# ╔═╡ 1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5d02
md"""
# Use Delay Distributions in Analysis

Epidemiological delay distributions — incubation period, generation time,
serial interval, reporting delay — are essential inputs for transmission
estimation, forecasting, and severity assessment.

This tutorial integrates the previous two: we use the parameter database from
Tutorial 01 together with the $R_t$ estimation framework from Tutorial 05. The
key workflow is to retrieve parameters from the literature, convert them to an
analysis-ready format, and plug them into estimation pipelines. To do this
effectively, we need distribution functions — density, CDF, quantile, and
random generation — to extract summary statistics and informative values from
our delay distributions.

This tutorial shows how to:
1. Retrieve delay distributions from `EpiParameters.jl`
2. Convert them to EpiNow2-compatible formats
3. Use them in $R_t$ estimation for different diseases

!!! warning "Computation time"
    Each `epinow()` call takes several minutes for MCMC sampling.
"""

# ╔═╡ 2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6d02
md"""
## The pipeline: EpiParameters → EpiNow2

`EpiParameters.jl` returns `Distributions.jl` objects. EpiNow2.jl accepts
these directly via `discretise()`:

```
EpiParameters.jl          Distributions.jl          EpiNow2.jl
─────────────────    →    ────────────────    →    ───────────
epiparameter(...)         LogNormal(μ, σ)          discretise(d; max=14)
  .distribution           Gamma(α, θ)              gt_opts(pmf)
                          Weibull(k, λ)            delay_opts(pmf)
```

Unlike the R `epiparameter` package, which defines its own distribution class
and requires explicit conversion steps, in Julia we benefit from the standard
`Distributions.jl` ecosystem directly. This means any distribution that
`Distributions.jl` supports — and it supports a very wide range — can be used
with EpiNow2.jl after discretisation, with no special conversion needed.
"""

# ╔═╡ 3a4b5c6d-7e8f-9a0b-1c2d-3e4f5a6b7d02
md"""
## Example 1: COVID-19

### Finding the serial interval
"""

# ╔═╡ 4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8d02
begin
    covid_si = epiparameter(disease="COVID-19", epi_name="serial interval")
    covid_si_fitted = filter(p -> !isnothing(p.distribution), covid_si)
end

# ╔═╡ 5a6b7c8d-9e0f-1a2b-3c4d-5e6f7a8b9d02
md"""
We have $(length(covid_si_fitted)) fitted serial interval distributions for
COVID-19. Let's use the first one as a proxy for the generation time:
"""

# ╔═╡ a1b2c3d4-e5f6-7a8b-9c0d-a1b2c3d4e5f6
md"""
We may get multiple entries for the same disease and parameter type because
different studies have estimated the distribution. Entries can differ in
distribution family (LogNormal, Gamma, Weibull), parameter values, sample
size, and study population. When multiple entries exist, it is good practice
to choose the one with the largest sample size or that best matches your
study context (region, variant, etc.).
"""

# ╔═╡ 6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0d02
begin
    covid_gt_dist = covid_si_fitted[1].distribution
    covid_gt = discretise(covid_gt_dist; max=14)
end

# ╔═╡ 7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1d02
md"""
### Finding the incubation period
"""

# ╔═╡ 8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2d02
begin
    covid_ip = epiparameter(disease="COVID-19", epi_name="incubation period")
    covid_ip_dist = filter(p -> !isnothing(p.distribution), covid_ip)[1].distribution
    covid_ip_pmf = discretise(covid_ip_dist; max=14)
end

# ╔═╡ 9a0b1c2d-3e4f-5a6b-7c8d-9e0f1a2b3d02
md"""
### Visualising the delay distributions
"""

# ╔═╡ aa1b2c3d-4e5f-6a7b-8c9d-ae1f2a3b4d02
let
    fig = Figure(size=(700, 350))
    ax = Axis(fig[1, 1]; xlabel="Days", ylabel="Probability",
              title="COVID-19 delay distributions (discretised)")

    barplot!(ax, 0:length(covid_gt.pmf)-1, covid_gt.pmf;
             width=0.4, offset=-0.2, color=(:steelblue, 0.7),
             label="Serial interval (μ=$(round(mean(covid_gt_dist); digits=1))d)")
    barplot!(ax, 0:length(covid_ip_pmf.pmf)-1, covid_ip_pmf.pmf;
             width=0.4, offset=0.2, color=(:firebrick, 0.7),
             label="Incubation period (μ=$(round(mean(covid_ip_dist); digits=1))d)")
    axislegend(ax; position=:rt)
    fig
end

# ╔═╡ b1c2d3e4-f5a6-7b8c-9d0e-b1c2d3e4f5a6
md"""
### Why discretisation matters

EpiNow2 operates on discrete-time count data (daily cases and deaths), so
continuous distributions must be discretised to probability mass functions
(PMFs). The `discretise()` function converts a continuous distribution to a
discrete one, computing the probability mass in each daily interval. The
`max` parameter sets the maximum delay to consider — values beyond this are
truncated. Choosing an appropriate `max` is important: too small and you lose
probability mass from the tail; too large and you add unnecessary computation.
"""

# ╔═╡ ba2b3c4d-5e6f-7a8b-9c0d-be2f3a4b5d02
md"""
### Estimating Rt for COVID-19
"""

# ╔═╡ ca3b4c5d-6e7f-8a9b-0c1d-ce3f4a5b6d02
begin
    covid_cases = example_confirmed()

    covid_result = epinow(covid_cases;
        generation_time = gt_opts(covid_gt),
        delays = delay_opts(covid_ip_pmf),
        inference = inference_opts(samples=1000, warmup=200, chains=2,
                                   progress=false),
        verbose = false)
end

# ╔═╡ da4b5c6d-7e8f-9a0b-1c2d-de4f5a6b7d02
let
    rt = covid_result.estimates.rt
    fig = Figure(size=(700, 350))
    ax = Axis(fig[1, 1]; xlabel="Date", ylabel="Rt",
              title="COVID-19 Rt estimate (using literature delays)")
    band!(ax, 1:nrow(rt), rt.lower_90, rt.upper_90; color=(:steelblue, 0.15))
    band!(ax, 1:nrow(rt), rt.lower_50, rt.upper_50; color=(:steelblue, 0.3))
    lines!(ax, 1:nrow(rt), rt.median; linewidth=2, color=:steelblue)
    hlines!(ax, [1.0]; color=:red, linestyle=:dash, linewidth=1.5)
    ax.xticks = (1:10:nrow(rt), string.(rt.date[1:10:end]))
    ax.xticklabelrotation = π/4
    fig
end

# ╔═╡ ea5b6c7d-8e9f-0a1b-2c3d-ee5f6a7b8d02
md"""
## Example 2: Influenza

For influenza, the generation time is typically modelled as a **Weibull**
distribution. Since EpiNow2 works with discretised PMFs, any distribution
family works — we just discretise it.

This example illustrates a key advantage of the discretisation approach: you
do not need to worry about which distribution families your analysis tool
supports. Whether the literature reports a Gamma, LogNormal, Weibull, or any
other distribution, discretising to a PMF creates a universal format that any
downstream tool can use.
"""

# ╔═╡ fa6b7c8d-9e0f-1a2b-3c4d-fe6f7a8b9d02
begin
    flu_gt_results = epiparameter(disease="Influenza", epi_name="generation time")
    flu_gt_fitted = filter(p -> !isnothing(p.distribution), flu_gt_results)
end

# ╔═╡ 0b7b8c9d-0e1f-2a3b-4c5d-0f7f8a9b0d02
md"""
The influenza generation time is a **$(typeof(flu_gt_fitted[1].distribution).name.name)** distribution —
not natively supported by many Rt estimation tools. But since we discretise
to a PMF, any distribution family works:

```julia
flu_gt = discretise(flu_gt_fitted[1].distribution; max=10)
# → NonParametricDist that EpiNow2 accepts directly
```
"""

# ╔═╡ 1b8b9c0d-1e2f-3a4b-5c6d-1f8f9a0b1d02
let
    d = flu_gt_fitted[1].distribution
    x = range(0, quantile(d, 0.99); length=200)
    fig = Figure(size=(600, 300))
    ax = Axis(fig[1, 1]; xlabel="Days", ylabel="Density",
              title="Influenza generation time ($(typeof(d).name.name), μ=$(round(mean(d); digits=1))d)")
    lines!(ax, x, pdf.(d, x); linewidth=2.5, color=:seagreen)
    fig
end

# ╔═╡ 2b9b0c1d-2e3f-4a5b-6c7d-2f9f0a1b2d02
md"""
## Composing delays

Estimating $R_t$ requires accounting for all delays between infection and
observation. If cases are reported based on confirmation date, the total delay
includes both the incubation period (infection → symptoms) and the reporting
delay (symptoms → confirmation). Rather than manually computing this combined
distribution, we can compose the individual PMFs using `+`, which performs
discrete convolution.

When cases are reported with a delay (incubation + reporting), we compose
the PMFs with `+`:

```julia
total_delay = discretise(incubation_dist; max=14) +
              discretise(reporting_dist; max=7)
```

This convolves the two PMFs — the resulting distribution represents the
total delay from infection to case report.

As a practical application, the serial interval is also useful for optimising
contact tracing — it provides a time window for containment. Knowing the
probability that an infector's symptoms appeared within a certain number of
days helps determine how far back to trace contacts.
"""

# ╔═╡ 3ba01c2d-3e4f-5a6b-7c8d-3fa01a2b3d02
let
    composed = covid_ip_pmf + discretise(LogNormal(0.5, 0.5); max=7)
    fig = Figure(size=(600, 300))
    ax = Axis(fig[1, 1]; xlabel="Days", ylabel="Probability",
              title="Composed delay: incubation + reporting")
    barplot!(ax, 0:length(composed.pmf)-1, composed.pmf; color=(:purple, 0.6))
    fig
end

# ╔═╡ 4ba12c3d-4e5f-6a7b-8c9d-4fa12a3b4d02
md"""
## Practical considerations

### Choosing the right delay

| Purpose | Delay to use | EpiNow2 option |
|---|---|---|
| $R_t$ estimation | Generation time or serial interval | `gt_opts()` |
| Accounting for reporting lag | Incubation + reporting delay | `delay_opts()` |
| CFR estimation | Onset-to-death | `CFR.cfr_static(delay_density=...)` |

### When parameters are uncertain

In real outbreaks, we rarely know delay distributions exactly. Sources of
uncertainty include: small sample sizes in the original studies, variation
between populations, and changes over time (e.g., reporting delays may
shorten as testing capacity increases). `UncertainDistribution` handles this
by placing priors on the distribution parameters, so that this uncertainty is
propagated through to the final estimates rather than being ignored.

For uncertain delay parameters, use `UncertainDistribution` which places
priors on the distribution parameters and marginalises over them during
inference:

```julia
uncertain_gt = UncertainDistribution(
    (μ, σ) -> LogNormal(μ, σ),
    [Normal(1.6, 0.2), truncated(Normal(0.5, 0.1); lower=0.01)],
    14.0
)
```
"""

# ╔═╡ 5ba23c4d-5e6f-7a8b-9c0d-5fa23a4b5d02
md"""
## Key points

- `EpiParameters.jl` returns `Distributions.jl` objects that plug directly into
  EpiNow2.jl via `discretise()`
- **Any distribution family** works — discretisation to PMF is universal
- **Compose delays** with `+` (convolution of PMFs)
- Use `UncertainDistribution` when delay parameters are themselves uncertain
- Always check the **source** of your delay distribution (citation, sample size,
  region) via the `EpiParam` metadata
"""

# ╔═╡ Cell order:
# ╟─1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5d02
# ╠═0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4d02
# ╟─2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6d02
# ╟─3a4b5c6d-7e8f-9a0b-1c2d-3e4f5a6b7d02
# ╠═4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8d02
# ╟─5a6b7c8d-9e0f-1a2b-3c4d-5e6f7a8b9d02
# ╟─a1b2c3d4-e5f6-7a8b-9c0d-a1b2c3d4e5f6
# ╠═6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0d02
# ╟─7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1d02
# ╠═8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2d02
# ╟─9a0b1c2d-3e4f-5a6b-7c8d-9e0f1a2b3d02
# ╠═aa1b2c3d-4e5f-6a7b-8c9d-ae1f2a3b4d02
# ╟─b1c2d3e4-f5a6-7b8c-9d0e-b1c2d3e4f5a6
# ╟─ba2b3c4d-5e6f-7a8b-9c0d-be2f3a4b5d02
# ╠═ca3b4c5d-6e7f-8a9b-0c1d-ce3f4a5b6d02
# ╠═da4b5c6d-7e8f-9a0b-1c2d-de4f5a6b7d02
# ╟─ea5b6c7d-8e9f-0a1b-2c3d-ee5f6a7b8d02
# ╠═fa6b7c8d-9e0f-1a2b-3c4d-fe6f7a8b9d02
# ╟─0b7b8c9d-0e1f-2a3b-4c5d-0f7f8a9b0d02
# ╠═1b8b9c0d-1e2f-3a4b-5c6d-1f8f9a0b1d02
# ╟─2b9b0c1d-2e3f-4a5b-6c7d-2f9f0a1b2d02
# ╠═3ba01c2d-3e4f-5a6b-7c8d-3fa01a2b3d02
# ╟─4ba12c3d-4e5f-6a7b-8c9d-4fa12a3b4d02
# ╟─5ba23c4d-5e6f-7a8b-9c0d-5fa23a4b5d02
