### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4c5e
begin
    import Pkg
    Pkg.develop(path=joinpath(@__DIR__, "..", "CFR"))
    Pkg.develop(path=joinpath(@__DIR__, "..", "EpiParameters"))
    using CFR
    using EpiParameters
    using Distributions
    using DataFrames
    using Dates
    using CSV
    using CairoMakie
end

# ╔═╡ 1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c6e
md"""
# Estimation of Outbreak Severity

During an outbreak, two properties of a pathogen jointly determine its pandemic
potential: **transmissibility** (how easily it spreads) and **clinical severity**
(how much harm it causes). Frameworks such as the CDC's Pandemic Severity
Assessment Framework explicitly use both axes to guide the public health
response. This tutorial focuses on severity.

The **case fatality ratio (CFR)** is the conditional probability of death given a
confirmed diagnosis. It is one of the most widely cited severity measures, yet
estimating it reliably during an ongoing epidemic is harder than it first
appears — simply dividing deaths by cases introduces systematic bias.

Severity estimation also helps determine whether an outbreak differs from
historical patterns. An unexpectedly high CFR may signal a new strain, reduced
population immunity, or overwhelmed health systems; an unexpectedly low CFR
could indicate improved treatment or a milder variant.

This tutorial covers:

1. The **naive CFR** and why it is biased
2. **Delay-adjusted CFR** using the method of Nishiura et al. (2009)
3. **Time-varying CFR** for monitoring severity over the course of an outbreak
"""

# ╔═╡ a1b2c3d4-e5f6-7a8b-9c0d-a1b2c3d4e5f6
md"""
### What does "severity" capture?

Different data sources capture different levels of disease severity. Severe and
critical cases are typically identified through hospital surveillance. Mildly
symptomatic cases may be picked up through routine case reporting, while
mild, subclinical, and asymptomatic infections are only detected through active
contact tracing or serological surveys.

Because the CFR is calculated among **confirmed cases** only, it reflects
severity conditional on being detected by the surveillance system. It therefore
exceeds the **infection fatality risk (IFR)**, which includes all infections
regardless of whether they were diagnosed. Interpreting CFR estimates requires
understanding the case ascertainment context.
"""

# ╔═╡ 2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6c7e
md"""
## The data: Ebola 1976

We'll use data from the first known Ebola outbreak in Yambuku, Zaire (now DRC)
in 1976.
"""

# ╔═╡ 3a4b5c6d-7e8f-9a0b-1c2d-3e4f5a6b7c8e
ebola = CSV.read(joinpath(@__DIR__, "ebola1976.csv"), DataFrame)

# ╔═╡ 4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8c9e
let
    fig = Figure(size=(700, 350))
    ax = Axis(fig[1, 1];
              xlabel="Date", ylabel="Count",
              title="Ebola 1976 outbreak — daily cases and deaths")
    barplot!(ax, 1:nrow(ebola), ebola.cases; color=(:steelblue, 0.7), label="Cases")
    barplot!(ax, 1:nrow(ebola), ebola.deaths; color=(:firebrick, 0.7), label="Deaths")
    ax.xticks = (1:14:nrow(ebola), string.(ebola.date[1:14:end]))
    ax.xticklabelrotation = π/4
    axislegend(ax; position=:rt)
    fig
end

# ╔═╡ 5a6b7c8d-9e0f-1a2b-3c4d-5e6f7a8b9c0e
md"""
Total cases: **$(sum(ebola.cases))**, total deaths: **$(sum(ebola.deaths))**
"""

# ╔═╡ 6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0c1e
md"""
## Naive CFR

The simplest estimate of the CFR at time ``t`` is the **naive CFR**:

```math
b_t = \frac{D_t}{C_t}
```

where ``D_t`` is the cumulative number of deaths and ``C_t`` is the cumulative
number of confirmed cases up to time ``t``.

This estimator is **biased downward** during an ongoing epidemic. The reason is
straightforward: recently confirmed cases have not yet had time for their
outcome — death or recovery — to be observed. The denominator includes these
unresolved cases, but the numerator does not yet include any deaths among them.
The faster the epidemic is growing and the longer the delay from onset to death,
the greater this bias.
"""

# ╔═╡ 7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1c2e
naive_cfr = cfr_static(ebola)

# ╔═╡ 8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2c3e
md"""
Naive CFR: **$(round(naive_cfr.estimate * 100; digits=1))%**
(95% CI: $(round(naive_cfr.lower * 100; digits=1))–$(round(naive_cfr.upper * 100; digits=1))%)

During an ongoing epidemic, this **underestimates** the true CFR because cases
confirmed recently haven't had time to die yet.
"""

# ╔═╡ 9a0b1c2d-3e4f-5a6b-7c8d-9e0f1a2b3c4e
md"""
## Delay-adjusted CFR

The method of **Nishiura et al. (2009)** corrects for the right-censoring bias
by accounting for the fact that recently confirmed cases have not yet had enough
time for their outcome to be known. It uses the **onset-to-death delay
distribution** to estimate the proportion of cases at each time point that have
"had enough time" for death to have occurred, had it been going to occur.

In a real-time outbreak, the onset-to-death distribution may not be directly
estimable from the current data, especially early on. In practice, analysts
often rely on delay distributions from previous outbreaks of the same pathogen
or from literature databases such as `EpiParameters.jl`.
"""

# ╔═╡ aa1b2c3d-4e5f-6a7b-8c9d-ae1f2a3b4c5e
begin
    ebola_otd = epiparameter(disease="Ebola", epi_name="onset to death")
    delay_param = filter(p -> !isnothing(p.distribution), ebola_otd)[1]
    delay = delay_param.distribution
end

# ╔═╡ ba2b3c4d-5e6f-7a8b-9c0d-be2f3a4b5c6e
md"""
Onset-to-death: **$(typeof(delay).name.name)** with mean **$(round(mean(delay); digits=1)) days**
"""

# ╔═╡ ca3b4c5d-6e7f-8a9b-0c1d-ce3f4a5b6c7e
adjusted_cfr = cfr_static(ebola; delay_density=delay)

# ╔═╡ da4b5c6d-7e8f-9a0b-1c2d-de4f5a6b7c8e
md"""
Delay-adjusted CFR: **$(round(adjusted_cfr.estimate * 100; digits=1))%**
(95% CI: $(round(adjusted_cfr.lower * 100; digits=1))–$(round(adjusted_cfr.upper * 100; digits=1))%)
"""

# ╔═╡ ea5b6c7d-8e9f-0a1b-2c3d-ee5f6a7b8c9e
md"""
## Time-varying CFR

Rather than producing a single summary estimate, we can track how the CFR
estimate evolves over time using `cfr_time_varying`. This function computes what
the estimated CFR would have been on each day of the outbreak, using only the
data available up to that point.

The `burn_in` period (here 14 days) excludes the earliest days when both case
and death counts are very small, leading to highly unstable estimates. After the
burn-in, the time-varying estimate should gradually converge, with the final
value matching the corresponding static estimate.
"""

# ╔═╡ fa6b7c8d-9e0f-1a2b-3c4d-fe6f7a8b9c0e
begin
    tv_naive = cfr_time_varying(ebola; burn_in=14)
    tv_adjusted = cfr_time_varying(ebola; delay_density=delay, burn_in=14)
end

# ╔═╡ 0b7b8c9d-0e1f-2a3b-4c5d-0f7f8a9b0c1e
let
    fig = Figure(size=(700, 400))
    ax = Axis(fig[1, 1];
              xlabel="Day of outbreak", ylabel="CFR",
              title="Time-varying CFR — naive vs delay-adjusted")

    # Day index (burn_in already applied)
    days_n = 15:nrow(ebola)
    days_a = 15:nrow(ebola)

    # Naive
    band!(ax, days_n, tv_naive.lower, tv_naive.upper; color=(:steelblue, 0.2))
    lines!(ax, days_n, tv_naive.cfr; linewidth=2, color=:steelblue, label="Naive")

    # Adjusted
    valid = .!isnan.(tv_adjusted.cfr)
    band!(ax, days_a[valid], tv_adjusted.lower[valid], tv_adjusted.upper[valid];
          color=(:firebrick, 0.2))
    lines!(ax, days_a[valid], tv_adjusted.cfr[valid];
           linewidth=2, color=:firebrick, label="Delay-adjusted")

    axislegend(ax; position=:rb)
    fig
end

# ╔═╡ 1b8b9c0d-1e2f-3a4b-5c6d-1f8f9a0b1c2e
md"""
The delay-adjusted estimate converges faster to the true CFR, especially early
in the outbreak when most cases have not yet had time for their outcome to be
known.
"""

# ╔═╡ b2c3d4e5-f6a7-8b9c-0d1e-b2c3d4e5f6a7
md"""
### Why does this matter in practice?

The difference between naive and delay-adjusted estimates is not merely
academic. During the 2009 H1N1 pandemic, early CFR estimates varied dramatically
between countries, partly because of this right-censoring bias interacting with
different epidemic growth rates and reporting delays. Having a method that
corrects for the delay between onset and death enables earlier, more reliable
severity assessment — which in turn supports timelier decisions about
interventions, resource allocation, and public communication.

As the outbreak progresses and more cases reach their final outcome, the naive
and adjusted estimates converge. The value of delay adjustment is greatest
precisely when it is most needed: early in the outbreak, when decisions are most
uncertain and most consequential.
"""

# ╔═╡ 2b9b0c1d-2e3f-4a5b-6c7d-2f9f0a1b2c3e
md"""
## Sensitivity to the delay distribution

The choice of delay distribution is one of the most important modelling
decisions when estimating a delay-adjusted CFR. What happens if the assumed
onset-to-death delay is wrong?

If the assumed delay is **too short**, cases are treated as having resolved
sooner than they actually have, and the adjusted CFR will be closer to the naive
estimate. If the assumed delay is **too long**, the method overestimates the
degree of right-censoring and may overestimate severity.

This is why using well-characterised delay distributions from the literature —
for instance via `EpiParameters.jl` — is important, and why sensitivity analyses
exploring a range of plausible delay distributions should be standard practice.
"""

# ╔═╡ 3ba01c2d-3e4f-5a6b-7c8d-3fa01a2b3c4e
let
    delays = [
        ("Short (μ=4d)", Gamma(2.0, 2.0)),
        ("Estimated (μ=$(round(mean(delay); digits=0))d)", delay),
        ("Long (μ=16d)", Gamma(2.0, 8.0)),
    ]

    fig = Figure(size=(600, 350))
    ax = Axis(fig[1, 1];
              xlabel="Assumed delay distribution", ylabel="CFR estimate",
              title="Sensitivity of delay-adjusted CFR to delay assumption",
              xticks=(1:3, [d[1] for d in delays]))

    for (i, (label, d)) in enumerate(delays)
        est = cfr_static(ebola; delay_density=d)
        scatter!(ax, [i], [est.estimate]; markersize=12, color=:steelblue)
        rangebars!(ax, [i], [est.lower], [est.upper]; color=:steelblue, linewidth=2)
    end

    hlines!(ax, [naive_cfr.estimate]; color=:grey, linestyle=:dash,
            label="Naive CFR")
    axislegend(ax; position=:rb)
    fig
end

# ╔═╡ 4ba12c3d-4e5f-6a7b-8c9d-4fa12a3b4c5e
md"""
## Key points

- The **naive CFR** (deaths/cases) underestimates true severity during an ongoing
  outbreak because of the delay between onset and death
- **Delay-adjusted CFR** (Nishiura et al. 2009) corrects for right-censoring
- `EpiParameters.jl` provides delay distributions that plug directly into `CFR.jl`
- **Time-varying CFR** tracks how severity estimates evolve over time
- The estimate is **sensitive** to the assumed delay distribution — always report
  this choice and consider sensitivity analyses
- The CFR only captures severity among **confirmed cases**. The **infection
  fatality risk (IFR)** — which includes all infections, including those never
  diagnosed — requires additional information about case ascertainment, typically
  from serological surveys or capture-recapture methods
"""

# ╔═╡ Cell order:
# ╟─1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c6e
# ╠═0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4c5e
# ╟─a1b2c3d4-e5f6-7a8b-9c0d-a1b2c3d4e5f6
# ╟─2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6c7e
# ╠═3a4b5c6d-7e8f-9a0b-1c2d-3e4f5a6b7c8e
# ╠═4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8c9e
# ╟─5a6b7c8d-9e0f-1a2b-3c4d-5e6f7a8b9c0e
# ╟─6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0c1e
# ╠═7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1c2e
# ╟─8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2c3e
# ╟─9a0b1c2d-3e4f-5a6b-7c8d-9e0f1a2b3c4e
# ╠═aa1b2c3d-4e5f-6a7b-8c9d-ae1f2a3b4c5e
# ╟─ba2b3c4d-5e6f-7a8b-9c0d-be2f3a4b5c6e
# ╠═ca3b4c5d-6e7f-8a9b-0c1d-ce3f4a5b6c7e
# ╟─da4b5c6d-7e8f-9a0b-1c2d-de4f5a6b7c8e
# ╟─ea5b6c7d-8e9f-0a1b-2c3d-ee5f6a7b8c9e
# ╠═fa6b7c8d-9e0f-1a2b-3c4d-fe6f7a8b9c0e
# ╠═0b7b8c9d-0e1f-2a3b-4c5d-0f7f8a9b0c1e
# ╟─1b8b9c0d-1e2f-3a4b-5c6d-1f8f9a0b1c2e
# ╟─b2c3d4e5-f6a7-8b9c-0d1e-b2c3d4e5f6a7
# ╟─2b9b0c1d-2e3f-4a5b-6c7d-2f9f0a1b2c3e
# ╠═3ba01c2d-3e4f-5a6b-7c8d-3fa01a2b3c4e
# ╟─4ba12c3d-4e5f-6a7b-8c9d-4fa12a3b4c5e
