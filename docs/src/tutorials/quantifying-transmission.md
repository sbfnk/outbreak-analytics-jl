# Quantifying Transmission

A key question during an outbreak is: **is the epidemic growing or shrinking?**
The **time-varying reproduction number** ``R_t`` answers this — it estimates how
many secondary infections each case generates at time ``t``.

- ``R_t > 1``: epidemic is growing
- ``R_t < 1``: epidemic is shrinking
- ``R_t = 1``: epidemic is stable

The transmission intensity of an outbreak is quantified using two key metrics: the **reproduction number**, which informs on the *strength* of transmission (how many new cases are expected from each existing case); and the **growth rate**, which informs on the *speed* of transmission (how rapidly the outbreak is spreading or declining).

The basic reproduction number, R₀, is the average number of secondary cases caused by one infectious individual in an entirely susceptible population. But in an ongoing outbreak, the population does not remain entirely susceptible — those who recover are typically immune, and behaviour changes or interventions may reduce transmission. The *effective* reproduction number ``R_t`` accounts for these factors at time ``t``.

To estimate ``R_t`` from case data, we must account for delays between the date of infection and the date cases are reported. Observations at time ``t`` reflect transmission events from some time in the past.

### Bayesian inference

EpiNow2.jl uses a Bayesian inference framework: it combines prior knowledge (prior distributions) with observed data (likelihood) to estimate posterior distributions for ``R_t`` and other quantities. This naturally quantifies uncertainty in our estimates.

### Reproduction number and growth rate

The reproduction number tells us the *strength* of transmission (how many secondary cases per primary case), whilst the growth rate tells us the *speed* (how fast case numbers are doubling or halving). These two metrics are related through the generation time distribution — for a given ``R_t``, a shorter generation time produces a faster growth rate.

The doubling time (when growing) or halving time (when declining) can be calculated from the growth rate ``r`` as ``\ln(2) / |r|``.

!!! warning "Computation time"
    MCMC sampling takes several minutes. The examples below use lighter
    inference settings for speed.

## The data

In an outbreak situation, data are usually available on reported dates only. We must use estimation methods that account for delays when trying to understand how transmission changes over time.

We use synthetic case data built into EpiNow2.jl:

```@example rt
using EpiNow2
using EpiParameters
using Distributions
using DataFrames
using Dates
using CairoMakie

cases = example_confirmed()
first(cases, 5)
```

```@example rt
fig = Figure(size=(700, 300))
ax = Axis(fig[1, 1]; xlabel="Date", ylabel="Confirmed cases",
          title="Daily confirmed cases")
barplot!(ax, 1:nrow(cases), cases.confirm; color=(:steelblue, 0.7))
ax.xticks = (1:10:nrow(cases), string.(cases.date[1:10:end]))
ax.xticklabelrotation = π/4
fig
```

## Specifying delay distributions

To estimate ``R_t``, EpiNow2 needs to know about the delays between infection and observation. The delay from infection to case report typically involves multiple processes:

- **Infection → symptom onset** (incubation period)
- **Symptom onset → case notification** (reporting delay)

Each of these has individual-level variation, which we capture with probability distributions. We must use distributions for positive values only (Gamma, LogNormal) since delays cannot be negative.

| Data source | Delay(s) to specify |
|---|---|
| Symptom onset dates | Incubation period |
| Case report dates | Incubation period + reporting delay |
| Hospitalisation dates | Incubation period + onset-to-hospitalisation |

### Why use positive-valued distributions?

Delays between events (incubation period, reporting delay, generation time) must be non-negative — we cannot have a negative number of days between infection and symptom onset. Therefore, we specify distributions that only take positive values. `Gamma()` and `LogNormal()` are the most commonly used families for this purpose. The exception is the serial interval, which *can* be negative when there is pre-symptomatic transmission, but this is usually handled by truncation at zero.

EpiNow2 also requires a **generation time** distribution — the delay between infection of a primary case and infection of a secondary case:

```@example rt
generation_time = discretise(LogNormal(1.2, 0.5); max=14)
incubation = discretise(LogNormal(1.4, 0.4); max=14)
reporting = discretise(LogNormal(0.5, 0.5); max=7)
total_delay = incubation + reporting
nothing # hide
```

!!! details "Coming from R?"
    In R's EpiNow2, you would write `EpiNow2::LogNormal(mean = 4, sd = 2, max = 20)` or `EpiNow2::Gamma(mean = 4, sd = 2, max = 20)`. In Julia, you create a standard `Distributions.jl` object (`LogNormal(1.2, 0.5)`) and then discretise it as a separate step. This is more composable — the same distribution object can be used with any package, not just EpiNow2.

Delays compose with `+` (convolution of PMFs).

Composing delays with `+` performs discrete convolution of the PMFs — the resulting distribution represents the total delay from infection to case report.

### Fixed vs uncertain delay distributions

In the example above, we specified delays with fixed parameter values. In practice, there is often uncertainty about the true values. EpiNow2.jl supports uncertain distributions where the parameters themselves have prior distributions — this is covered in the [Using Delay Distributions in Analysis](@ref) tutorial with `UncertainDistribution`.

For this tutorial, we use fixed distributions to keep computation time manageable. Using uncertain distributions is generally recommended for production analyses as it accounts for additional uncertainty.

## Running the estimation

```@example rt
result = epinow(cases;
    generation_time = gt_opts(generation_time),
    delays = delay_opts(total_delay),
    rt = rt_opts(prior=LogNormal(log(1.0), 1.0)),
    inference = inference_opts(samples=1000, warmup=250, chains=3,
                               progress=false),
    verbose = false)
nothing # hide
```

EpiNow2 performs Bayesian inference using MCMC methods via Turing.jl. The model uses a renewal equation to link infections over time through the generation time distribution, and a Gaussian process to smooth ``R_t``. Running multiple chains helps assess convergence.

## Reproduction number over time

```@example rt
rt = result.estimates.rt
fig = Figure(size=(700, 350))
ax = Axis(fig[1, 1]; xlabel="Date", ylabel="Rt",
          title="Time-varying reproduction number")

band!(ax, 1:nrow(rt), rt.lower_90, rt.upper_90; color=(:steelblue, 0.15))
band!(ax, 1:nrow(rt), rt.lower_50, rt.upper_50; color=(:steelblue, 0.3))
lines!(ax, 1:nrow(rt), rt.median; linewidth=2, color=:steelblue)
hlines!(ax, [1.0]; color=:red, linestyle=:dash, linewidth=1.5)

ax.xticks = (1:10:nrow(rt), string.(rt.date[1:10:end]))
ax.xticklabelrotation = π/4
fig
```

In the plot, the shaded bands represent credible intervals — the darker 50% band contains the most likely values, whilst the lighter 90% band captures most of the uncertainty. The uncertainty increases towards the end of the data because the most recent infections are least likely to have been observed yet.

Estimates are categorised as: **Estimate** (uses all available data) and **Estimate based on partial data** (uses less data for the most recent dates, hence wider intervals).

### Summary metrics

From the estimation, we can extract key summary metrics at the latest date:

```@example rt
latest_rt = result.estimates.rt[end, :]
println("Latest Rt estimate: $(round(latest_rt.median; digits=2)) " *
        "(90% CrI: $(round(latest_rt.lower_90; digits=2))–$(round(latest_rt.upper_90; digits=2)))")
```

The **expected change in reports** can be assessed from the posterior probability that ``R_t < 1``:

| Probability (``R_t < 1``) | Expected change |
|---|---|
| p < 0.05 | Increasing |
| 0.05 ≤ p < 0.4 | Likely increasing |
| 0.4 ≤ p < 0.6 | Stable |
| 0.6 ≤ p < 0.95 | Likely decreasing |
| 0.95 ≤ p ≤ 1 | Decreasing |

### Estimated infections

We can also visualise the estimated infection curve alongside the observed cases:

```@example rt
inf = result.estimates.infections
fig = Figure(size=(700, 350))
ax = Axis(fig[1, 1]; xlabel="Date", ylabel="Count",
          title="Estimated infections vs observed cases")

band!(ax, 1:nrow(inf), inf.lower_90, inf.upper_90; color=(:steelblue, 0.15))
lines!(ax, 1:nrow(inf), inf.median; linewidth=2, color=:steelblue,
       label="Infections (estimated)")

obs_idx = 1:min(nrow(cases), nrow(inf))
scatter!(ax, obs_idx, cases.confirm[obs_idx]; markersize=5,
         color=:black, label="Observed cases")

ax.xticks = (1:10:nrow(inf), string.(inf.date[1:10:end]))
ax.xticklabelrotation = π/4
axislegend(ax; position=:lt)
fig
```

The estimated infections precede the observed cases in time, reflecting the delays between infection and reporting.

### How does EpiNow2 work?

EpiNow2 uses the **renewal equation** to model the relationship between infections over time:

```math
I_t = R_t \sum_{s=1}^{T} I_{t-s} \, w_s
```

where ``I_t`` is the number of infections at time ``t``, ``R_t`` is the effective reproduction number, and ``w_s`` is the probability mass function of the generation time at lag ``s``. The current infections are the product of ``R_t`` and the sum of past infections weighted by the generation time distribution.

``R_t`` is modelled as a smooth function over time using a Gaussian process prior, which allows it to vary but penalises rapid fluctuations. Delay distributions are then applied to map from infections to observations (accounting for incubation, reporting, etc.).

## Using literature delays

Using literature-sourced delay distributions is particularly valuable early in an outbreak, when local data on delays may not yet be available.

Instead of specifying delays manually, pull them from `EpiParameters.jl`:

```@example rt
using EpiParameters

# COVID-19 serial interval is Normal (can be negative), so truncate at 0
gt_results = epiparameter(disease="COVID-19", epi_name="serial interval")
gt_dist = filter(p -> !isnothing(p.distribution), gt_results)[1].distribution
generation_time_lit = discretise(truncated(gt_dist; lower=0); max=14)
```

This pipeline is covered in detail in the
[Using Delay Distributions in Analysis](@ref) tutorial.

## Key points

- **``R_t``** measures whether an epidemic is growing or shrinking
- **EpiNow2.jl** estimates ``R_t`` via a Bayesian renewal equation model
- You must specify a **generation time**; reporting delays are optional but improve estimates
- Delays can be **composed** with `+` (convolution)
- Output includes **credible intervals** — always report uncertainty
- The **growth rate** complements ``R_t`` by quantifying the speed of transmission
- EpiNow2 uses a **renewal equation** with Gaussian process smoothing of ``R_t``
- Summary metrics (latest ``R_t``, expected change) help communicate findings to decision makers


---

*Adapted from the [Epiverse-TRACE tutorials](https://epiverse-trace.github.io/tutorials/), © Epiverse-TRACE contributors, licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).*
