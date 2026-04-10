# Using Delay Distributions in Analysis

Epidemiological delay distributions — incubation period, generation time,
serial interval, reporting delay — are essential inputs for transmission
estimation, forecasting, and severity assessment.

This tutorial shows how to:
1. Retrieve delay distributions from `EpiParameters.jl`
2. Convert them to EpiNow2-compatible formats
3. Use them in ``R_t`` estimation for different diseases

This tutorial integrates the previous two: using the parameter database from the [Access Epidemiological Delay Distributions](@ref) tutorial with the ``R_t`` estimation framework from the [Quantifying Transmission](@ref) tutorial. The key workflow is: retrieve parameters from the literature, convert them to analysis-ready format, and plug them into estimation pipelines.

Unlike the R `epiparameter` package which has its own distribution class, in Julia we benefit from the standard `Distributions.jl` ecosystem directly. This means any distribution that `Distributions.jl` supports can be used with EpiNow2.jl after discretisation — no special conversion needed.

!!! warning "Computation time"
    Each `epinow()` call takes several minutes for MCMC sampling.

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

The `discretise()` function converts a continuous distribution to a discrete probability mass function (PMF), computing the probability mass in each daily interval. The `max` parameter sets the maximum delay to consider — values beyond this are truncated. EpiNow2 operates on discrete-time count data (daily cases and deaths), which is why this discretisation step is necessary.

## Distribution functions

In Julia, the `Distributions.jl` package provides a consistent interface for working with any probability distribution. The four key functions are:

- `pdf(d, x)` — probability density at point `x`
- `cdf(d, x)` — cumulative probability up to `x`
- `quantile(d, p)` — the value below which probability `p` of the distribution lies
- `rand(d, n)` — generate `n` random samples

These functions work with *any* distribution stored in an `EpiParam` object's `.distribution` field. Let's demonstrate with the COVID-19 serial interval:

```@example delays
using EpiNow2
using EpiParameters
using Distributions
using DataFrames
using Dates
using CairoMakie

covid_si_example = filter(p -> !isnothing(p.distribution), 
    epiparameter(disease="COVID-19", epi_name="serial interval"))[1]
d_si = covid_si_example.distribution
```

```@example delays
# What is the probability density at day 5?
println("PDF at day 5: $(round(pdf(d_si, 5.0); digits=4))")

# What proportion of serial intervals are ≤ 7 days?
println("CDF at day 7: $(round(cdf(d_si, 7.0); digits=3))")

# Below what value do 95% of serial intervals fall?
println("95th percentile: $(round(quantile(d_si, 0.95); digits=1)) days")

# Generate 5 random serial interval values
println("Random samples: $(round.(rand(d_si, 5); digits=1))")
```

!!! details "Coming from R?"
    R's `epiparameter` package has its own `density()`, `cdf()`, `quantile()`, and `generate()` functions for `<epiparameter>` objects. In Julia, the standard `Distributions.jl` functions (`pdf`, `cdf`, `quantile`, `rand`) work directly on any distribution — including those retrieved from `EpiParameters.jl`. No special wrapper functions needed.

### Contact tracing window and the serial interval

The serial interval is important for optimising contact tracing, since it provides a time window for the containment of disease spread. We can use the CDF to quantify how expanding the tracing window improves detection of potential infectors.

```@example delays
# How much of the backward serial interval is captured by tracing 2 vs 6 days back?
p_2days = cdf(d_si, 2.0)
p_6days = cdf(d_si, 6.0)
println("Tracing 2 days back captures $(round(p_2days * 100; digits=1))% of infectors")
println("Tracing 6 days back captures $(round(p_6days * 100; digits=1))% of infectors")
println("Extending from 2→6 days captures $(round((p_6days - p_2days) * 100; digits=1))% more")
```

### Quarantine length and the incubation period

The incubation period distribution helps assess the appropriate length of quarantine or active monitoring. We can ask: after how many days following exposure will a given percentage of people who will develop symptoms actually show them?

```@example delays
covid_ip_example = filter(p -> !isnothing(p.distribution),
    epiparameter(disease="COVID-19", epi_name="incubation period"))[1]
d_ip = covid_ip_example.distribution
q99 = quantile(d_ip, 0.99)
println("99% of symptomatic cases show symptoms within $(round(q99; digits=0)) days of infection")
```

Is this consistent with quarantine recommendations during the COVID-19 pandemic?

## Finding COVID-19 delays

```@example delays
covid_ip = epiparameter(disease="COVID-19", epi_name="incubation period")
covid_ip_fitted = filter(p -> !isnothing(p.distribution), covid_ip)
```

The database may contain multiple entries for the same disease and parameter type from different studies. These entries may differ in distribution family (LogNormal, Gamma, Weibull), parameter values, sample size, and study population. When multiple entries exist, it is good practice to choose the one with the largest sample size or that best matches your study context (region, variant, etc.).

The incubation period is a LogNormal — well suited for discretisation:

```@example delays
covid_ip_dist = covid_ip_fitted[1].distribution
covid_ip_pmf = discretise(covid_ip_dist; max=14)
```

The database also has a COVID-19 serial interval, but it is a Normal distribution
which can take negative values. We truncate it at zero before discretising:

```@example delays
covid_si = epiparameter(disease="COVID-19", epi_name="serial interval")
covid_si_dist = filter(p -> !isnothing(p.distribution), covid_si)[1].distribution
covid_gt_dist = truncated(covid_si_dist; lower=0)
covid_gt = discretise(covid_gt_dist; max=14)
```

!!! details "Coming from R?"
    In R, `epiparameter::discretise()` converts an `<epiparameter>` object to a discrete distribution. In Julia, `discretise()` from EpiNow2.jl works on any `Distributions.jl` distribution directly — you don't need to wrap it in a special class first.

### Visualising the discretised delays

```@example delays
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
```

### Continuous vs discrete

When we discretise a continuous distribution, we convert the smooth probability density function (PDF) into a probability mass function (PMF) — a set of probabilities for each discrete day. This is the appropriate representation for daily count data:

```@example delays
fig = Figure(size=(700, 300))

ax1 = Axis(fig[1, 1]; xlabel="Days", ylabel="Density / Probability",
           title="Continuous PDF")
x_cont = range(0, 14; length=200)
lines!(ax1, x_cont, pdf.(covid_ip_dist, x_cont); linewidth=2, color=:steelblue)

ax2 = Axis(fig[1, 2]; xlabel="Days", ylabel="Probability",
           title="Discretised PMF")
barplot!(ax2, 0:length(covid_ip_pmf.pmf)-1, covid_ip_pmf.pmf; color=(:steelblue, 0.7))

fig
```

The PMF values sum to 1 (within the truncation range), and each bar represents the probability that the delay falls within that specific day.

### Choosing the maximum delay

The `max` parameter in `discretise()` sets the upper truncation point. A good default is the 99th percentile of the continuous distribution:

```@example delays
max_val = ceil(Int, quantile(covid_ip_dist, 0.99))
println("99th percentile of COVID-19 incubation period: $max_val days")
```

Setting `max` too low truncates meaningful probability mass; setting it too high wastes computation on negligible probabilities.

### Estimating Rt with literature delays

```@example delays
cases = example_confirmed()

covid_result = epinow(cases;
    generation_time = gt_opts(covid_gt),
    delays = delay_opts(covid_ip_pmf),
    inference = inference_opts(samples=1000, warmup=250, chains=3,
                               progress=false),
    verbose = false)
nothing # hide
```

```@example delays
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
```

## Any distribution family works

This illustrates a key advantage of the discretisation approach: you do not need to worry about which distribution families your downstream analysis tools support. Whether the literature reports a Gamma, LogNormal, Weibull, or any other distribution, discretising to a PMF creates a universal format.

!!! details "Coming from R?"
    The reason `pdf(d, x)`, `cdf(d, x)`, and `quantile(d, p)` work identically whether `d` is a `LogNormal`, `Gamma`, `Weibull`, or any other distribution is **multiple dispatch** — Julia's core design principle. The same function name dispatches to different implementations based on the *type* of its arguments. This is similar to R's S4 methods or S3 generic functions (e.g. `print.lm`, `print.glm`), but in Julia it is pervasive and central to the language. It's why the `EpiParameters.jl` → `Distributions.jl` → `EpiNow2.jl` pipeline works so seamlessly — every package speaks the same type language.

For influenza, the generation time is a **Weibull** — not natively supported
by many Rt tools. Since we discretise to a PMF, it works seamlessly:

```@example delays
flu_gt = epiparameter(disease="Influenza", epi_name="generation time")
flu_gt_fitted = filter(p -> !isnothing(p.distribution), flu_gt)
d = flu_gt_fitted[1].distribution
typeof(d)
```

```@example delays
x = range(0, quantile(d, 0.99); length=200)
fig = Figure(size=(600, 300))
ax = Axis(fig[1, 1]; xlabel="Days", ylabel="Density",
          title="Influenza generation time ($(typeof(d).name.name), μ=$(round(mean(d); digits=1))d)")
lines!(ax, x, pdf.(d, x); linewidth=2.5, color=:seagreen)
fig
```

## Composing delays

Estimating ``R_t`` requires accounting for all delays between infection and observation. If cases are reported based on confirmation date, the total delay includes both the incubation period (infection → symptoms) and the reporting delay (symptoms → confirmation). Rather than manually computing this combined distribution, we compose the individual PMFs using `+`, which performs discrete convolution.

The serial interval is also useful for optimising contact tracing, since it provides a time window for containment. Knowing the probability that an infector's symptoms appeared within a certain number of days helps determine how far back to trace contacts.

When cases are reported with a delay, compose the PMFs with `+`:

```@example delays
composed = covid_ip_pmf + discretise(LogNormal(0.5, 0.5); max=7)
fig = Figure(size=(600, 300))
ax = Axis(fig[1, 1]; xlabel="Days", ylabel="Probability",
          title="Composed delay: incubation + reporting")
barplot!(ax, 0:length(composed.pmf)-1, composed.pmf; color=(:purple, 0.6))
fig
```

## Uncertain delay parameters

In real outbreaks, we rarely know delay distributions with certainty. Sources of uncertainty include: small sample sizes in the original studies, variation between populations, and changes over time (e.g. reporting delays may shorten as testing capacity increases).

When delay parameters are themselves uncertain, use `UncertainDistribution` to place priors on them:

```@example delays
uncertain_gt = UncertainDistribution(
    (μ, σ) -> LogNormal(μ, σ),
    [Normal(1.6, 0.2), truncated(Normal(0.5, 0.1); lower=0.01)],
    14.0
)
```

This marginalises over parameter uncertainty during inference.

## Challenge: Ebola Rt with reporting delays

Estimate the effective reproduction number for Ebola using delay distributions from the literature:

1. Find the Ebola serial interval or generation time in `EpiParameters.jl`
2. Find the Ebola incubation period for use as a reporting delay
3. Run `epinow()` with both distributions
4. Compare the result to an estimate using only the generation time (no reporting delay) — does adding the delay change the Rt estimate? Does it change the uncertainty?

!!! hint
    Look at the distribution family for the Ebola serial interval — it may not be LogNormal. Since `discretise()` works with any `Distributions.jl` distribution, any family will work.

## Challenge: Influenza with a Weibull distribution

Some delay distributions in the literature use less common families like the Weibull. Try estimating Rt for influenza:

1. Retrieve the influenza generation time from `EpiParameters.jl`
2. Note the distribution family — it should be a Weibull
3. Discretise it and use it with `epinow()` — it works exactly the same way as LogNormal or Gamma

This demonstrates the universality of the discretisation approach: you never need to worry about which distribution families your analysis tools support.

## Key points

- `EpiParameters.jl` → `discretise()` → `gt_opts()` / `delay_opts()` is the
  standard pipeline
- **Any distribution family** works — discretisation to PMF is universal
- **Compose delays** with `+` (convolution of PMFs)
- Use `UncertainDistribution` when delay parameters are themselves uncertain
- Always check the **source** of your delay distribution via the metadata
- The serial interval can inform contact tracing windows — use `cdf()` to calculate the probability of capturing infectors within a given time window
- Use `pdf`, `cdf`, `quantile`, and `rand` to extract summary statistics and informative values from delay distributions
- The CDF of the serial interval informs contact tracing window design
- The quantile function of the incubation period informs quarantine length
- Set `max` in `discretise()` to the 99th percentile of the continuous distribution


---

*Adapted from the [Epiverse-TRACE tutorials](https://epiverse-trace.github.io/tutorials/), © Epiverse-TRACE contributors, licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).*
