# Estimation of Outbreak Severity

During an outbreak, a key question is: **how severe is this disease?** The case
fatality ratio (CFR) — the proportion of confirmed cases who die — is a critical
measure for public health response.

However, calculating CFR during an ongoing epidemic is not as simple as dividing
deaths by cases. This tutorial covers:

1. The **naive CFR** and why it is biased
2. **Delay-adjusted CFR** using the method of Nishiura et al. (2009)
3. **Time-varying CFR** for monitoring severity over the course of an outbreak

Common questions at the start of an epidemic include: What is the likely public health impact? Which groups are most severely affected? Does this outbreak have pandemic potential? Frameworks such as the CDC's Pandemic Severity Assessment Framework consider both transmissibility and clinical severity when assessing pandemic potential.

The CFR is defined as the conditional probability of death given confirmed diagnosis — the ratio of cumulative deaths to cumulative confirmed cases. However, calculating this directly during an active epidemic yields biased estimates due to time delays between onset and death. Values vary substantially early on and stabilise only in later outbreak stages.

Estimating severity is also valuable beyond pandemic planning. Knowing whether an outbreak differs from historical severity patterns can motivate investigation into whether the difference is intrinsic to the pathogen (e.g. a new strain) or population-based (e.g. reduced immunity or increased comorbidities).

### Data sources for severity

Different data sources capture different levels of severity. Severe and critical cases are identified in hospital settings. Symptomatic cases may be found through surveillance systems. Mild and asymptomatic cases are typically only identified through contact tracing or serological surveys. The CFR only captures severity among *confirmed* cases, so it will typically exceed the infection fatality risk (IFR), which includes all infections.

## The data: Ebola 1976

We use data from the first known Ebola outbreak in Yambuku, Zaire (1976):

```@example severity
using CFR
using EpiParameters
using Distributions
using DataFrames
using Dates
using CSV
using CairoMakie

ebola = CSV.read(joinpath(@__DIR__, "ebola1976.csv"), DataFrame)
first(ebola, 5)
```

!!! details "Coming from R?"
    Julia uses `CSV.read(file, DataFrame)` where R uses `read.csv(file)` or `readr::read_csv(file)`. The `DataFrame` argument tells Julia what type to return — Julia's type system means the same `CSV.read` function can produce different output types.

```@example severity
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
```

```@example severity
println("Total cases: $(sum(ebola.cases)), total deaths: $(sum(ebola.deaths))")
```

!!! details "Coming from R?"
    Julia uses `$()` inside strings for interpolation: `"Total: $(sum(x))"` is equivalent to R's `glue::glue("Total: {sum(x)}")` or `paste0("Total: ", sum(x))`. For concatenating strings, Julia uses the `*` operator: `"hello" * " " * "world"`. This surprises R users who expect `paste0()`, but the convention comes from mathematics (string concatenation is a monoid operation, like multiplication).

!!! details "Coming from R?"
    Julia's `CairoMakie` replaces `ggplot2` for plotting. The key difference is that Makie builds plots imperatively (`lines!()`, `scatter!()` add to an existing axis) rather than declaratively (ggplot2's `+` layers). The `!` in function names like `barplot!()` indicates the function modifies an existing object.

## Naive CFR

The naive (or crude) CFR at time $t$ is calculated as:

$$b_t = \frac{D_t}{C_t}$$

where $D_t$ is cumulative deaths and $C_t$ is cumulative confirmed cases at time $t$.

The simplest estimate divides total deaths by total cases:

```@example severity
naive_cfr = cfr_static(ebola)
```

```@example severity
println("Naive CFR: $(round(naive_cfr.estimate * 100; digits=1))% " *
        "(95% CI: $(round(naive_cfr.lower * 100; digits=1))–$(round(naive_cfr.upper * 100; digits=1))%)")
```

During an ongoing epidemic, this **underestimates** the true CFR because cases
confirmed recently haven't had time to die yet.

The magnitude of this bias depends on two factors: the epidemic growth rate and the onset-to-death delay distribution. Faster growth and longer delays create greater bias, because a larger proportion of recent cases will not yet have had time for their outcome to be observed.

### Sources of bias in CFR estimation

There are two main biases that affect CFR estimation (Lipsitch et al. 2015):

**1. Preferential ascertainment of severe cases:** For diseases with a spectrum of clinical presentations, cases reaching public health authorities typically represent individuals with the most severe symptoms — those seeking medical care, being hospitalised, or dying. The CFR among detected cases will therefore typically exceed the CFR across the entire infected population, which includes mild, subclinical, and asymptomatic infections.

**2. Delayed reporting of deaths:** During ongoing epidemics, there is a time lag between case confirmation and death. At any point, the case list includes people who will die but whose deaths have not yet occurred or been reported. Dividing reported deaths by reported cases at a given time point underestimates the true CFR. The magnitude of this bias depends on the epidemic growth rate and the delay distribution from case reporting to death — longer delays and faster growth create greater underestimation.

This tutorial focuses on addressing the second bias using `CFR.jl`.

### Case study: Influenza A (H1N1), Mexico, 2009

The importance of delay-adjusted CFR estimation was starkly illustrated during the 2009 H1N1 pandemic. Early in the outbreak, Mexico's CFR estimates appeared alarmingly high, suggesting a highly virulent virus. Meanwhile, other countries perceived the same virus as relatively mild. Much of this discrepancy arose from differences in case ascertainment and the right-censoring bias — in countries with few deaths in the first days following identification, naive CFR estimates were misleadingly low, whilst in Mexico, where more time had elapsed and deaths had accumulated, they appeared higher. Having methods that correct for reporting delays enables more consistent and reliable cross-country severity comparisons.

## Delay-adjusted CFR

The method of **Nishiura et al. (2009)** corrects for this bias using the
**onset-to-death delay distribution**:

The key insight is that at any point during an outbreak, only a fraction of confirmed cases have had enough time for their outcome (death or recovery) to be observed. By estimating this fraction using the onset-to-death delay distribution, we can correct the denominator to include only cases with known outcomes.

In real-time outbreaks, we may not have sufficient data to estimate the onset-to-death distribution directly. Instead, we use estimates from historical outbreaks or literature databases like `EpiParameters.jl`.

```@example severity
ebola_otd = epiparameter(disease="Ebola", epi_name="onset to death")
delay_param = filter(p -> !isnothing(p.distribution), ebola_otd)[1]
delay = delay_param.distribution
```

```@example severity
adjusted_cfr = cfr_static(ebola; delay_density=delay)
```

### How does the delay adjustment work?

The method relates cumulative deaths ($D_t$) to daily case incidence ($c_s$) and the probability density function of the onset-to-death delay ($f$):

$$D_t = p \sum_{s=0}^{t} c_s F(t - s)$$

where $p$ is the true (unbiased) CFR and $F(t-s)$ is the cumulative distribution function of the onset-to-death delay evaluated at time $t-s$. The term $\sum_s c_s F(t-s)$ represents the total expected number of cases with known outcomes by time $t$ — each case is weighted by the probability that their outcome has been observed given the time since their onset.

The method then estimates $p$ using maximum likelihood, treating the observed deaths as a binomial sample from the cases with known outcomes.

The delay-adjusted CFR is higher than the naive estimate, because it accounts
for cases whose outcomes were not yet known.

The delay-adjusted CFR is higher than the naive estimate because it correctly accounts for cases whose outcomes were not yet known at the time of analysis.

## Time-varying CFR

We can track how the CFR estimate evolves over time:

The `cfr_time_varying` function calculates what the CFR estimate would have been on each day of the outbreak, using only data available up to that point. The burn-in period (14 days here) excludes the very early phase when few cases and deaths make estimates highly unstable.

```@example severity
tv_naive = cfr_time_varying(ebola; burn_in=14)
tv_adjusted = cfr_time_varying(ebola; delay_density=delay, burn_in=14)

days = 15:nrow(ebola)
valid = .!isnan.(tv_adjusted.cfr)

fig = Figure(size=(700, 400))
ax = Axis(fig[1, 1];
          xlabel="Day of outbreak", ylabel="CFR",
          title="Time-varying CFR — naive vs delay-adjusted")

band!(ax, days, tv_naive.lower, tv_naive.upper; color=(:steelblue, 0.2))
lines!(ax, days, tv_naive.cfr; linewidth=2, color=:steelblue, label="Naive")

band!(ax, days[valid], tv_adjusted.lower[valid], tv_adjusted.upper[valid];
      color=(:firebrick, 0.2))
lines!(ax, days[valid], tv_adjusted.cfr[valid];
       linewidth=2, color=:firebrick, label="Delay-adjusted")

axislegend(ax; position=:rb)
fig
```

The delay-adjusted estimate converges to the true CFR earlier than the naive estimate. This is crucial for real-time decision making — during the 2009 H1N1 pandemic, for example, early CFR estimates varied dramatically between countries, partly because of this right-censoring bias. Methods that correct for reporting delays enable earlier and more reliable severity assessment.

!!! note "Rolling vs time-varying CFR"
    The time-varying CFR computed by `cfr_time_varying` shows what the CFR estimate would have been on each day using all data available up to that point. This is useful for assessing how quickly estimates stabilise. It differs from a *moving-window* approach, which would compute the CFR over a sliding window of fixed width — that approach is more sensitive to changes in severity over time (e.g. due to new variants or changes in treatment) but requires more data per window.

### Case study: SARS, Hong Kong, 2003

Nishiura et al. (2009) demonstrated the delay-adjusted method using SARS outbreak data from Hong Kong in 2003. Naive CFR estimates during the outbreak substantially underestimated the realised CFR at the outbreak's conclusion (302/1755 = 17.2%). However, the delay-adjusted CFR computed as early as 27 March 2003 was 18.1% (95% CI: 10.5–28.1%) — an overestimate at that early stage, but one whose 95% confidence interval already included the eventual realised value. By contrast, the naive estimate did not reach the true value until much later in the outbreak. This demonstrates the practical value of delay adjustment for early-stage severity assessment.

## Sensitivity to the delay distribution

The choice of onset-to-death delay distribution is one of the most important modelling decisions in CFR estimation. If the assumed delay is too short, the adjusted estimate will be closer to the naive (biased) value. If too long, it may overestimate severity. This underscores the importance of using well-characterised delay distributions from the published literature.

What happens if we assume a shorter or longer onset-to-death delay?

```@example severity
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

hlines!(ax, [naive_cfr.estimate]; color=:grey, linestyle=:dash, label="Naive CFR")
axislegend(ax; position=:rb)
fig
```

## Beyond CFR: other severity measures

The same delay-adjustment framework can be applied to other severity measures by changing the numerator and denominator:

| Measure | Numerator | Denominator | Delay distribution |
|---|---|---|---|
| **Case Fatality Risk (CFR)** | Deaths | Confirmed cases | Onset-to-death |
| **Infection Fatality Risk (IFR)** | Deaths | All infections | Exposure-to-death |
| **Hospitalisation Fatality Risk (HFR)** | Deaths | Hospitalisations | Hospitalisation-to-death |

To estimate the IFR, for example, you would need data on total infections (not just confirmed cases), which typically requires serological surveys or capture-recapture methods. The HFR uses hospitalisation data and the delay from hospital admission to death.

CFR may also differ across populations — by age, geography, or treatment regimen. Quantifying these heterogeneities helps target resources appropriately and compare care quality across settings.

## Key points

- The **naive CFR** (deaths/cases) underestimates true severity during an ongoing
  outbreak because of the delay between onset and death
- **Delay-adjusted CFR** (Nishiura et al. 2009) corrects for right-censoring
- `EpiParameters.jl` provides delay distributions that plug directly into `CFR.jl`
- **Time-varying CFR** tracks how severity estimates evolve over time
- The estimate is **sensitive** to the assumed delay distribution — always report
  this choice and consider sensitivity analyses
- CFR estimates severity among *confirmed cases* only — the infection fatality risk (IFR) would require data on all infections, including undetected ones
- Use `cfr_time_varying` to assess how quickly severity estimates stabilise during an outbreak
- The same delay-adjustment framework applies to other severity measures (IFR, HFR) by changing the input data and delay distribution
- Severity may vary across subgroups — consider stratified analyses


---

*Adapted from the [Epiverse-TRACE tutorials](https://epiverse-trace.github.io/tutorials/), © Epiverse-TRACE contributors, licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).*
