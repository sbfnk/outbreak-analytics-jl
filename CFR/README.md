# CFR.jl

A Julia package for estimating case fatality ratios (CFR) from epidemic data, with support for delay-adjusted and time-varying estimates.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/sbfnk/outbreak-analytics-jl", subdir="CFR")
```

## Quick start

```julia
using CFR
using DataFrames
using Distributions

data = DataFrame(
    date = Date(2020, 1, 1):Day(1):Date(2020, 3, 31),
    cases = [...],
    deaths = [...]
)

# Naive CFR
cfr_static(data)

# Delay-adjusted CFR (Nishiura et al. 2009)
delay = LogNormal(2.5, 0.8)
cfr_static(data; delay_density=delay)
```

## Features

### Static (overall) CFR

```julia
# Naive: deaths / cases
result = cfr_static(data)
# CfrEstimate(0.032, 95% CI: 0.028–0.037)

# Delay-adjusted: corrects for right-censoring
result = cfr_static(data; delay_density=LogNormal(2.5, 0.8))
```

The delay-adjusted method (Nishiura et al. 2009) accounts for cases whose outcomes are not yet known. During a growing epidemic, the naive CFR underestimates the true value because recent cases have not yet had time to die. The adjustment convolves case counts with the onset-to-death delay distribution to estimate cases with known outcomes.

### Time-varying CFR

Track how the CFR evolves over the course of an outbreak:

```julia
tv = cfr_time_varying(data; delay_density=LogNormal(2.5, 0.8))
# DataFrame with columns: date, cfr, lower, upper, known_outcomes
```

Options:
- `burn_in=7` — discard the first *n* time points (unstable with few cases)
- `smoothing_window=7` — rolling sum of cases and deaths before estimation

### Estimating known outcomes

The convolution step is available separately:

```julia
data_with_outcomes = estimate_outcomes(data, LogNormal(2.5, 0.8))
# Adds :known_outcomes column
```

### Confidence intervals

All estimates include exact binomial (Clopper-Pearson) confidence intervals:

```julia
result = cfr_static(data; confidence_level=0.99)
result.estimate  # point estimate
result.lower     # lower bound
result.upper     # upper bound
```

## API reference

| Function | Description |
|----------|-------------|
| `cfr_static(data; delay_density, confidence_level)` | Overall CFR estimate |
| `cfr_time_varying(data; delay_density, burn_in, smoothing_window)` | Time-varying CFR |
| `estimate_outcomes(data, delay_density)` | Known outcomes via convolution |

## Input format

Data must be a DataFrame with columns:
- `:date` — date of observation
- `:cases` — number of new cases
- `:deaths` — number of new deaths

## Related packages

- [EpiParameters.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/EpiParameters) — retrieve onset-to-death delay distributions from the literature
- R equivalent: [cfr](https://github.com/epiverse-trace/cfr)

## References

Nishiura, H., Klinkenberg, D., Roberts, M., & Heesterbeek, J. A. P. (2009). Early epidemiological assessment of the virulence of emerging infectious diseases: a case study of an influenza pandemic. *PLoS ONE*, 4(8), e6852.
