# CFR.jl

A Julia package for estimating case fatality ratios from epidemic data, with support for delay-adjusted and time-varying estimates.

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

# Time-varying CFR
cfr_time_varying(data; delay_density=delay)
```

## Delay adjustment

During a growing epidemic, the naive CFR (deaths / cases) underestimates the true value because recent cases have not yet had time to die. The delay-adjusted method (Nishiura et al. 2009) convolves case counts with the onset-to-death delay distribution to estimate the number of cases with known outcomes, correcting for this right-censoring bias.

## Input format

Data must be a DataFrame with columns:
- `:date` — date of observation
- `:cases` — number of new cases
- `:deaths` — number of new deaths

## References

Nishiura, H., Klinkenberg, D., Roberts, M., & Heesterbeek, J. A. P. (2009). Early epidemiological assessment of the virulence of emerging infectious diseases: a case study of an influenza pandemic. *PLoS ONE*, 4(8), e6852.
