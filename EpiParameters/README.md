# EpiParameters.jl

A Julia package providing a curated database of epidemiological delay distributions (incubation periods, serial intervals, onset-to-death, generation times, and more) from the published literature.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/sbfnk/outbreak-analytics-jl", subdir="EpiParameters")
```

## Quick start

```julia
using EpiParameters
using Distributions

# Find COVID-19 incubation period estimates
results = epiparameter(disease="COVID-19", epi_name="incubation period")

# Get the fitted distribution
d = results[1].distribution  # e.g. LogNormal(1.62, 0.42)
mean(d)       # mean incubation period
quantile(d, 0.95)  # 95th percentile
```

## Features

### Querying the database

```julia
# By disease
epiparameter(disease="Ebola")

# By parameter type
epiparameter(epi_name="serial interval")

# By author
epiparameter(author="Linton")

# Combined (AND logic)
epiparameter(disease="COVID-19", epi_name="incubation period", author="Lauer")
```

All filters are case-insensitive substring matches.

### Browsing available data

```julia
list_diseases()     # all diseases in the database
list_parameters()   # all parameter types
```

### Working with distributions

Parameters with fitted distributions return standard `Distributions.jl` objects, so all the usual functions work:

```julia
d = results[1].distribution

pdf(d, 5.0)         # probability density at day 5
cdf(d, 7.0)         # P(delay <= 7 days)
quantile(d, 0.99)   # 99th percentile
rand(d, 1000)       # random samples
mean(d)             # mean
std(d)              # standard deviation
```

Supported distribution families include LogNormal, Gamma, Weibull, Normal, NegativeBinomial, Geometric, and Poisson.

### Metadata and citations

Each `EpiParam` entry carries full metadata:

```julia
p = results[1]

p.disease          # "COVID-19"
p.epi_name         # "incubation period"
p.summary_stats    # Dict with mean, sd, quantiles, ...
p.citation         # Dict with author, title, year, doi
p.metadata         # Dict with sample_size, region, units, ...
p.method_assessment  # Dict with truncation, censoring flags
```

### Filtering for fitted distributions

Not all entries have fitted distributions. Filter for those that do:

```julia
fitted = filter(p -> !isnothing(p.distribution), results)
```

## Ecosystem integration

EpiParameters.jl provides inputs for other packages in the ecosystem:

```julia
using EpiParameters
using CFR

# Get onset-to-death delay for CFR estimation
delay = epiparameter(disease="COVID-19", epi_name="onset to death")
d = filter(p -> !isnothing(p.distribution), delay)[1].distribution
cfr_static(data; delay_density=d)
```

```julia
using EpiParameters
using EpiNow2

# Get generation time for Rt estimation
gt = epiparameter(disease="COVID-19", epi_name="serial interval")
d = filter(p -> !isnothing(p.distribution), gt)[1].distribution
pmf = discretise(truncated(d; lower=0); max=14)
```

## API reference

| Function | Description |
|----------|-------------|
| `epiparameter(; disease, epi_name, author)` | Query the database |
| `list_diseases()` | List all diseases |
| `list_parameters()` | List all parameter types |

## Related packages

- R equivalent: [epiparameter](https://github.com/epiverse-trace/epiparameter)
