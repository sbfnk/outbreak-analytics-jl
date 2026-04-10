# EpiParameters.jl

A Julia package providing a curated database of epidemiological delay distributions from the published literature. Covers incubation periods, serial intervals, generation times, onset-to-death delays, offspring distributions, and more.

## Quick start

```julia
using EpiParameters
using Distributions

# Find COVID-19 incubation period estimates
results = epiparameter(disease="COVID-19", epi_name="incubation period")

# Get the fitted distribution
fitted = filter(p -> !isnothing(p.distribution), results)
d = fitted[1].distribution  # e.g. LogNormal(1.62, 0.42)

# Use standard Distributions.jl functions
mean(d)
quantile(d, 0.95)
rand(d, 1000)
```

## Querying

```julia
# By disease
epiparameter(disease="Ebola")

# By parameter type
epiparameter(epi_name="serial interval")

# By author
epiparameter(author="Linton")

# Combined (AND logic, case-insensitive)
epiparameter(disease="COVID-19", epi_name="incubation period", author="Lauer")
```

## Browsing

```julia
list_diseases()     # all diseases in the database
list_parameters()   # all parameter types
```

## Distribution families

The database includes LogNormal, Gamma, Weibull, Normal, NegativeBinomial, Geometric, and Poisson distributions. Since all are standard `Distributions.jl` types, they integrate directly with other packages:

```julia
using CFR
cfr_static(data; delay_density=d)

using EpiNow2
discretise(d; max=14)
```

## Metadata

Each `EpiParam` entry carries citation information, summary statistics, and method assessment flags:

```julia
p = results[1]
p.citation         # author, title, year, doi
p.summary_stats    # mean, sd, quantiles
p.metadata         # sample_size, region, units
p.method_assessment  # truncation, censoring flags
```
