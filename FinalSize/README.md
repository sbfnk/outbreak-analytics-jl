# FinalSize.jl

A Julia package for computing the final size of epidemics — the proportion of the population ultimately infected — given a basic reproduction number and (optionally) heterogeneous contact patterns and susceptibility.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/sbfnk/outbreak-analytics-jl", subdir="FinalSize")
```

## Quick start

### Homogeneous population

```julia
using FinalSize

# What proportion is infected when R₀ = 1.5?
final_size(1.5)  # ≈ 0.583
```

### Age-structured population

```julia
using FinalSize
using ContactMatrices
using SocialMixr

# Build a contact matrix from POLYMOD data
result = polymod() |>
    s -> filter_survey(s; countries=["United Kingdom"]) |>
    s -> assign_age_groups(s; age_limits=[0, 15, 45, 65]) |>
    compute_matrix

uk_pop = polymod_population(countries=["United Kingdom"])
sym = symmetrise(result, uk_pop)

population = pop_age(uk_pop, [0, 15, 45, 65]).population

# Final size by age group
fs = final_size(1.5, sym.matrix; demography=population)
```

Returns a DataFrame:

| group | susc_group | susceptibility | p_infected |
|-------|-----------|---------------|------------|
| [0,15) | 1 | 1.0 | 0.67 |
| [15,45) | 1 | 1.0 | 0.59 |
| ... | ... | ... | ... |

## Features

### Heterogeneous susceptibility

Model differential susceptibility across groups (e.g. partial immunity from prior exposure or vaccination):

```julia
# Two susceptibility classes per age group:
# 80% fully susceptible, 20% with 50% reduced susceptibility
susceptibility = [1.0 0.5; 1.0 0.5; 1.0 0.5; 1.0 0.5]
p_susceptibility = [0.8 0.2; 0.8 0.2; 0.8 0.2; 0.8 0.2]

fs = final_size(1.5, cm;
    demography = population,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility
)
```

### Works with ContactMatrix or raw matrices

```julia
# Using ContactMatrix (labels preserved in output)
final_size(1.5, cm; demography=population)

# Using a plain matrix
final_size(1.5, [10.0 3.0; 3.0 5.0]; demography=[1000.0, 2000.0])
```

## API reference

| Function | Description |
|----------|-------------|
| `final_size(R0)` | Homogeneous final size equation |
| `final_size(R0, cm; demography, susceptibility, p_susceptibility)` | Heterogeneous final size with contact matrix |

## Related packages

- [ContactMatrices.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/ContactMatrices) — contact matrix type
- [Epidemics.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/Epidemics) — dynamic SEIR simulation (complements static final size)
- R equivalent: [finalsize](https://github.com/epiverse-trace/finalsize)
