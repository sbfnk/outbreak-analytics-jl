# FinalSize.jl

A Julia package for computing the final size of epidemics — the proportion of the population ultimately infected — given a basic reproduction number and (optionally) heterogeneous contact patterns and susceptibility.

## Quick start

### Homogeneous population

```julia
using FinalSize

final_size(1.5)  # ≈ 0.583
```

### Age-structured population

```julia
using FinalSize
using ContactMatrices
using SocialMixr

result = polymod() |>
    s -> filter_survey(s; countries=["United Kingdom"]) |>
    s -> assign_age_groups(s; age_limits=[0, 15, 45, 65]) |>
    compute_matrix

uk_pop = polymod_population(countries=["United Kingdom"])
sym = symmetrise(result, uk_pop)
population = pop_age(uk_pop, [0, 15, 45, 65]).population

fs = final_size(1.5, sym.matrix; demography=population)
```

### Heterogeneous susceptibility

```julia
# Two susceptibility classes: 80% fully susceptible, 20% with 50% reduction
susceptibility = [1.0 0.5; 1.0 0.5; 1.0 0.5; 1.0 0.5]
p_susceptibility = [0.8 0.2; 0.8 0.2; 0.8 0.2; 0.8 0.2]

fs = final_size(1.5, cm;
    demography = population,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility
)
```
