# Epidemics.jl

A Julia package for simulating age-structured SEIR epidemic models with interventions and vaccination. Wraps [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) with a clean teaching-oriented API.

## Quick start

```julia
using Epidemics
using ContactMatrices
using SocialMixr

# Set up contact matrix and population
result = polymod() |>
    s -> filter_survey(s; countries=["United Kingdom"]) |>
    s -> assign_age_groups(s; age_limits=[0, 15, 45, 65]) |>
    compute_matrix

uk_pop = polymod_population(countries=["United Kingdom"])
sym = symmetrise(result, uk_pop)
population = pop_age(uk_pop, [0, 15, 45, 65]).population

# Define and run the model
model = SEIR(beta=0.025, sigma=1/2, gamma=1/5)

I0 = zeros(length(population))
I0[1] = 10.0

output = simulate(model, sym.matrix;
    demography = population,
    initial_infected = I0,
    tspan = (0.0, 365.0)
)
```

Returns a tidy DataFrame with columns `:time`, `:group`, `:S`, `:E`, `:I`, `:R`.

## Interventions

Model social distancing or other contact-reducing measures:

```julia
simulate(model, cm;
    demography = population,
    initial_infected = I0,
    tspan = (0.0, 365.0),
    interventions = [
        Intervention(time_begin=30.0, time_end=90.0, reduction=0.3),
    ]
)
```

## Vaccination

Target vaccination campaigns to specific age groups:

```julia
simulate(model, cm;
    demography = population,
    initial_infected = I0,
    tspan = (0.0, 365.0),
    vaccination = [
        Vaccination(time_begin=0.0, time_end=30.0,
                    rate=0.03, groups=[true, false, false, false]),
    ]
)
```
