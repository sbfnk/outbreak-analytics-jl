# Epidemics.jl

A Julia package for simulating age-structured SEIR epidemic models with interventions and vaccination. Wraps [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) with a clean teaching-oriented API.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/sbfnk/outbreak-analytics-jl", subdir="Epidemics")
```

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
I0[1] = 10.0  # seed 10 infections in youngest group

output = simulate(model, sym.matrix;
    demography = population,
    initial_infected = I0,
    tspan = (0.0, 365.0)
)
```

Returns a tidy DataFrame with columns `:time`, `:group`, `:S`, `:E`, `:I`, `:R`.

## Features

### SEIR model

The model parameters are specified as rates:

```julia
model = SEIR(
    beta = 0.025,   # transmission rate per contact per day
    sigma = 1/2,     # 1 / latent period (days)
    gamma = 1/5      # 1 / infectious period (days)
)
```

The SEIR differential equations for each age group *i* are:

$$\frac{dS_i}{dt} = -\beta S_i \sum_j C_{ij} \frac{I_j}{N_j}, \quad
\frac{dE_i}{dt} = \beta S_i \sum_j C_{ij} \frac{I_j}{N_j} - \sigma E_i, \quad
\frac{dI_i}{dt} = \sigma E_i - \gamma I_i, \quad
\frac{dR_i}{dt} = \gamma I_i$$

### Interventions

Model social distancing or other contact-reducing measures:

```julia
output = simulate(model, cm;
    demography = population,
    initial_infected = I0,
    tspan = (0.0, 365.0),
    interventions = [
        Intervention(time_begin=30.0, time_end=90.0, reduction=0.3),
        Intervention(time_begin=120.0, time_end=180.0, reduction=0.5),
    ]
)
```

The `reduction` parameter scales down the contact matrix during the intervention window (0.3 = 30% fewer contacts).

### Vaccination

Target vaccination campaigns to specific age groups:

```julia
output = simulate(model, cm;
    demography = population,
    initial_infected = I0,
    tspan = (0.0, 365.0),
    vaccination = [
        # Vaccinate children (group 1) at 3% per day for the first 30 days
        Vaccination(time_begin=0.0, time_end=30.0,
                    rate=0.03, groups=[true, false, false, false]),
    ]
)
```

Vaccination moves individuals from S to R at the specified rate.

### Works with ContactMatrix or raw matrices

```julia
# Using ContactMatrix (group labels preserved in output)
simulate(model, cm; demography=population, ...)

# Using a plain matrix
simulate(model, [10.0 3.0; 3.0 5.0]; demography=[1000.0, 2000.0], ...)
```

## API reference

| Type / Function | Description |
|----------------|-------------|
| `SEIR(; beta, sigma, gamma)` | SEIR model parameters |
| `Intervention(; time_begin, time_end, reduction)` | Contact reduction intervention |
| `Vaccination(; time_begin, time_end, rate, groups)` | Vaccination campaign |
| `simulate(model, cm; demography, initial_infected, tspan, interventions, vaccination, saveat)` | Run the simulation |

## Related packages

- [ContactMatrices.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/ContactMatrices) — contact matrix type
- [FinalSize.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/FinalSize) — analytical final size (complements dynamic simulation)
- R equivalent: [epidemics](https://github.com/epiverse-trace/epidemics)
