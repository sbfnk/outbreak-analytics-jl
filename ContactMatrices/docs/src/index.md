# ContactMatrices.jl

A Julia package for representing and manipulating contact matrices used in epidemiological modelling.

`ContactMatrix` is a shared type designed to be dispatched on across the Julia epidemiology ecosystem — used by [SocialMixr.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/SocialMixr), [FinalSize.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/FinalSize), and [Epidemics.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/Epidemics).

## Quick start

```julia
using ContactMatrices

# From a matrix and labels
cm = ContactMatrix([10.0 3.0; 3.0 5.0], ["0-14", "15+"]; setting=:all)

# From tidy-format data
cm = ContactMatrix(
    of   = ["0-14", "0-14", "15+", "15+"],
    with = ["0-14", "15+", "0-14", "15+"],
    value = [10.0, 3.0, 3.0, 5.0]
)

# Symmetrise using population data
cm_sym = make_symmetric(cm, [5e6, 10e6])

# Combine settings
cm_total = cm_home + cm_work

# Extract the raw matrix
Matrix(cm)
```

## Ecosystem integration

```julia
using FinalSize
final_size(1.5, cm; demography=population)

using Epidemics
simulate(SEIR(beta=0.025, sigma=0.5, gamma=0.2), cm; demography=population, ...)
```
