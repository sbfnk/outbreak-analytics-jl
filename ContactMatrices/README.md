# ContactMatrices.jl

A Julia package for representing and manipulating contact matrices used in epidemiological modelling. Provides a shared `ContactMatrix` type that other packages in the ecosystem can dispatch on.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/sbfnk/outbreak-analytics-jl", subdir="ContactMatrices")
```

## Quick start

```julia
using ContactMatrices

# Construct from a matrix and age group labels
cm = ContactMatrix(
    [10.0 3.0; 3.0 5.0],
    ["0-14", "15+"]
)

# Construct from tidy-format data
cm = ContactMatrix(
    of = ["0-14", "0-14", "15+", "15+"],
    with = ["0-14", "15+", "0-14", "15+"],
    value = [10.0, 3.0, 3.0, 5.0]
)
```

## Features

### Core type

`ContactMatrix{T, N}` is an `AbstractArray` with grouping metadata and a setting label (e.g. `:home`, `:work`, `:school`, `:all`). It supports standard array operations including indexing, iteration, and broadcasting.

```julia
cm = ContactMatrix([10.0 3.0; 3.0 5.0], ["0-14", "15+"]; setting=:home)
setting(cm)     # :home
groupings(cm)   # Groupings with labels ["0-14", "15+"]
Matrix(cm)      # Extract the underlying matrix
```

### Symmetrisation

Ensure reciprocity of contacts using population data:

```julia
cm_sym = make_symmetric(cm, [1000.0, 2000.0])
```

This applies the formula `c_sym[i,j] = (N_i * c[i,j] + N_j * c[j,i]) / (2 * N_i)`, ensuring that the total number of contacts from group *i* to group *j* equals the reverse.

### Group aggregation

Combine fine age groups into coarser ones:

```julia
mapping = Dict("0-4" => "0-14", "5-9" => "0-14", "10-14" => "0-14",
               "15-19" => "15+", "20-24" => "15+")
cm_coarse = reduce_groups(cm, mapping; population=pop)
```

### Arithmetic

```julia
cm_total = cm_home + cm_work + cm_school   # combine settings
cm_reduced = 0.5 * cm                      # scale contact rates
```

## Ecosystem integration

`ContactMatrix` is designed as a shared type across Julia epidemiology packages:

```julia
using FinalSize
final_size(1.5, cm; demography=population)

using Epidemics
simulate(SEIR(beta=0.025, sigma=0.5, gamma=0.2), cm; demography=population, ...)
```

## API reference

| Function | Description |
|----------|-------------|
| `ContactMatrix(matrix, labels; setting)` | Construct from a matrix and labels |
| `ContactMatrix(; of, with, value, setting)` | Construct from tidy-format vectors |
| `setting(cm)` | Get the setting label |
| `groupings(cm)` | Get the `Groupings` metadata |
| `ndimgroups(g)` | Number of grouping dimensions |
| `make_symmetric(cm, population)` | Symmetrise using population weights |
| `reduce_groups(cm, mapping; population)` | Aggregate to coarser groups |
| `Matrix(cm)` | Extract the underlying matrix |

## Related packages

- [SocialMixr.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/SocialMixr) — compute contact matrices from survey data
- R equivalent: [contactmatrix](https://github.com/epiverse-trace/contactmatrix)
