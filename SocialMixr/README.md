# SocialMixr.jl

A Julia package for processing social contact survey data into contact matrices for epidemiological modelling. Julia port of the R [socialmixr](https://github.com/epiforecasts/socialmixr) package.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/sbfnk/outbreak-analytics-jl", subdir="SocialMixr")
```

## Quick start

```julia
using SocialMixr
using ContactMatrices

# Load POLYMOD data, filter to UK, compute a contact matrix
result = polymod() |>
    s -> filter_survey(s; countries=["United Kingdom"]) |>
    s -> assign_age_groups(s; age_limits=[0, 5, 15, 45, 65]) |>
    compute_matrix

cm = result.matrix  # a ContactMatrix
```

## Features

### Bundled data

The package ships with the [POLYMOD](https://doi.org/10.1371/journal.pmed.0050074) contact survey (Mossong et al. 2008) and World Population Prospects demographic data:

```julia
survey = polymod()                                        # ContactSurvey
pop = polymod_population(countries=["United Kingdom"])     # DataFrame
```

### Processing pipeline

The typical workflow mirrors the R package, using Julia's `|>` pipe:

```julia
result = polymod() |>
    s -> filter_survey(s; countries=["Germany"]) |>
    s -> assign_age_groups(s; age_limits=[0, 5, 15, 45, 65]) |>
    compute_matrix
```

**1. Load survey data** — `polymod()` returns a `ContactSurvey` containing participant and contact DataFrames.

**2. Filter** — `filter_survey(survey; countries=["..."])` selects participants by country. Custom filters on contact-level columns are also supported.

**3. Assign age groups** — `assign_age_groups(survey; age_limits=[0, 15, 65])` bins participant and contact ages into groups. Handles missing ages and age ranges (e.g. "estimated minimum/maximum age" columns) via imputation.

**4. Compute matrix** — `compute_matrix(survey)` returns a named tuple with a `ContactMatrix` and the processed participant DataFrame.

### Weighting

Apply survey weights (e.g. to correct for day-of-week sampling bias):

```julia
weighted = weigh(survey, :dayofweek)
```

### Post-processing

```julia
uk_pop = polymod_population(countries=["United Kingdom"])

# Symmetrise contacts using population data
sym_result = symmetrise(result, uk_pop)

# Adjust population to match age groups
pop_df = pop_age(uk_pop, [0, 5, 15, 45, 65])

# Per-capita contact rates
pc = per_capita(result, uk_pop)

# Decompose into mean contacts, normalisation, and assortativity
decomposed = split_matrix(result, uk_pop)
```

## API reference

| Function | Description |
|----------|-------------|
| `polymod()` | Load bundled POLYMOD survey |
| `polymod_population(; countries)` | Load WPP population data |
| `load_population(path)` | Load population from CSV |
| `filter_survey(survey; countries, filter)` | Filter participants |
| `assign_age_groups(survey; age_limits, ...)` | Bin ages into groups |
| `weigh(survey, by; target, groups)` | Apply survey weights |
| `compute_matrix(survey; counts)` | Compute contact matrix |
| `symmetrise(result, population)` | Symmetrise matrix |
| `pop_age(pop, age_limits)` | Adjust population to age groups |
| `per_capita(result, population)` | Per-capita contact rates |
| `split_matrix(result, population)` | Decompose matrix |

## Related packages

- [ContactMatrices.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/ContactMatrices) — the `ContactMatrix` type this package produces
- R equivalent: [socialmixr](https://github.com/epiforecasts/socialmixr)
