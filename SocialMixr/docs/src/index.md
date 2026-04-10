# SocialMixr.jl

A Julia package for processing social contact survey data into contact matrices for epidemiological modelling. Julia port of the R [socialmixr](https://github.com/epiforecasts/socialmixr) package.

## Quick start

```julia
using SocialMixr
using ContactMatrices

# Load POLYMOD, filter to UK, compute contact matrix
result = polymod() |>
    s -> filter_survey(s; countries=["United Kingdom"]) |>
    s -> assign_age_groups(s; age_limits=[0, 5, 15, 45, 65]) |>
    compute_matrix

cm = result.matrix  # ContactMatrix

# Symmetrise with population data
uk_pop = polymod_population(countries=["United Kingdom"])
sym = symmetrise(result, uk_pop)
```

## Pipeline

The typical workflow mirrors the R package:

1. **Load** — `polymod()` returns a `ContactSurvey`
2. **Filter** — `filter_survey(survey; countries=["..."])` selects by country
3. **Age groups** — `assign_age_groups(survey; age_limits=[0, 15, 65])` bins ages
4. **Compute** — `compute_matrix(survey)` returns a `ContactMatrix`
5. **Post-process** — `symmetrise`, `per_capita`, `split_matrix`

## Bundled data

The package ships with the [POLYMOD](https://doi.org/10.1371/journal.pmed.0050074) contact survey (Mossong et al. 2008) and World Population Prospects demographic data for POLYMOD countries.
