# Outbreak Analytics in Julia

Teaching materials for outbreak analytics using the Julia epidemiology ecosystem.
These tutorials are Julia translations of the
[Epiverse-TRACE tutorials](https://epiverse-trace.github.io/tutorials/),
covering the same material with Julia packages.

## Packages

| Package | Description | R equivalent |
|---------|-------------|--------------|
| [ContactMatrices.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/ContactMatrices) | Contact matrix type | contactmatrix |
| [SocialMixr.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/SocialMixr) | Social contact survey data | socialmixr |
| [FinalSize.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/FinalSize) | Final size of epidemics | finalsize |
| [CFR.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/CFR) | Case fatality ratio estimation | cfr |
| [EpiParameters.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/EpiParameters) | Epidemiological parameter database | epiparameter |
| [Epidemics.jl](https://github.com/sbfnk/outbreak-analytics-jl/tree/main/Epidemics) | Compartmental models (SEIR) | epidemics |
| [EpiNow2.jl](https://github.com/sbfnk/EpiNow2.jl) | Rt estimation, forecasting, secondary outcomes | EpiNow2 |
| [EpiBranch.jl](https://github.com/epiforecasts/EpiBranch.jl) | Branching processes, line lists, chain statistics | epichains, simulist |

## Tutorials

### Early tasks

```@contents
Pages = [
    "tutorials/read-case-data.md",
    "tutorials/clean-case-data.md",
    "tutorials/validate-case-data.md",
    "tutorials/aggregate-visualise.md",
]
```

### Middle tasks

```@contents
Pages = [
    "tutorials/epidemiological-parameters.md",
    "tutorials/outbreak-severity.md",
    "tutorials/quantifying-transmission.md",
    "tutorials/delays-in-analysis.md",
    "tutorials/short-term-forecast.md",
    "tutorials/superspreading.md",
    "tutorials/transmission-chains.md",
]
```

### Late tasks

```@contents
Pages = [
    "tutorials/contact-matrices.md",
    "tutorials/simulating-transmission.md",
    "tutorials/choosing-model.md",
    "tutorials/disease-burden.md",
]
```
