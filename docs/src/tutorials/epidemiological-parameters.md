# Access Epidemiological Delay Distributions

In outbreak analytics, we frequently need **epidemiological delay distributions** —
probability distributions describing the time between key events in the course of
an infection (e.g. time from exposure to symptom onset, or from onset to death).

**EpiParameters.jl** provides a curated database of previously estimated parameters
that can be queried programmatically and used directly as `Distributions.jl` objects.

Infectious diseases follow an infection cycle, which usually includes a presymptomatic period, a symptomatic period, and a recovery period. The time periods describing these stages can be used to understand transmission dynamics and inform disease prevention and control interventions. However, early in an epidemic, efforts to understand the outbreak can be delayed by the lack of an easy way to access key parameters. Manually searching the literature for these values is slow and error-prone. `EpiParameters.jl` addresses this by providing a curated, queryable database of parameters from the published literature.

## The problem

Suppose we want to estimate the transmissibility of an infection from case data using a tool like EpiNow2.jl. To do this, we need delay distributions as inputs — for example, the serial interval or incubation period. Without a tool like `EpiParameters.jl`, we would have to manually search the published literature, identify the right paper for our disease and setting, locate the fitted distribution parameters in the text or supplementary material, and then transcribe them into our code by hand.

This manual workflow is slow, error-prone, and not reproducible. Different analysts may choose different source papers, copy parameters incorrectly, or fail to record which study their estimates came from. When the analysis needs to be updated or reviewed, retracing these steps is difficult.

Our goal in this tutorial is to access a specific set of epidemiological parameters **programmatically** — querying a curated database rather than extracting values by hand — and then to show how the retrieved distributions can be plugged directly into downstream analysis tools.

## Learning objectives

- Query the epidemiological parameter database
- Extract distribution objects for use in analysis
- Visualise and compare delay distributions

## Generation time vs serial interval

Two commonly confused delay distributions:

- **Generation time**: time from infection of a primary case to infection of a
  secondary case (infection → infection)
- **Serial interval**: time from symptom onset of a primary case to symptom onset
  of a secondary case (symptoms → symptoms)

The serial interval can be **negative** for diseases with pre-symptomatic
transmission (e.g. COVID-19), while the generation time is always positive.

The generation time, jointly with the reproduction number (R), provides valuable insight into the likely growth rate of an epidemic. The larger the value of R and/or the shorter the generation time, the more new infections we would expect per unit of time, and hence the faster the incidence of cases will grow.

In calculating the effective reproduction number, the generation time distribution is often approximated by the serial interval distribution. This approximation is used because it is easier to observe and record the onset of symptoms than the exact time of infection. However, this approximation is most appropriate for diseases in which infectiousness starts after symptom onset.

### From mean delays to probability distributions

If we measure the serial interval in real data, we typically see that not all case pairs have the same delay. To summarise this variability, we fit statistical distributions (e.g. log-normal, gamma, Weibull) to the observed data. These distributions are characterised by their summary statistics (mean, standard deviation) or their distribution parameters (shape, scale). The fitted distribution captures the full range of variability, rather than just a single average value.

## Exploring the database

The `epiparameter()` function works as a filtering function. We can query by disease name and parameter type. Note that some database entries have summary statistics but no fitted probability distribution — we need to filter for entries with fitted distributions to get `Distributions.jl` objects we can use in analysis.

Let's start by seeing what diseases and parameter types are available:

```@example epiparams
using EpiParameters
using Distributions
using CairoMakie

list_diseases()
```

```@example epiparams
list_parameters()
```

!!! details "Coming from R?"
    `using X` in Julia is equivalent to `library(X)` in R — it loads a package and makes its exports available. Unlike R, Julia packages are installed separately using `Pkg.add("X")` (similar to `install.packages("X")`), and the `using` statement only loads them.

## Querying for a specific disease

The database may contain multiple entries for the same disease and parameter type from different studies. These entries may differ in distribution family, parameter values, sample size, and study population. We filter for entries with fitted distributions:

Let's find all COVID-19 incubation periods with fitted distributions:

```@example epiparams
covid_ip = epiparameter(disease="COVID-19", epi_name="incubation period")
```

```@example epiparams
covid_ip_fitted = filter(p -> !isnothing(p.distribution), covid_ip)
```

!!! details "Coming from R?"
    `filter(p -> !isnothing(p.distribution), covid_ip)` is Julia's equivalent of R's `Filter()` or `dplyr::filter()`. The `p -> ...` syntax creates an anonymous function, like `\(p) ...` or `function(p) ...` in R. `isnothing()` is Julia's equivalent of `is.null()`.

## Visualising a distribution

The distribution is a standard `Distributions.jl` object — we can plot its
probability density function:

```@example epiparams
best = covid_ip_fitted[1]
(mean=mean(best.distribution), std=std(best.distribution), median=median(best.distribution))
```

```@example epiparams
d = best.distribution
x = range(0, quantile(d, 0.99); length=200)
fig = Figure(size=(600, 350))
ax = Axis(fig[1, 1];
          xlabel="Days", ylabel="Density",
          title="COVID-19 incubation period ($(typeof(d).name.name))")
lines!(ax, x, pdf.(d, x); linewidth=2.5)
vlines!(ax, [mean(d)]; color=:red, linestyle=:dash, linewidth=1.5,
        label="Mean = $(round(mean(d); digits=1))d")
axislegend(ax; position=:rt)
fig
```

!!! details "Coming from R?"
    Notice `pdf.(d, x)` — the dot before the parentheses is **broadcasting**, Julia's way of applying a function element-wise over an array. In R, most functions like `dnorm(x, ...)` are automatically vectorised. In Julia, you must opt in by adding `.` to the function call. This applies to operators too: `x .+ 1`, `x .== 0`, `.!mask`. Once you get used to it, the explicitness is helpful — you always know when an operation is element-wise.

We can also plot the **cumulative distribution function (CDF)**, which shows the probability of the delay being less than or equal to a given value. This is useful for answering questions like "what proportion of cases develop symptoms within 7 days?":

```@example epiparams
fig = Figure(size=(600, 350))
ax = Axis(fig[1, 1];
          xlabel="Days", ylabel="Cumulative probability",
          title="COVID-19 incubation period — CDF")
lines!(ax, x, cdf.(d, x); linewidth=2.5, color=:steelblue)
hlines!(ax, [0.95]; color=:grey, linestyle=:dash, linewidth=1)
text!(ax, 0.5, 0.96; text="95%", fontsize=10, color=:grey)
fig
```

## Extracting summary statistics

Once we have a distribution, we can extract summary statistics and distribution parameters for use in downstream analysis tools. In Julia, these come from the standard `Distributions.jl` interface:

```@example epiparams
println("Mean: $(round(mean(d); digits=2)) days")
println("Std:  $(round(std(d); digits=2)) days")
println("Median: $(round(median(d); digits=2)) days")
println("95th percentile: $(round(quantile(d, 0.95); digits=2)) days")
```

The distribution parameters (e.g. `meanlog` and `sdlog` for a log-normal) can be accessed directly:

```@example epiparams
println("Distribution type: $(typeof(d).name.name)")
println("Parameters: $(params(d))")
```

!!! details "Coming from R?"
    R uses function families for each distribution: `rnorm`/`dnorm`/`pnorm`/`qnorm` for random draws, density, CDF, and quantiles. In Julia, you create a single distribution object and use generic functions:

    | Operation | R | Julia |
    |-----------|---|-------|
    | Random draw | `rnorm(10, 5, 2)` | `rand(Normal(5, 2), 10)` |
    | Density | `dnorm(x, 5, 2)` | `pdf(Normal(5, 2), x)` |
    | CDF | `pnorm(x, 5, 2)` | `cdf(Normal(5, 2), x)` |
    | Quantile | `qnorm(p, 5, 2)` | `quantile(Normal(5, 2), p)` |

    The same pattern works for all distributions: `LogNormal`, `Gamma`, `Weibull`, etc.

These values can be plugged directly into analysis tools like EpiNow2.jl — see the [Using Delay Distributions in Analysis](@ref) tutorial for the complete pipeline.

## Comparing distributions across diseases

Let's compare incubation periods for different diseases:

```@example epiparams
diseases = ["COVID-19", "Influenza", "SARS", "Ebola Virus Disease", "Mpox"]
fig = Figure(size=(700, 400))
ax = Axis(fig[1, 1];
          xlabel="Days", ylabel="Density",
          title="Incubation periods by disease")
x = range(0, 25; length=300)

for disease in diseases
    results = epiparameter(disease=disease, epi_name="incubation period")
    fitted = filter(p -> !isnothing(p.distribution), results)
    isempty(fitted) && continue
    d = fitted[1].distribution
    lines!(ax, x, pdf.(d, x); linewidth=2,
           label="$(disease) (μ=$(round(mean(d); digits=1))d)")
end

axislegend(ax; position=:rt)
fig
```

Notice the range of incubation periods across diseases. COVID-19 and influenza have relatively short incubation periods (a few days), whilst Ebola and Mpox have longer ones (up to two weeks or more). These differences have important implications for surveillance and control: diseases with shorter incubation periods leave less time for contact tracing and isolation to be effective.

## Citation and metadata

Always check the source when selecting a parameter for your analysis. The reliability of an estimate depends on the sample size, the study population (region, age groups), and the inference method used. A parameter estimated from a small convenience sample may not be appropriate for a different setting:

```@example epiparams
best.citation["year"], best.citation["doi"]
```

```@example epiparams
best.metadata
```

## Onset-to-death distributions

These are critical for CFR estimation:

```@example epiparams
results = epiparameter(disease="Ebola", epi_name="onset to death")
fitted = filter(p -> !isnothing(p.distribution), results)
d = fitted[1].distribution
x = range(0, quantile(d, 0.99); length=200)

fig = Figure(size=(600, 350))
ax = Axis(fig[1, 1];
          xlabel="Days from onset", ylabel="Density",
          title="Ebola onset-to-death delay ($(typeof(d).name.name))")
band!(ax, x, zeros(length(x)), pdf.(d, x); color=(:steelblue, 0.3))
lines!(ax, x, pdf.(d, x); linewidth=2.5, color=:steelblue)
vlines!(ax, [mean(d)]; color=:red, linestyle=:dash, linewidth=1.5,
        label="Mean = $(round(mean(d); digits=1))d")
axislegend(ax; position=:rt)
fig
```

## Exercise: Explore the database

Take a few minutes to explore the `EpiParameters.jl` database:

1. Choose a disease of interest (e.g. influenza, measles, MERS)
2. How many delay distributions are available for that disease?
3. How many different types of probability distribution (gamma, log-normal, etc.) are used for a given delay?
4. Do you recognise the papers? Are there studies that should be included but aren't?

```@example epiparams
# Example: exploring Ebola
ebola_all = epiparameter(disease="Ebola")
println("Total Ebola entries: $(length(ebola_all))")
```

```@example epiparams
# Ebola incubation periods specifically
ebola_ip = epiparameter(disease="Ebola", epi_name="incubation period")
ebola_ip_fitted = filter(p -> !isnothing(p.distribution), ebola_ip)
println("Fitted Ebola incubation period entries: $(length(ebola_ip_fitted))")
for ep in ebola_ip_fitted
    println("  $(typeof(ep.distribution).name.name) — $(ep.citation["author"]) ($(ep.citation["year"]))")
end
```

## Challenge: Exploring Ebola parameters

1. How many Ebola parameter entries are in the database?
2. Find the serial interval for Ebola — what is the mean?
3. If COVID-19 and SARS have similar reproduction numbers, which would be harder to control given their serial intervals? Why?

!!! hint
    A disease with a shorter serial interval will produce new generations of infection more quickly, requiring faster response to control spread.

!!! details "Solution"
    ```@example epiparams
    ebola_si = epiparameter(disease="Ebola", epi_name="serial interval")
    ebola_si_fitted = filter(p -> !isnothing(p.distribution), ebola_si)
    if !isempty(ebola_si_fitted)
        println("Ebola serial interval mean: $(round(mean(ebola_si_fitted[1].distribution); digits=1)) days")
    end
    ```
    COVID-19 would be harder to control than SARS given similar R values, because its shorter serial interval means new generations of infection arise more quickly, requiring faster response times for contact tracing and isolation.

## Key points

- **EpiParameters.jl** provides access to a curated database of epidemiological
  delay distributions from the literature
- Use `epiparameter(disease=..., epi_name=...)` to query by disease and parameter type
- Entries with fitted parameters return `Distributions.jl` objects that can be
  plotted and used directly in analysis pipelines
- Use `list_diseases()` and `list_parameters()` to explore what's available
- Always check the citation and metadata (sample size, region, inference method)
  when selecting a parameter for your analysis
- Use the `Distributions.jl` interface (`mean`, `std`, `quantile`, `cdf`, `pdf`) to work with any retrieved distribution
- Reuse known estimates for diseases in the early stages of an outbreak when local data is not yet available


---

*Adapted from the [Epiverse-TRACE tutorials](https://epiverse-trace.github.io/tutorials/), © Epiverse-TRACE contributors, licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).*
