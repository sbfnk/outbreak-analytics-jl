# Simulating Disease Transmission

The **SEIR model** divides a population into Susceptible → Exposed → Infectious →
Recovered compartments. Combined with age-structured contact matrices, it captures
realistic epidemic dynamics.

This tutorial covers:
1. Age-structured SEIR simulation with real contact data
2. Modelling **interventions** (social distancing)
3. Comparing **vaccination** strategies

Mathematical models are useful tools for generating future trajectories of disease spread. The SEIR model categorises individuals based on their infection status:

- **S** (Susceptible): not yet infected
- **E** (Exposed): infected but not yet infectious (latent period)
- **I** (Infectious): infected and able to transmit
- **R** (Recovered): no longer infectious, assumed immune

Adding age structure via contact matrices captures how different age groups interact and transmit infections differently.

!!! note "Exposed vs infected"
    In modelling, "exposed" means infected but not yet infectious — this differs from everyday usage. The distinction between E (infected, not yet infectious) and I (infected and infectious) is important for modelling the timing of transmission.

## Setting up the model

To simulate an epidemic, we need to specify:

1. **Contact matrix** — from survey data (see the [Contact Matrices and Final Size](@ref) tutorial)
2. **Population structure** — the number of individuals in each age group
3. **Model parameters** — rates governing transitions between compartments
4. **Initial conditions** — where and how the epidemic starts

### Model parameters as rates

In population-level models defined by differential equations, parameters are specified as *rates* — the inverse of the average time until an event. The key parameters are:

- **Transmission rate β**: derived from R₀ and the recovery rate
- **Infectiousness rate σ** (= 1/latent period): if the average latent period is 2 days, then σ = 1/2 = 0.5 per day
- **Recovery rate γ** (= 1/infectious period): if the average infectious period is 5 days, then γ = 1/5 = 0.2 per day

### The SEIR differential equations

For each age group $i$, the SEIR model is described by:

$$\frac{dS_i}{dt} = -\beta S_i \sum_j C_{ij} \frac{I_j}{N_j}$$

$$\frac{dE_i}{dt} = \beta S_i \sum_j C_{ij} \frac{I_j}{N_j} - \sigma E_i$$

$$\frac{dI_i}{dt} = \sigma E_i - \gamma I_i$$

$$\frac{dR_i}{dt} = \gamma I_i$$

where $\beta$ is the transmission rate (derived from R₀ and γ), $C_{ij}$ is the contact matrix, $\sigma$ is the rate of becoming infectious (inverse of latent period), and $\gamma$ is the recovery rate (inverse of infectious period). The basic reproduction number is $R_0 = \beta / \gamma$.

!!! details "Coming from R?"
    In R's `deSolve`, ODE functions return `list(c(dS, dI, dR))` and a single `ode()` call defines and solves the system. Julia's `DifferentialEquations.jl` (wrapped by `Epidemics.jl`) splits this into steps — define the problem, then solve it — and ODE functions modify their output argument in place rather than returning values.

```@example seir
using ContactMatrices
using SocialMixr
using Epidemics
using DataFrames
using CairoMakie

survey_result = polymod() |>
    s -> filter_survey(s; countries=["United Kingdom"]) |>
    s -> assign_age_groups(s; age_limits=[0, 15, 45, 65]) |>
    compute_matrix

uk_pop = polymod_population(countries=["United Kingdom"])
sym_result = symmetrise(survey_result, uk_pop)
cm = sym_result.matrix

population = pop_age(uk_pop, [0, 15, 45, 65]).population
```

```@example seir
model = SEIR(beta=0.025, sigma=1/2, gamma=1/5)
```

```@example seir
I0 = zeros(length(population))
I0[1] = 10.0

output = simulate(model, cm;
                  demography=population,
                  initial_infected=I0,
                  tspan=(0.0, 365.0))
first(output, 8)
```

We seed the epidemic with 10 infectious individuals in the youngest age group, with all other groups initially infection-free. For models described by differential equations, "running" the model means numerically solving the system of ODEs. Julia's `DifferentialEquations.jl` ecosystem (which `Epidemics.jl` wraps) provides efficient and accurate ODE solvers for this.

### Time increments

The ODE solver uses a time step (increment) to advance the solution. When parameters are specified on a daily time scale, the default step is one day. The choice of increment should be smaller than the fastest process in the model. For example:

- If parameters are on a *daily* time scale and events are reported daily, a one-day increment is appropriate
- If some events occur on a sub-daily time scale, a smaller increment may be needed for accuracy

## Epidemic curves by age group

```@example seir
groups = sort(unique(output.group))
colors = Makie.wong_colors()

fig = Figure(size=(800, 500))

ax1 = Axis(fig[1, 1]; xlabel="Day", ylabel="Infectious",
           title="Infectious individuals by age group")
for (i, g) in enumerate(groups)
    rows = output[output.group .== g, :]
    lines!(ax1, rows.time, rows.I; linewidth=2, color=colors[i], label=g)
end
axislegend(ax1; position=:rt)

ax2 = Axis(fig[1, 2]; xlabel="Day", ylabel="Proportion",
           title="Total epidemic trajectory")
total = combine(groupby(output, :time),
                :S => sum => :S, :E => sum => :E,
                :I => sum => :I, :R => sum => :R)
N_total = sum(population)
band!(ax2, total.time, zeros(nrow(total)), total.R / N_total;
      color=(:seagreen, 0.5), label="Recovered")
band!(ax2, total.time, total.R / N_total,
      (total.R .+ total.I) / N_total;
      color=(:firebrick, 0.5), label="Infectious")
band!(ax2, total.time, (total.R .+ total.I) / N_total,
      (total.R .+ total.I .+ total.E) / N_total;
      color=(:orange, 0.5), label="Exposed")
band!(ax2, total.time, (total.R .+ total.I .+ total.E) / N_total,
      ones(nrow(total));
      color=(:steelblue, 0.5), label="Susceptible")
axislegend(ax2; position=:rt)

fig
```

!!! details "Coming from R?"
    `groupby(df, :col) |> combine` in Julia's DataFrames.jl is equivalent to `df %>% group_by(col) %>% summarise()` in dplyr. The `:col` syntax creates a `Symbol`, used for column references (like R's `quo()` or unquoted column names in dplyr). The `=>` arrow defines transformations: `:I => sum => :I` means "take column I, apply sum, name the result I".

!!! details "Coming from R?"
    Julia's `df.col` accesses a DataFrame column, equivalent to R's `df$col`. Reassuringly, many DataFrames.jl functions mirror dplyr: `nrow(df)`, `select(df, :col)`, `filter(row -> ..., df)`, `rename!(df, :old => :new)`. The `!` in `rename!` means it modifies the DataFrame in place rather than returning a copy — R's `dplyr::rename()` always returns a new data frame.

### Summary statistics

The `Epidemics.jl` package provides helper functions to extract key summary statistics:

```@example seir
# Peak timing and size for each age group
peak_data = combine(groupby(output, :group)) do gdf
    idx = argmax(gdf.I)
    (peak_day = gdf.time[idx], peak_size = gdf.I[idx])
end
peak_data
```

```@example seir
# Total epidemic size (final number recovered)
final_state = output[output.time .== 365.0, :]
total_infected = sum(final_state.R)
println("Total infected: $(round(Int, total_infected)) ($(round(total_infected/sum(population)*100; digits=1))% of population)")
```

The left panel above shows the cumulative number of infectious individuals at each time. If you want to estimate the *daily burden* (new infections per day), you can difference the cumulative recovered counts.

## Modelling interventions

Interventions such as social distancing reduce the rate of contact between groups during a specified time window. The model captures this by scaling down the contact matrix during the intervention period. This is a simplified representation — in reality, interventions affect different groups differently and compliance varies over time.

Let's compare scenarios with no intervention, 30% contact reduction, and 60% contact reduction:

```@example seir
scenarios = [
    ("No intervention", Intervention[]),
    ("30% reduction (d30–90)", [Intervention(time_begin=30.0, time_end=90.0, reduction=0.3)]),
    ("60% reduction (d30–90)", [Intervention(time_begin=30.0, time_end=90.0, reduction=0.6)]),
]

scenario_outputs = map(scenarios) do (name, ivs)
    simulate(model, cm;
             demography=population,
             initial_infected=I0,
             tspan=(0.0, 365.0),
             interventions=ivs)
end
nothing # hide
```

!!! details "Coming from R?"
    The `map(scenarios) do (name, ivs) ... end` pattern is Julia's equivalent of R's `lapply(scenarios, function(x) { ... })`. The `do` keyword creates an anonymous function that is passed as the *first* argument to `map`. The `(name, ivs)` destructures each tuple element — something R would require manual unpacking for. This `do` block syntax is used throughout Julia wherever you'd use `lapply`/`purrr::map` in R.

```@example seir
intervention_colors = [:steelblue, :orange, :firebrick]
fig = Figure(size=(700, 400))
ax = Axis(fig[1, 1]; xlabel="Day", ylabel="Total infectious",
          title="Impact of interventions on epidemic peak")

for (i, (name, _)) in enumerate(scenarios)
    total_I = combine(groupby(scenario_outputs[i], :time), :I => sum => :I)
    lines!(ax, total_I.time, total_I.I; linewidth=2,
           color=intervention_colors[i], label=name)
end

vspan!(ax, 30.0, 90.0; color=(:grey, 0.1))
axislegend(ax; position=:rt)
fig
```

The stronger intervention flattens and delays the epidemic peak — the classic "flatten the curve" effect. Note, however, that interventions don't necessarily reduce total infections unless maintained long enough. The epidemic can rebound once interventions are lifted if a sufficient proportion of the population remains susceptible.

## Comparing vaccination strategies

Who should we vaccinate? Children (who drive transmission) or the elderly
(who face higher risk)?

This is a fundamental question in vaccination policy. Children typically have the highest contact rates and drive transmission across the population. The elderly face the highest risk of severe outcomes from infection. The optimal strategy depends on whether the primary goal is to reduce overall transmission (indirect/herd protection) or to directly protect the most vulnerable.

```@example seir
vacc_scenarios = [
    ("No vaccination", Vaccination[]),
    ("Vaccinate children", [Vaccination(time_begin=0.0, time_end=30.0,
                                        rate=0.03, groups=[true, false, false, false])]),
    ("Vaccinate elderly", [Vaccination(time_begin=0.0, time_end=30.0,
                                       rate=0.03, groups=[false, false, false, true])]),
]

vacc_outputs = map(vacc_scenarios) do (name, vaccs)
    simulate(model, cm;
             demography=population,
             initial_infected=I0,
             tspan=(0.0, 365.0),
             vaccination=vaccs)
end

vacc_colors = [:steelblue, :seagreen, :firebrick]

fig = Figure(size=(700, 400))
ax = Axis(fig[1, 1]; xlabel="Day", ylabel="Total infectious",
          title="Vaccination strategy comparison")

for (i, (name, _)) in enumerate(vacc_scenarios)
    total_I = combine(groupby(vacc_outputs[i], :time), :I => sum => :I)
    lines!(ax, total_I.time, total_I.I; linewidth=2,
           color=vacc_colors[i], label=name)
end

axislegend(ax; position=:rt)
fig
```

Vaccinating the high-contact group (children) has a larger impact on **overall
transmission**, while vaccinating the elderly better protects **that specific
group**.

### Accounting for uncertainty

The epidemic model is deterministic — the same parameters always produce the same trajectory. However, in reality we may not know exact parameter values. To explore the impact of parameter uncertainty, we can run the model for multiple values drawn from a distribution.

For example, suppose R₀ follows a Normal distribution with mean 1.5 and standard deviation 0.05:

```@example seir
using Distributions: Normal as Norm
using LinearAlgebra: eigvals

R0_samples = rand(Norm(1.5, 0.05), 100)
beta_samples = R0_samples .* (1/5) ./ maximum(real.(eigvals(Matrix(cm))))

uncertainty_outputs = map(beta_samples) do β
    m = SEIR(beta=β, sigma=1/2, gamma=1/5)
    out = simulate(m, cm;
                   demography=population,
                   initial_infected=I0,
                   tspan=(0.0, 365.0))
    total_I = combine(groupby(out, :time), :I => sum => :I)
    total_I
end

fig = Figure(size=(700, 400))
ax = Axis(fig[1, 1]; xlabel="Day", ylabel="Total infectious",
          title="Uncertainty in epidemic trajectory (100 R₀ samples)")

for total_I in uncertainty_outputs
    lines!(ax, total_I.time, total_I.I; color=(:steelblue, 0.1), linewidth=0.5)
end

fig
```

Each line represents a single simulation with a different R₀ value. The spread of trajectories reflects our uncertainty about the transmission rate. Deciding which parameters to include uncertainty in depends on how well-informed a parameter value is, how sensitive outputs are to parameter changes, and the purpose of the modelling exercise.

## Challenge

From the uncertainty plot above:

1. How do the time and size of the epidemic peak change as R₀ varies within its uncertainty range?
2. Based on the definition of R₀, are these changes expected? Explain briefly.

!!! hint
    A higher R₀ means each infectious individual generates more secondary cases, leading to a faster-growing epidemic with an earlier and larger peak.

## Key points

- **SEIR models** track disease progression through compartments
- **Age structure** via contact matrices gives realistic heterogeneous dynamics
- **Interventions** flatten and delay the epidemic peak
- **Vaccination** strategies involve trade-offs between direct and indirect protection
- `Epidemics.jl` wraps `DifferentialEquations.jl` with a clean teaching API


---

*Adapted from the [Epiverse-TRACE tutorials](https://epiverse-trace.github.io/tutorials/), © Epiverse-TRACE contributors, licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).*
