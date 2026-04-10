# Choosing an Appropriate Model

Different epidemiological questions call for different modelling approaches.
This tutorial compares the three main frameworks available in the Julia
ecosystem and discusses when to use each.

## Model frameworks

| Framework | Package | Best for |
|---|---|---|
| **Branching process** | EpiBranch.jl | Early outbreak, individual chains, stochastic extinction, superspreading |
| **Compartmental (SEIR)** | Epidemics.jl | Population-level dynamics, interventions at scale, vaccination strategies |
| **Renewal equation** | EpiNow2.jl | Rt estimation from surveillance data, forecasting, real-time analysis |

## When to use a branching process

Branching processes model **individual transmission chains**. They are most
appropriate when:

- The outbreak is **small** (tens to hundreds of cases)
- You care about **stochastic extinction** — will the outbreak die out?
- **Superspreading** matters (negative binomial offspring)
- You want to generate **line lists** or **contact networks**
- You need to model **individual-level interventions** (isolation, contact
  tracing, ring vaccination)

```@example choosing_bp
using EpiBranch
using Distributions
using StableRNGs

# Early Ebola-like outbreak: will it take off?
model = BranchingProcess(NegBin(1.5, 0.3), LogNormal(2.0, 0.8))

rng = StableRNG(42)
states = simulate_batch(model, 1000;
    sim_opts = SimOpts(max_cases = 500),
    rng = rng)

p_extinct = count(s -> s.extinct, states) / length(states)
println("Extinction probability: $(round(p_extinct; digits=2))")
println("Median outbreak size (conditional on non-extinction): " *
        "$(median([s.cumulative_cases for s in states if !s.extinct]))")
```

**Limitation**: branching processes assume an unlimited susceptible population.
They break down once a significant fraction of the population is infected.

## When to use a compartmental model

Compartmental models (SEIR) track **population-level flows** between disease
states. Use them when:

- The outbreak is **large** (thousands+ of cases)
- You need **age structure** via contact matrices
- You want to compare **population-level interventions** (lockdowns,
  vaccination campaigns)
- **Depletion of susceptibles** matters
- You need the **final size** of an epidemic

```@example choosing_seir
using Epidemics
using ContactMatrices
using SocialMixr
using DataFrames

# UK flu-like scenario
cm = polymod() |>
    s -> filter_survey(s; countries=["United Kingdom"]) |>
    s -> assign_age_groups(s; age_limits=[0, 18, 65]) |>
    compute_matrix |>
    r -> symmetrise(r, polymod_population(countries=["United Kingdom"])) |>
    r -> r.matrix

pop = pop_age(polymod_population(countries=["United Kingdom"]), [0, 18, 65]).population

result = simulate(
    SEIR(beta=0.025, sigma=1/2, gamma=1/5), cm;
    demography=pop,
    initial_infected=[10.0, 0.0, 0.0],
    tspan=(0.0, 365.0))

final = result[result.time .== 365.0, :]
println("Final attack rates:")
for row in eachrow(final)
    rate = round(row.R / (row.S + row.E + row.I + row.R) * 100; digits=1)
    println("  $(row.group): $(rate)%")
end
```

**Limitation**: deterministic — no stochastic extinction, no individual-level
detail. Best for large populations where the law of large numbers applies.

## When to use a renewal equation

The renewal equation (EpiNow2.jl) is a **statistical model** fitted to observed
data. Use it when:

- You have **surveillance data** (daily case counts) and want to estimate Rt
- You need **real-time situational awareness** — is the epidemic growing?
- You want a **short-term forecast**
- You need to account for **reporting delays** and **observation noise**

```julia
using EpiNow2

cases = example_confirmed()
result = epinow(cases;
    generation_time = gt_opts(discretise(LogNormal(1.2, 0.5); max=14)),
    inference = inference_opts(samples=1000, chains=3))

# Is Rt above 1?
latest_rt = result.estimates.rt[end, :]
println("Latest Rt: $(round(latest_rt.median; digits=2)) " *
        "(90% CrI: $(round(latest_rt.lower_90; digits=2))–$(round(latest_rt.upper_90; digits=2)))")
```

**Limitation**: requires observed data — cannot simulate hypothetical scenarios
like compartmental models can.

## Decision tree

```
Is this a hypothetical scenario or real data?
├─ Real data → EpiNow2.jl (Rt estimation, forecasting)
└─ Hypothetical scenario
   ├─ Small outbreak / individual level?
   │  └─ EpiBranch.jl (branching process)
   └─ Large population / age structure?
      └─ Epidemics.jl (compartmental SEIR)
```

## Combining approaches

In practice, models are often combined:

- **EpiNow2** estimates Rt from data → feeds into **Epidemics.jl** scenarios
- **EpiBranch** assesses early containment → if containment fails,
  **Epidemics.jl** models the large-scale epidemic
- **EpiParameters.jl** provides delay distributions to all three frameworks
- **ContactMatrices.jl** from **SocialMixr.jl** informs both **Epidemics.jl**
  and **FinalSize.jl**

## Key points

- **Branching processes** (EpiBranch.jl): stochastic, individual-level, early
  outbreak, superspreading, containment
- **Compartmental models** (Epidemics.jl): deterministic, population-level,
  age-structured, interventions at scale
- **Renewal equation** (EpiNow2.jl): statistical, data-driven, real-time Rt,
  forecasting
- Choose based on your **question**, not your data — different questions about
  the same outbreak may need different models
- The Julia ecosystem packages share types (`ContactMatrix`, `Distributions.jl`
  objects) enabling smooth hand-offs between frameworks
