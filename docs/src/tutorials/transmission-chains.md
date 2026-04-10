# Simulate Transmission Chains

In the early stages of an outbreak, we often observe **transmission chains** —
clusters of linked cases. Understanding chain statistics (sizes, lengths,
extinction probability) is critical for assessing containment potential.

**EpiBranch.jl** provides a branching process framework for simulating
transmission chains, generating line lists, modelling interventions, and
computing analytical chain statistics.

## Learning objectives

- Simulate transmission chains with offspring and generation time distributions
- Generate realistic **line lists** and **contact tables**
- Model **interventions** (isolation, contact tracing, vaccination)
- Compute analytical **extinction probabilities** and **chain size distributions**

## Basic simulation

A branching process needs an offspring distribution (how many secondary cases
each case generates) and a generation time distribution (when they generate
them):

```@example chains
using EpiBranch
using Distributions
using StableRNGs
using DataFrames
using CairoMakie

# R₀ = 2.5, dispersion k = 0.16 (SARS-like superspreading)
model = BranchingProcess(
    NegBin(2.5, 0.16),
    LogNormal(1.6, 0.5)
)
```

```@example chains
rng = StableRNG(42)
state = simulate(model;
    sim_opts = SimOpts(max_cases = 200),
    rng = rng,
)
println("Cases: $(state.cumulative_cases), Extinct: $(state.extinct)")
```

## Adding clinical presentation

To model symptom onset (needed for isolation-based interventions), provide
an incubation period via `Disease()`:

```@example chains
rng = StableRNG(42)
state = simulate(model;
    attributes = Disease(incubation_period = LogNormal(1.5, 0.5)),
    sim_opts = SimOpts(max_cases = 200),
    rng = rng,
)

ind = state.individuals[1]
println("Infection: $(round(ind.infection_time; digits=1)), " *
        "Onset: $(round(onset_time(ind); digits=1))")
```

## Generating a line list

Convert simulation output to an epidemiological line list:

```@example chains
using Dates

rng = StableRNG(42)
ll = linelist(state;
    reference_date = Date(2024, 1, 1),
    delays = DelayOpts(onset_to_reporting = Exponential(3.0)),
    outcomes = OutcomeOpts(prob_hospitalisation = 0.3, prob_death = 0.1),
    demographics = DemographicOpts(age_distribution = Normal(35, 15)),
    rng = rng,
)
first(ll, 5)
```

## Contact table

```@example chains
ct = contacts(state; reference_date = Date(2024, 1, 1))
first(ct, 5)
```

## Interventions

Interventions are composable — pass them as a vector. Let's model isolation
of symptomatic cases and contact tracing:

```@example chains
iso = Isolation(delay = Exponential(2.0))
ct_intervention = ContactTracing(probability = 0.5, delay = Exponential(1.5))

rng = StableRNG(42)
state_iv = simulate(model;
    interventions = [iso, ct_intervention],
    attributes = Disease(
        incubation_period = LogNormal(1.5, 0.5),
        prob_asymptomatic = 0.1,
    ),
    sim_opts = SimOpts(max_cases = 500),
    rng = rng,
)

println("Cases: $(state_iv.cumulative_cases)")
println("Isolated: $(count(is_isolated, state_iv.individuals))")
println("Traced: $(count(is_traced, state_iv.individuals))")
```

### Ring vaccination

Vaccinate traced contacts (post-exposure prophylaxis):

```@example chains
vacc = RingVaccination(efficacy = 0.8, delay_to_immunity = 0.0)

rng = StableRNG(42)
state_vacc = simulate(model;
    interventions = [iso, ct_intervention, vacc],
    attributes = Disease(
        incubation_period = LogNormal(1.5, 0.5),
        prob_asymptomatic = 0.1,
    ),
    sim_opts = SimOpts(max_cases = 500),
    rng = rng,
)

println("Cases (with PEP): $(state_vacc.cumulative_cases)")
println("Vaccinated: $(count(is_vaccinated, state_vacc.individuals))")
```

## Batch simulation and containment

Run many simulations to estimate containment probability:

```@example chains
rng = StableRNG(42)
states = simulate_batch(model, 500;
    attributes = Disease(incubation_period = LogNormal(1.5, 0.5)),
    sim_opts = SimOpts(max_cases = 500),
    rng = rng,
)

p_contain = containment_probability(states; max_cases = 500)
println("Containment probability: $(round(p_contain; digits=2))")
```

## Chain statistics

```@example chains
cs = chain_statistics(states)
first(cs, 5)
```

```@example chains
fig = Figure(size=(600, 300))
ax = Axis(fig[1, 1]; xlabel="Chain size", ylabel="Count",
          title="Distribution of chain sizes ($(length(states)) simulations)")
hist!(ax, cs.size; bins=50, color=(:steelblue, 0.7))
fig
```

## Analytical functions

For many quantities, closed-form solutions exist without needing simulation:

```@example chains
R0 = 2.5
k = 0.16

p_ext = extinction_probability(R0, k)
println("Extinction probability: $(round(p_ext; digits=3))")
println("Epidemic probability: $(round(1 - p_ext; digits=3))")
```

```@example chains
# What proportion of transmission comes from the top 20% of cases?
prop = proportion_transmission(R0, k; prop_cases = 0.2)
println("Top 20% of cases cause $(round(prop * 100; digits=0))% of transmission")
```

This is the classic "superspreading" result — with high overdispersion
(low k), a small fraction of cases drives most transmission.

## Comparing intervention scenarios

```@example chains
scenarios = [
    ("No intervention", AbstractIntervention[]),
    ("Isolation only", [Isolation(delay = Exponential(2.0))]),
    ("Isolation + CT", [Isolation(delay = Exponential(2.0)),
                        ContactTracing(probability = 0.5, delay = Exponential(1.5))]),
    ("Isolation + CT + PEP", [Isolation(delay = Exponential(2.0)),
                              ContactTracing(probability = 0.5, delay = Exponential(1.5)),
                              RingVaccination(efficacy = 0.8)]),
]

disease = Disease(incubation_period = LogNormal(1.5, 0.5), prob_asymptomatic = 0.1)
results = DataFrame(scenario = String[], containment = Float64[])

for (name, ivs) in scenarios
    rng = StableRNG(42)
    states = simulate_batch(model, 500;
        interventions = ivs,
        attributes = disease,
        sim_opts = SimOpts(max_cases = 500),
        rng = rng,
    )
    push!(results, (scenario = name, containment = containment_probability(states; max_cases = 500)))
end

results
```

```@example chains
fig = Figure(size=(600, 350))
ax = Axis(fig[1, 1]; xlabel="Scenario", ylabel="Containment probability",
          title="Intervention effectiveness (R₀=2.5, k=0.16)",
          xticks=(1:nrow(results), results.scenario),
          xticklabelrotation=π/6)
barplot!(ax, 1:nrow(results), results.containment; color=:steelblue)
fig
```

## Key points

- **EpiBranch.jl** simulates stochastic branching processes with offspring
  distributions, generation times, and interventions
- **Line lists** and **contact tables** are generated from simulation output,
  matching the format of real surveillance data
- **Interventions** (isolation, contact tracing, ring vaccination) are composable
  and can be scheduled with population-level triggers
- **Analytical functions** provide closed-form extinction probabilities, chain
  size distributions, and superspreading metrics without simulation
- Use `NegBin(R, k)` for the epidemiological parameterisation of offspring
  (mean R, dispersion k)
