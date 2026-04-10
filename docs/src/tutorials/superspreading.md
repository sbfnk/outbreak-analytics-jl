# Account for Superspreading

Not all infected individuals transmit equally. **Superspreading** — where a
small fraction of cases causes a disproportionate share of transmission — has
been documented in SARS, MERS, COVID-19, and many other diseases.

This tutorial covers:
1. Quantifying overdispersion with the **dispersion parameter k**
2. The "80/20 rule" of superspreading
3. **Fitting offspring distributions** to data
4. Impact of superspreading on **containment**

## Learning objectives

- Understand the negative binomial offspring distribution and dispersion parameter k
- Compute the proportion of transmission from top spreaders
- Fit offspring distributions to observed data
- Assess how superspreading affects intervention effectiveness

## The dispersion parameter k

The **negative binomial** offspring distribution has two parameters:
- **R** (mean): average number of secondary cases
- **k** (dispersion): how variable transmission is

Lower k = more overdispersion = more superspreading:

```@example superspreading
using EpiBranch
using Distributions
using DataFrames
using StableRNGs
using CairoMakie

fig = Figure(size=(700, 400))
ax = Axis(fig[1, 1]; xlabel="Secondary cases", ylabel="Probability",
          title="Offspring distributions with R₀ = 2.5",
          xticks=0:15)

for (k, col) in [(0.1, :firebrick), (0.5, :orange), (1.0, :steelblue), (Inf, :seagreen)]
    d = k == Inf ? Poisson(2.5) : NegBin(2.5, k)
    x = 0:15
    label = k == Inf ? "Poisson (k=∞)" : "NegBin (k=$k)"
    barplot!(ax, x, pdf.(d, x); width=0.2,
             offset=(k == 0.1 ? -0.3 : k == 0.5 ? -0.1 : k == 1.0 ? 0.1 : 0.3),
             color=(col, 0.7), label=label)
end

axislegend(ax; position=:rt)
fig
```

With low k (e.g. 0.1), most cases infect nobody while a few cause large
clusters — the hallmark of superspreading.

## The "80/20 rule"

What proportion of transmission comes from the top 20% of cases?

```@example superspreading
R0 = 2.5

results = DataFrame(k = Float64[], top_20_pct = Float64[])
for k in [0.01, 0.05, 0.1, 0.16, 0.3, 0.5, 1.0, 5.0]
    prop = proportion_transmission(R0, k; prop_cases = 0.2)
    push!(results, (k = k, top_20_pct = round(prop * 100; digits=1)))
end

results
```

```@example superspreading
fig = Figure(size=(600, 350))
ax = Axis(fig[1, 1]; xlabel="Dispersion parameter k", ylabel="% transmission from top 20%",
          title="Superspreading: concentration of transmission",
          xscale=log10)
scatterlines!(ax, results.k, results.top_20_pct; linewidth=2, markersize=8,
              color=:steelblue)
hlines!(ax, [80]; color=:red, linestyle=:dash, label="80% threshold")
hlines!(ax, [20]; color=:grey, linestyle=:dot, label="Homogeneous (no superspreading)")
axislegend(ax; position=:rb)
fig
```

For SARS (k ≈ 0.16), the top 20% of cases cause ~90% of transmission.

## Extinction probability

Superspreading has a silver lining: it increases the probability that
outbreaks die out on their own.

```@example superspreading
fig = Figure(size=(600, 350))
ax = Axis(fig[1, 1]; xlabel="R₀", ylabel="Extinction probability",
          title="Extinction probability by dispersion")

R_values = range(0.5, 5.0; length=100)

for (k, col, lab) in [(0.1, :firebrick, "k=0.1"), (0.5, :orange, "k=0.5"),
                       (1.0, :steelblue, "k=1.0")]
    p_ext = [extinction_probability(R, k) for R in R_values]
    lines!(ax, R_values, p_ext; linewidth=2, color=col, label=lab)
end

# Poisson (k → ∞)
p_ext_pois = [R <= 1 ? 1.0 : extinction_probability(Poisson(R)) for R in R_values]
lines!(ax, R_values, p_ext_pois; linewidth=2, color=:seagreen, label="Poisson")

axislegend(ax; position=:rt)
fig
```

Even with R₀ = 2.5, an outbreak with k = 0.1 has a ~75% chance of dying out
without intervention.

## Fitting offspring distributions to data

Given observed secondary case counts (e.g. from contact tracing), fit R and k:

```@example superspreading
# Simulated contact tracing data: observed secondary cases per index case
rng = StableRNG(42)
true_dist = NegBin(1.8, 0.3)
observed_offspring = rand(rng, true_dist, 50)

println("Observed offspring counts (first 20): $(observed_offspring[1:20])")
println("Mean: $(round(mean(observed_offspring); digits=2)), " *
        "Variance: $(round(var(observed_offspring); digits=2))")
```

```@example superspreading
data = OffspringCounts(observed_offspring)
fitted = fit(NegativeBinomial, data)

println("Fitted: R = $(round(mean(fitted); digits=2)), " *
        "k = $(round(fitted.r; digits=2))")
println("True:   R = 1.8, k = 0.3")
```

## Impact on containment

Superspreading affects intervention effectiveness. Let's compare containment
with and without overdispersion:

```@example superspreading
disease = Disease(incubation_period = LogNormal(1.5, 0.5), prob_asymptomatic = 0.1)
iso = Isolation(delay = Exponential(2.0))
ct = ContactTracing(probability = 0.5, delay = Exponential(1.5))

scenarios = DataFrame(k = Float64[], containment_no_iv = Float64[], containment_iv = Float64[])

for k in [0.1, 0.3, 0.5, 1.0, 5.0]
    model = BranchingProcess(NegBin(2.0, k), LogNormal(1.6, 0.5))

    rng = StableRNG(42)
    states_base = simulate_batch(model, 300;
        attributes = disease,
        sim_opts = SimOpts(max_cases = 300),
        rng = rng)

    rng = StableRNG(42)
    states_iv = simulate_batch(model, 300;
        interventions = [iso, ct],
        attributes = disease,
        sim_opts = SimOpts(max_cases = 300),
        rng = rng)

    push!(scenarios, (
        k = k,
        containment_no_iv = containment_probability(states_base; max_cases = 300),
        containment_iv = containment_probability(states_iv; max_cases = 300),
    ))
end

scenarios
```

```@example superspreading
fig = Figure(size=(600, 350))
ax = Axis(fig[1, 1]; xlabel="Dispersion k", ylabel="Containment probability",
          title="Superspreading and containment (R₀ = 2.0)",
          xscale=log10)
scatterlines!(ax, scenarios.k, scenarios.containment_no_iv; linewidth=2,
              markersize=8, color=:grey, label="No intervention")
scatterlines!(ax, scenarios.k, scenarios.containment_iv; linewidth=2,
              markersize=8, color=:steelblue, label="Isolation + contact tracing")
axislegend(ax; position=:rb)
fig
```

With high overdispersion (low k), interventions are *more effective* at
containment because most chains die out naturally — interventions only need
to catch the rare superspreading events.

## Key points

- **Superspreading** means a small fraction of cases drives most transmission
- The **dispersion parameter k** quantifies overdispersion: lower k = more
  superspreading
- Use `proportion_transmission(R, k)` to compute the "80/20 rule"
- Superspreading **increases** natural extinction probability but makes large
  outbreaks harder to control once established
- **Fit offspring distributions** to contact tracing data with
  `fit(NegativeBinomial, OffspringCounts(data))`
- `EpiBranch.jl` provides both simulation and analytical tools for
  superspreading analysis
