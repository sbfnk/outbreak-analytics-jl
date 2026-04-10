# Modelling Disease Burden

Understanding **disease burden** — how many people will be hospitalised,
admitted to ICU, or die — is essential for healthcare capacity planning.

This tutorial takes **population-level infection trajectories** from a
compartmental model (Epidemics.jl) and convolves them with **severity
distributions** from EpiParameters.jl to estimate hospitalisations and
deaths over time.

## Learning objectives

- Take infection output from Epidemics.jl
- Convolve with onset-to-hospitalisation and onset-to-death delays
- Estimate hospitalisations and deaths by age group
- Assess the impact of interventions on disease burden

## Step 1: Simulate an epidemic

```@example burden
using Epidemics
using ContactMatrices
using SocialMixr
using EpiParameters
using Distributions
using DataFrames
using CairoMakie

# Set up age-structured SEIR
cm = polymod() |>
    s -> filter_survey(s; countries=["United Kingdom"]) |>
    s -> assign_age_groups(s; age_limits=[0, 15, 45, 65]) |>
    compute_matrix |>
    r -> symmetrise(r, polymod_population(countries=["United Kingdom"])) |>
    r -> r.matrix

pop = pop_age(polymod_population(countries=["United Kingdom"]), [0, 15, 45, 65])
population = pop.population

model = SEIR(beta=0.025, sigma=1/2, gamma=1/5)
I0 = [10.0, 0.0, 0.0, 0.0]

result = simulate(model, cm;
    demography=population,
    initial_infected=I0,
    tspan=(0.0, 365.0))

nothing # hide
```

## Step 2: Extract daily new infections

The SEIR output gives compartment sizes. New infections per day are the
daily increase in cumulative infections (= decrease in S):

```@example burden
groups = sort(unique(result.group))

new_infections = DataFrame(time=Float64[], group=String[], new_infections=Float64[])
for g in groups
    rows = sort(result[result.group .== g, :], :time)
    daily_new = max.(0.0, -diff(rows.S))
    append!(new_infections, DataFrame(
        time = rows.time[2:end],
        group = fill(g, length(daily_new)),
        new_infections = daily_new,
    ))
end

first(new_infections, 5)
```

## Step 3: Convolve with severity delays

To estimate **hospitalisations**, we convolve daily infections with:
1. An **infection-to-hospitalisation probability** (age-specific)
2. An **onset-to-hospitalisation delay** distribution

```@example burden
# Age-specific severity parameters
severity = DataFrame(
    group = groups,
    p_hospitalisation = [0.01, 0.03, 0.08, 0.20],
    p_death = [0.001, 0.005, 0.02, 0.08],
)

# Delay distributions (could come from EpiParameters.jl)
onset_to_hosp = Gamma(3.0, 2.0)    # mean ~6 days
onset_to_death = Gamma(2.4, 3.33)  # mean ~8 days (Ebola-like, for illustration)

function convolve_delay(infections::Vector{Float64}, delay::Distribution;
                        max_delay::Int=30)
    n = length(infections)
    pmf = [cdf(delay, t) - cdf(delay, t-1) for t in 0:max_delay]
    output = zeros(n)
    for t in 1:n
        for s in 0:min(t-1, max_delay)
            output[t] += infections[t-s] * pmf[s+1]
        end
    end
    output
end
nothing # hide
```

```@example burden
burden = DataFrame(time=Float64[], group=String[], new_infections=Float64[],
                   hospitalisations=Float64[], deaths=Float64[])

for g in groups
    inf = new_infections[new_infections.group .== g, :]
    sev = severity[severity.group .== g, :]

    hosp = convolve_delay(inf.new_infections, onset_to_hosp) .* sev.p_hospitalisation[1]
    deaths = convolve_delay(inf.new_infections, onset_to_death) .* sev.p_death[1]

    append!(burden, DataFrame(
        time = inf.time,
        group = fill(g, nrow(inf)),
        new_infections = inf.new_infections,
        hospitalisations = hosp,
        deaths = deaths,
    ))
end

nothing # hide
```

## Step 4: Visualise burden

```@example burden
total = combine(groupby(burden, :time),
    :new_infections => sum => :infections,
    :hospitalisations => sum => :hospitalisations,
    :deaths => sum => :deaths,
)

fig = Figure(size=(700, 400))
ax = Axis(fig[1, 1]; xlabel="Day", ylabel="Daily count",
          title="Epidemic trajectory: infections, hospitalisations, deaths")
lines!(ax, total.time, total.infections; linewidth=2, color=:steelblue,
       label="Infections")
lines!(ax, total.time, total.hospitalisations; linewidth=2, color=:orange,
       label="Hospitalisations")
lines!(ax, total.time, total.deaths; linewidth=2, color=:firebrick,
       label="Deaths")
axislegend(ax; position=:rt)
fig
```

The hospitalisation and death curves are delayed and attenuated versions of
the infection curve — shifted by the delay distributions and scaled by
severity probabilities.

## Step 5: Burden by age group

```@example burden
# Cumulative burden per group
cumulative = combine(groupby(burden, :group),
    :new_infections => sum => :total_infections,
    :hospitalisations => sum => :total_hospitalisations,
    :deaths => sum => :total_deaths,
)
cumulative
```

```@example burden
fig = Figure(size=(700, 350))
ax = Axis(fig[1, 1]; xlabel="Age group", ylabel="Cumulative count",
          title="Total disease burden by age group",
          xticks=(1:length(groups), groups), xticklabelrotation=π/4)

x = 1:length(groups)
barplot!(ax, x .- 0.25, cumulative.total_infections; width=0.22,
         color=:steelblue, label="Infections")
barplot!(ax, x, cumulative.total_hospitalisations; width=0.22,
         color=:orange, label="Hospitalisations")
barplot!(ax, x .+ 0.25, cumulative.total_deaths; width=0.22,
         color=:firebrick, label="Deaths")
axislegend(ax; position=:lt)
fig
```

The elderly (65+) have fewer total infections (lower contact rates) but
disproportionate hospitalisations and deaths (higher severity).

## Step 6: Impact of interventions on burden

```@example burden
# Baseline vs intervention
iv = Intervention(time_begin=30.0, time_end=90.0, reduction=0.5)

result_iv = simulate(model, cm;
    demography=population,
    initial_infected=I0,
    tspan=(0.0, 365.0),
    interventions=[iv])

# Extract total infections for intervention scenario
total_iv = combine(groupby(result_iv, :time), :S => sum => :S)
total_iv.daily_inf = max.(0.0, vcat(0.0, -diff(total_iv.S)))

# Apply average death probability and delay
avg_p_death = sum(severity.p_death .* (population ./ sum(population)))
total_iv.deaths = convolve_delay(total_iv.daily_inf .* avg_p_death, onset_to_death)

fig = Figure(size=(700, 350))
ax = Axis(fig[1, 1]; xlabel="Day", ylabel="Daily deaths",
          title="Deaths: baseline vs 50% transmission reduction (days 30–90)")
lines!(ax, total.time, total.deaths; linewidth=2, color=:firebrick, label="No intervention")
lines!(ax, total_iv.time, total_iv.deaths; linewidth=2, color=:seagreen, label="50% reduction")
vspan!(ax, 30.0, 90.0; color=(:grey, 0.1))
axislegend(ax; position=:rt)
fig
```

## Key points

- **Disease burden** = infections convolved with severity delays and probabilities
- Use **Epidemics.jl** for population-level infection trajectories (appropriate
  for large outbreaks and capacity planning)
- **Age-specific severity** captures the disproportionate impact on older groups
- **Delay distributions** from EpiParameters.jl give realistic timing of
  hospitalisations and deaths relative to infections
- Interventions that reduce transmission also reduce burden — the reduction
  is amplified through the severity cascade
