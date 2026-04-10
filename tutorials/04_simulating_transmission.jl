### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4c60
begin
    import Pkg
    Pkg.develop(path=joinpath(@__DIR__, "..", "ContactMatrices"))
    Pkg.develop(path=joinpath(@__DIR__, "..", "SocialMixr"))
    Pkg.develop(path=joinpath(@__DIR__, "..", "Epidemics"))
    using ContactMatrices
    using SocialMixr
    using Epidemics
    using DataFrames
    using CairoMakie
end

# ╔═╡ 1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c60
md"""
# Simulating Disease Transmission

Mathematical models are useful tools for generating future trajectories of disease spread. By encoding our understanding of how infections transmit between individuals, we can explore "what if" scenarios and evaluate potential interventions before they are deployed.

The **SEIR model** is one of the most widely used frameworks in infectious disease epidemiology. It categorises individuals by their infection status into four compartments:

- **S**usceptible — not yet infected and at risk of infection
- **E**xposed — infected but not yet infectious (in a latent period)
- **I**nfectious — infected and capable of transmitting the pathogen
- **R**ecovered — no longer infectious (assumed immune)

Individuals progress through these compartments in sequence: S → E → I → R. Adding **age structure** via contact matrices captures how different age groups interact and transmit infections differently — for example, school-age children tend to have many more daily contacts than older adults, which has important implications for epidemic dynamics and control strategies.

This tutorial covers:
1. Age-structured SEIR simulation with real contact data
2. Modelling **interventions** (social distancing)
3. Comparing **vaccination** strategies
"""

# ╔═╡ 2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6c70
md"""
## Setting up the model

An SEIR model is governed by several key parameters:

- **Transmission rate β**: controls how quickly the disease spreads. This is often derived from the basic reproduction number R₀ and the recovery rate.
- **Contact matrix C**: encodes the frequency of contacts between age groups, typically estimated from survey data such as POLYMOD.
- **Infectiousness rate σ** (sigma): the rate at which exposed individuals become infectious. The average latent period (time spent in the E compartment) is 1/σ. For example, if the pre-infectious period is 2 days on average, then σ = 1/2 = 0.5 per day.
- **Recovery rate γ** (gamma): the rate at which infectious individuals recover. The average infectious period is 1/γ. For example, if people are infectious for 5 days on average, then γ = 1/5 = 0.2 per day.

Note that these **rates are the inverse of average durations** — a common convention in compartmental modelling that can be counterintuitive at first.

A note on terminology: "Exposed" in epidemiological modelling means infected but not yet infectious, which differs from everyday usage where "exposed" might simply mean having been in contact with an infected person. The distinction between the E compartment (infected, not yet infectious) and the I compartment (infected and infectious) is important because it determines the timing of onward transmission.

We begin by constructing the contact matrix from POLYMOD survey data and obtaining the UK population structure.
"""

# ╔═╡ 3a4b5c6d-7e8f-9a0b-1c2d-3e4f5a6b7c80
begin
    survey_result = polymod() |>
        s -> filter_survey(s; countries=["United Kingdom"]) |>
        s -> assign_age_groups(s; age_limits=[0, 15, 45, 65]) |>
        compute_matrix

    uk_pop = polymod_population(countries=["United Kingdom"])
    sym_result = symmetrise(survey_result, uk_pop)
    cm = sym_result.matrix

    population = pop_age(uk_pop, [0, 15, 45, 65]).population
end

# ╔═╡ 4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8c90
model = SEIR(beta=0.025, sigma=1/2, gamma=1/5)

# ╔═╡ a1b2c3d4-e5f6-7a8b-9c0d-a1b2c3d4e5f6
md"""
We now need to set the **initial conditions** — the state of the population at the start of the simulation. We seed the epidemic with 10 infectious individuals in the youngest age group (0–14), with all other age groups initially entirely susceptible and infection-free.

The choice of initial conditions affects the early dynamics of the epidemic (e.g., which age groups are affected first) but generally has less impact on the final epidemic size and overall dynamics, which are primarily determined by the model parameters and the contact structure.
"""

# ╔═╡ 5a6b7c8d-9e0f-1a2b-3c4d-5e6f7a8b9ca0
begin
    I0 = zeros(length(population))
    I0[1] = 10.0

    output = simulate(model, cm;
                      demography=population,
                      initial_infected=I0,
                      tspan=(0.0, 365.0))
end

# ╔═╡ b2c3d4e5-f6a7-8b9c-0d1e-b2c3d4e5f6a7
md"""
**What does "running the model" mean?** The SEIR model is described by a system of ordinary differential equations (ODEs) that specify how the number of individuals in each compartment changes over time. "Running the model" means numerically solving this system of ODEs to find the trajectory of compartment sizes from the initial conditions forward in time.

Julia's `DifferentialEquations.jl` ecosystem (which `Epidemics.jl` wraps) is particularly well-suited for this task, providing efficient and accurate ODE solvers. The output is a DataFrame with the number of individuals in each compartment at each time point, broken down by age group.
"""

# ╔═╡ 6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0cb0
md"""
## Epidemic curves by age group
"""

# ╔═╡ 7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1cc0
let
    groups = sort(unique(output.group))
    colors = Makie.wong_colors()

    fig = Figure(size=(800, 500))

    # Infectious
    ax1 = Axis(fig[1, 1]; xlabel="Day", ylabel="Infectious",
               title="Infectious individuals by age group")
    for (i, g) in enumerate(groups)
        rows = output[output.group .== g, :]
        lines!(ax1, rows.time, rows.I; linewidth=2, color=colors[i], label=g)
    end
    axislegend(ax1; position=:rt)

    # All compartments (total)
    ax2 = Axis(fig[1, 2]; xlabel="Day", ylabel="Population",
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
end

# ╔═╡ 8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2cd0
md"""
## Modelling interventions

Interventions such as social distancing, school closures, or stay-at-home orders reduce the rate of contact between groups during a specified time window. In the model, we capture this by **scaling down the contact matrix** during the intervention period — for example, a 30% reduction means contacts are reduced to 70% of their baseline level.

This is a simplified representation of reality. In practice, interventions affect different groups differently (e.g., school closures primarily reduce contacts among children), compliance varies over time, and there may be behavioural adaptation. Nevertheless, this approach captures the essential mechanism: fewer contacts mean fewer opportunities for transmission.

An `Intervention` reduces transmission during a time window. Let's compare scenarios with different levels of contact reduction:
"""

# ╔═╡ 9a0b1c2d-3e4f-5a6b-7c8d-9e0f1a2b3ce0
begin
    scenarios = [
        ("No intervention", Intervention[]),
        ("30% reduction (d30–90)", [Intervention(time_begin=30.0, time_end=90.0, reduction=0.3)]),
        ("60% reduction (d30–90)", [Intervention(time_begin=30.0, time_end=90.0, reduction=0.6)]),
    ]

    scenario_outputs = map(scenarios) do (name, ivs)
        out = simulate(model, cm;
                       demography=population,
                       initial_infected=I0,
                       tspan=(0.0, 365.0),
                       interventions=ivs)
        out.scenario .= name
        out
    end
end

# ╔═╡ aa1b2c3d-4e5f-6a7b-8c9d-ae1f2a3b4cf0
let
    colors = [:steelblue, :orange, :firebrick]
    fig = Figure(size=(700, 400))
    ax = Axis(fig[1, 1]; xlabel="Day", ylabel="Total infectious",
              title="Impact of interventions on epidemic peak")

    for (i, (name, _)) in enumerate(scenarios)
        df = scenario_outputs[i]
        total_I = combine(groupby(df, :time), :I => sum => :I)
        lines!(ax, total_I.time, total_I.I; linewidth=2, color=colors[i], label=name)
    end

    # Intervention period
    vspan!(ax, 30.0, 90.0; color=(:grey, 0.1))
    text!(ax, 60, maximum(combine(groupby(scenario_outputs[1], :time),
          :I => sum => :I).I) * 0.95;
          text="Intervention\nperiod", align=(:center, :top), fontsize=10, color=:grey)

    axislegend(ax; position=:rt)
    fig
end

# ╔═╡ ba2b3c4d-5e6f-7a8b-9c0d-be2f3a4b5d00
md"""
The stronger intervention flattens and delays the peak — the classic "flatten the curve" effect. However, it is important to note that while interventions flatten and delay the peak, they do not necessarily reduce the total number of infections unless they are maintained for long enough. The trade-off between intervention strength and duration is a key consideration in public health planning.

Notice that the epidemic can **rebound** once interventions are lifted if a sufficient number of susceptible individuals remain in the population. This means that a short, sharp intervention may simply postpone the epidemic rather than prevent it. Sustaining interventions until enough of the population has been infected (or vaccinated) to achieve herd immunity is often necessary to avoid this rebound effect.
"""

# ╔═╡ ca3b4c5d-6e7f-8a9b-0c1d-ce3f4a5b6d10
md"""
## Comparing vaccination strategies

Who should we vaccinate? This is a fundamental question in vaccination policy that arises with every new vaccine. Children typically have the highest contact rates (as seen in the contact matrices from earlier tutorials), making them key drivers of transmission. The elderly, on the other hand, generally face the highest risk of severe outcomes such as hospitalisation and death.

The optimal strategy depends on whether the primary goal is to **reduce overall transmission** (providing indirect or herd protection to the whole population) or to **directly protect the most vulnerable** groups. In practice, vaccination programmes often attempt to balance both objectives.

Let's compare two strategies: vaccinating children versus vaccinating the elderly.
"""

# ╔═╡ da4b5c6d-7e8f-9a0b-1c2d-de4f5a6b7d20
begin
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
end

# ╔═╡ ea5b6c7d-8e9f-0a1b-2c3d-ee5f6a7b8d30
let
    colors = [:steelblue, :seagreen, :firebrick]
    groups = sort(unique(vacc_outputs[1].group))

    fig = Figure(size=(800, 400))

    # Total infectious
    ax1 = Axis(fig[1, 1]; xlabel="Day", ylabel="Total infectious",
               title="Total epidemic curve")
    for (i, (name, _)) in enumerate(vacc_scenarios)
        total_I = combine(groupby(vacc_outputs[i], :time), :I => sum => :I)
        lines!(ax1, total_I.time, total_I.I; linewidth=2, color=colors[i], label=name)
    end
    axislegend(ax1; position=:rt)

    # Final attack rate by group
    ax2 = Axis(fig[1, 2]; xlabel="Age group", ylabel="Proportion infected",
               title="Final attack rate by strategy",
               xticks=(1:length(groups), groups), xticklabelrotation=π/4)

    width = 0.25
    for (i, (name, _)) in enumerate(vacc_scenarios)
        final = vacc_outputs[i][vacc_outputs[i].time .== 365.0, :]
        N_grp = population[sortperm(groups)]  # sorted to match labels
        attack = final.R ./ (final.S .+ final.E .+ final.I .+ final.R)
        barplot!(ax2, (1:length(groups)) .+ (i-2)*width, attack;
                 width=width, color=colors[i], label=name)
    end
    axislegend(ax2; position=:rt)

    fig
end

# ╔═╡ fa6b7c8d-9e0f-1a2b-3c4d-fe6f7a8b9d40
md"""
**Key insight**: vaccinating children has a larger impact on **overall
transmission** (indirect protection), while vaccinating the elderly better
protects **that specific group** (direct protection). The optimal strategy
depends on the objective.
"""

# ╔═╡ c3d4e5f6-a7b8-9c0d-1e2f-c3d4e5f6a7b8
md"""
## A note on uncertainty

The model we have used here is **deterministic** — the same parameters always produce exactly the same trajectory. In reality, disease transmission involves considerable randomness (stochasticity), and many of the model parameters are uncertain. For example, the exact value of R₀ for a novel pathogen is rarely known precisely.

One approach to account for parameter uncertainty is to run the model across a range of plausible parameter values — for instance, by sampling R₀ from a probability distribution reflecting our uncertainty, running the model for each sampled value, and summarising the resulting range of trajectories. This gives a sense of how robust our conclusions are to uncertainty in the inputs.

Stochastic models, which explicitly incorporate randomness in the transmission process, provide another way to capture uncertainty — particularly important when case numbers are small (e.g., at the start of an outbreak) and chance events can have a large influence on the trajectory.
"""

# ╔═╡ 0b7b8c9d-0e1f-2a3b-4c5d-0f7f8a9b0d50
md"""
## Key points

- **SEIR models** track disease progression through compartments
- **Age structure** via contact matrices gives realistic heterogeneous dynamics
- **Interventions** flatten and delay the epidemic peak
- **Vaccination** strategies involve trade-offs between direct and indirect
  protection
- `Epidemics.jl` wraps `DifferentialEquations.jl` with a clean teaching API,
  returning tidy DataFrames for analysis and visualisation with `CairoMakie`
"""

# ╔═╡ Cell order:
# ╟─1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c60
# ╠═0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4c60
# ╟─2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6c70
# ╠═3a4b5c6d-7e8f-9a0b-1c2d-3e4f5a6b7c80
# ╠═4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8c90
# ╟─a1b2c3d4-e5f6-7a8b-9c0d-a1b2c3d4e5f6
# ╠═5a6b7c8d-9e0f-1a2b-3c4d-5e6f7a8b9ca0
# ╟─b2c3d4e5-f6a7-8b9c-0d1e-b2c3d4e5f6a7
# ╟─6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0cb0
# ╠═7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1cc0
# ╟─8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2cd0
# ╠═9a0b1c2d-3e4f-5a6b-7c8d-9e0f1a2b3ce0
# ╠═aa1b2c3d-4e5f-6a7b-8c9d-ae1f2a3b4cf0
# ╟─ba2b3c4d-5e6f-7a8b-9c0d-be2f3a4b5d00
# ╟─ca3b4c5d-6e7f-8a9b-0c1d-ce3f4a5b6d10
# ╠═da4b5c6d-7e8f-9a0b-1c2d-de4f5a6b7d20
# ╠═ea5b6c7d-8e9f-0a1b-2c3d-ee5f6a7b8d30
# ╟─fa6b7c8d-9e0f-1a2b-3c4d-fe6f7a8b9d40
# ╟─c3d4e5f6-a7b8-9c0d-1e2f-c3d4e5f6a7b8
# ╟─0b7b8c9d-0e1f-2a3b-4c5d-0f7f8a9b0d50
