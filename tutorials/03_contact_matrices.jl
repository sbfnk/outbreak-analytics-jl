### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4c5f
begin
    import Pkg
    Pkg.develop(path=joinpath(@__DIR__, "..", "ContactMatrices"))
    Pkg.develop(path=joinpath(@__DIR__, "..", "SocialMixr"))
    Pkg.develop(path=joinpath(@__DIR__, "..", "FinalSize"))
    using ContactMatrices
    using SocialMixr
    using FinalSize
    using DataFrames
    using LinearAlgebra
    using CairoMakie
end

# ╔═╡ 1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c6f
md"""
# Contact Matrices and Final Size

Infectious diseases spread through **contacts** between people. Not all groups
have the same number of contacts: schoolchildren tend to have many close
interactions each day, whilst the elderly typically have far fewer. This
heterogeneity in contact patterns affects disease transmission because certain
groups are more likely to transmit infection both within their own group and to
other groups.

A **contact matrix** captures who contacts whom by subgroup: entry `C[i,j]`
gives the average number of contacts a person in group `i` has with people in
group `j`. Subgroups are most often defined by age category, but the same
framework applies to geographic areas, risk groups, or social settings (e.g.
households, schools, workplaces). Contact matrices summarise these rates of
interaction in a form that can be plugged directly into epidemiological models.

This tutorial covers:
1. Loading social contact survey data with **SocialMixr.jl**
2. Visualising and interpreting **contact matrices**
3. Computing the **final size** of an age-structured epidemic
"""

# ╔═╡ 2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6c7f
md"""
## Loading and processing survey data

Contact matrices are commonly estimated from **diary studies** in which
participants record every person they had contact with over the course of a day,
along with the age and setting of each interaction. The best-known of these is
the **POLYMOD** survey, which measured contact patterns across 8 European
countries (Mossong et al., 2008). Many other surveys — covering additional
countries, pandemic periods, and specific populations — are available from the
[Zenodo Social Contact Data community](https://zenodo.org/communities/social_contact_data).

The `SocialMixr.jl` package provides functions to download and process these
survey data into contact matrices. In the code below we load the POLYMOD survey,
filter to United Kingdom participants, define age groups, and compute a contact
matrix.

We also **symmetrise** the matrix using population data. Setting `symmetric=true`
ensures that the total number of contacts from group `i` to group `j` equals the
total from `j` to `i`. This reciprocity may not hold in the raw survey data due
to recall bias (people forget contacts), reporting bias (e.g. adults may not
report brief contacts with children), and sampling uncertainty (some age groups
are over- or under-represented in the survey).
"""

# ╔═╡ 3a4b5c6d-7e8f-9a0b-1c2d-3e4f5a6b7c8f
begin
    uk_survey = polymod() |>
        s -> filter_survey(s; countries=["United Kingdom"]) |>
        s -> assign_age_groups(s; age_limits=[0, 5, 18, 40, 65])

    result = compute_matrix(uk_survey)
    uk_pop = polymod_population(countries=["United Kingdom"])
    sym = symmetrise(result, uk_pop)
    cm = sym.matrix
end

# ╔═╡ 4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8c9f
md"""
## Visualising the contact matrix

Contact matrices have a characteristic pattern: high contact rates along the
diagonal (people contact others of similar age) and off-diagonal features
reflecting parent-child and workplace contacts.
"""

# ╔═╡ 5a6b7c8d-9e0f-1a2b-3c4d-5e6f7a8b9c0f
let
    labels = groupings(cm).labels[1]
    M = Matrix(cm)
    n = length(labels)

    fig = Figure(size=(550, 450))
    ax = Axis(fig[1, 1];
              xlabel="Age of contact", ylabel="Age of participant",
              title="UK contact matrix (POLYMOD)",
              xticks=(1:n, labels), yticks=(1:n, labels),
              xticklabelrotation=π/4,
              yreversed=true)
    hm = heatmap!(ax, M; colormap=:YlOrRd)
    Colorbar(fig[1, 2], hm; label="Mean contacts per day")

    # Add text values
    for i in 1:n, j in 1:n
        text!(ax, j, i; text=string(round(M[i,j]; digits=1)),
              align=(:center, :center), fontsize=11,
              color=M[i,j] > maximum(M)/2 ? :white : :black)
    end

    fig
end

# ╔═╡ 6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0c1f
md"""
Key patterns:
- **School-age children** (5–18) have the highest within-group contact rates
- **Working-age adults** (18–40, 40–65) show moderate assortative mixing
- Off-diagonal contacts between children and adults reflect household mixing
"""

# ╔═╡ 1b2c3d4e-5f6a-7b8c-9d0e-1f2a3b4c5d6e
md"""
## From contact matrices to models

Contact matrices can be used in many epidemiological analyses: calculating the
basic reproduction number R₀ while accounting for heterogeneous contacts,
computing the final epidemic size, assessing the impact of interventions such as
school closures or vaccination, and as an input to dynamic transmission models.

All of these applications require additional calculations beyond the raw contact
matrix. A key step is **normalising** the contact matrix. When simulating an
epidemic, we often want to ensure that R₀ is consistent with a known or
estimated value. Rather than using the raw contact numbers directly, we normalise
the matrix so that its largest eigenvalue equals 1, then scale by R₀. This
preserves the *relative* contact patterns — which groups interact with which —
while giving the correct reproduction number for the overall epidemic.
"""

# ╔═╡ 7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1c2f
md"""
## Final size of an epidemic

The **final size** equation gives the proportion of each age group ultimately
infected in an epidemic, assuming no intervention and no replenishment of
susceptibles (i.e. a closed population). With heterogeneous mixing, groups that
have more contacts — such as school-age children — tend to have higher attack
rates, whilst groups with fewer contacts are partially shielded.

Given a contact matrix and R₀, we can compute these proportions:
"""

# ╔═╡ 8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2c3f
begin
    pop = pop_age(uk_pop, [0, 5, 18, 40, 65])
    demography = pop.population ./ sum(pop.population)
    fs = final_size(1.5, cm; demography=demography)
end

# ╔═╡ 9a0b1c2d-3e4f-5a6b-7c8d-9e0f1a2b3c4f
let
    fig = Figure(size=(600, 350))
    ax = Axis(fig[1, 1];
              xlabel="Age group", ylabel="Proportion infected",
              title="Final epidemic size by age group (R₀ = 1.5)",
              xticks=(1:nrow(fs), fs.group),
              xticklabelrotation=π/4)
    barplot!(ax, 1:nrow(fs), fs.p_infected; color=:steelblue)
    hlines!(ax, [final_size(1.5)]; color=:red, linestyle=:dash, linewidth=1.5,
            label="Homogeneous ($(round(final_size(1.5)*100; digits=0))%)")
    axislegend(ax; position=:rt)
    fig
end

# ╔═╡ aa1b2c3d-4e5f-6a7b-8c9d-ae1f2a3b4c5f
md"""
School-age children have the highest attack rate (~65%) due to their high
contact rates, whilst the elderly are partially shielded (~29%).
The red dashed line shows the homogeneous mixing assumption — age structure matters.

The difference between the homogeneous and age-structured final sizes
illustrates why ignoring contact structure can lead to misleading projections.
Under homogeneous mixing, every group has the same attack rate. In reality,
age-structured mixing means some groups bear a disproportionate burden of
infection whilst others experience lower rates than the population average.
Getting this right matters for planning interventions — for example, targeting
resources at the groups most likely to be infected or most likely to transmit.
"""

# ╔═╡ ba2b3c4d-5e6f-7a8b-9c0d-be2f3a4b5c6f
md"""
## Sensitivity to R₀

In practice, R₀ is often uncertain or varies between settings and time periods.
It is therefore useful to explore how the final size in each age group changes
across a range of plausible R₀ values. This sensitivity analysis helps
understand how robust our conclusions are to uncertainty in R₀.
"""

# ╔═╡ ca3b4c5d-6e7f-8a9b-0c1d-ce3f4a5b6c7f
let
    R0_values = [1.1, 1.3, 1.5, 2.0, 2.5, 3.0]
    fig = Figure(size=(700, 400))
    ax = Axis(fig[1, 1];
              xlabel="R₀", ylabel="Proportion infected",
              title="Final size by age group across R₀ values")

    labels = nothing
    for (i, row) in enumerate(eachrow(
        final_size(first(R0_values), cm; demography=demography)))
        # Get group label for legend
    end

    groups = sort(unique(fs.group))
    colors = Makie.wong_colors()

    for (gi, group) in enumerate(groups)
        attack_rates = [
            final_size(R0, cm; demography=demography) |>
                df -> df[df.group .== group, :p_infected][1]
            for R0 in R0_values
        ]
        lines!(ax, R0_values, attack_rates; linewidth=2,
               color=colors[gi], label=group)
        scatter!(ax, R0_values, attack_rates; markersize=8, color=colors[gi])
    end

    # Homogeneous
    homo = [final_size(R0) for R0 in R0_values]
    lines!(ax, R0_values, homo; linewidth=2, color=:black, linestyle=:dash,
           label="Homogeneous")

    axislegend(ax; position=:lt)
    fig
end

# ╔═╡ da4b5c6d-7e8f-9a0b-1c2d-de4f5a6b7c8f
md"""
## Vaccination scenario

Contact matrices are particularly valuable for evaluating vaccination strategies,
because different strategies involve fundamentally different trade-offs.
Vaccinating **high-contact groups** (such as children) reduces transmission
across the entire population through indirect protection — fewer infections in
children means fewer opportunities for the virus to reach other age groups.
Vaccinating **high-risk groups** (such as the elderly) directly protects those
most likely to have severe outcomes, even if it has less impact on overall
transmission. The optimal strategy depends on the policy objective: minimising
total infections favours vaccinating high-contact groups, whilst minimising
severe outcomes or deaths may favour targeting high-risk groups.

Below, we model a scenario in which 50% of each age group is vaccinated with a
vaccine that reduces susceptibility by 80%:
"""

# ╔═╡ ea5b6c7d-8e9f-0a1b-2c3d-ee5f6a7b8c9f
begin
    n_groups = nrow(fs)
    susceptibility = repeat([1.0 0.2], n_groups, 1)
    p_susceptibility = repeat([0.5 0.5], n_groups, 1)

    fs_vacc = final_size(1.5, cm;
                         demography=demography,
                         susceptibility=susceptibility,
                         p_susceptibility=p_susceptibility)
end

# ╔═╡ fa6b7c8d-9e0f-1a2b-3c4d-fe6f7a8b9c0f
let
    # Overall attack rate per group (weighted by susceptibility class proportion)
    groups = sort(unique(fs.group))
    no_vacc = fs.p_infected

    vacc_weighted = [
        sum(fs_vacc[fs_vacc.group .== g, :p_infected] .*
            p_susceptibility[findfirst(==(g), groups), :])
        for g in groups
    ]

    fig = Figure(size=(600, 350))
    ax = Axis(fig[1, 1];
              xlabel="Age group", ylabel="Proportion infected",
              title="Impact of 50% vaccination (80% efficacy, R₀ = 1.5)",
              xticks=(1:length(groups), groups),
              xticklabelrotation=π/4)

    x = 1:length(groups)
    barplot!(ax, x .- 0.2, no_vacc; width=0.35, color=:steelblue, label="No vaccination")
    barplot!(ax, x .+ 0.2, vacc_weighted; width=0.35, color=:seagreen, label="50% vaccinated")
    axislegend(ax; position=:rt)
    fig
end

# ╔═╡ 0b7b8c9d-0e1f-2a3b-4c5d-0f7f8a9b0c1f
md"""
## Key points

- **Contact matrices** capture who-contacts-whom patterns by age group
- **POLYMOD** is the standard social contact survey; `SocialMixr.jl` provides
  easy access
- Matrices should be **symmetrised** using population data to enforce reciprocity
- **Final size** calculations reveal heterogeneous attack rates — age structure
  matters substantially
- **Heatmaps** are the natural visualisation for contact matrices
"""

# ╔═╡ Cell order:
# ╟─1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c6f
# ╠═0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4c5f
# ╟─2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6c7f
# ╠═3a4b5c6d-7e8f-9a0b-1c2d-3e4f5a6b7c8f
# ╟─4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8c9f
# ╠═5a6b7c8d-9e0f-1a2b-3c4d-5e6f7a8b9c0f
# ╟─6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0c1f
# ╟─1b2c3d4e-5f6a-7b8c-9d0e-1f2a3b4c5d6e
# ╟─7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1c2f
# ╠═8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2c3f
# ╠═9a0b1c2d-3e4f-5a6b-7c8d-9e0f1a2b3c4f
# ╟─aa1b2c3d-4e5f-6a7b-8c9d-ae1f2a3b4c5f
# ╟─ba2b3c4d-5e6f-7a8b-9c0d-be2f3a4b5c6f
# ╠═ca3b4c5d-6e7f-8a9b-0c1d-ce3f4a5b6c7f
# ╟─da4b5c6d-7e8f-9a0b-1c2d-de4f5a6b7c8f
# ╠═ea5b6c7d-8e9f-0a1b-2c3d-ee5f6a7b8c9f
# ╠═fa6b7c8d-9e0f-1a2b-3c4d-fe6f7a8b9c0f
# ╟─0b7b8c9d-0e1f-2a3b-4c5d-0f7f8a9b0c1f
