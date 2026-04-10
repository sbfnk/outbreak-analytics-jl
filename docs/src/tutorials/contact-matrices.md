# Contact Matrices and Final Size

Infectious diseases spread through **contacts** between people. A **contact matrix**
captures who contacts whom by age group: entry `C[i,j]` gives the average number
of contacts a person in group `i` has with people in group `j`.

This tutorial covers:
1. Loading social contact survey data with **SocialMixr.jl**
2. Visualising and interpreting **contact matrices**
3. Computing the **final size** of an age-structured epidemic

Some groups of individuals have more contacts than others — the average schoolchild has many more daily contacts than the average elderly person. This heterogeneity in contact patterns affects disease transmission, because certain groups are more likely to transmit to others within their own group as well as to other groups.

The subgroups used in a contact matrix are often age categories, but they can also represent geographic areas, risk groups (e.g. occupational risk), or social settings (e.g. household, workplace, school).

### A note on notation

In a contact matrix, the entry `C[i,j]` at row `i` and column `j` represents the average number of contacts an individual in group `i` has with individuals in group `j`. This is calculated by dividing the total number of contacts between groups `i` and `j` by the size of group `i`.

## Loading and processing survey data

Contact matrices are commonly estimated from diary studies where participants record their daily interactions. The POLYMOD survey, one of the most widely used contact studies, measured contact patterns in 8 European countries using data on the location and duration of contacts reported by study participants.

The `SocialMixr.jl` package provides functions to load and process these survey data. Other surveys are available from the Zenodo Social Contact Data community.

When computing a contact matrix, we symmetrise it to ensure that the total number of contacts from group `i` to group `j` equals the total from `j` to `i`. Raw survey data may not satisfy this property due to recall bias (different age groups remember contacts differently), reporting bias, and sampling uncertainty.

```@example contacts
using ContactMatrices
using SocialMixr
using FinalSize
using DataFrames
using LinearAlgebra
using CairoMakie

uk_survey = polymod() |>
    s -> filter_survey(s; countries=["United Kingdom"]) |>
    s -> assign_age_groups(s; age_limits=[0, 5, 18, 40, 65])

result = compute_matrix(uk_survey)
uk_pop = polymod_population(countries=["United Kingdom"])
sym = symmetrise(result, uk_pop)
cm = sym.matrix
```

!!! details "Coming from R?"
    Julia's `|>` pipe operator works like R's `|>` (base pipe) or `%>%` (magrittr). In Julia, `x |> f` passes `x` as the argument to `f`. For functions that need additional arguments, use anonymous functions: `x |> s -> f(s; kwarg=value)`, which is analogous to R's `x |> (\(s) f(s, kwarg = value))()`.

!!! details "Coming from R?"
    In Julia, a semicolon `;` separates positional arguments from keyword arguments: `assign_age_groups(s; age_limits=[0, 5, 18, 40, 65])`. R uses commas for everything: `assign_age_groups(s, age_limits = c(0, 5, 18, 40, 65))`. You'll see this pattern throughout Julia code — `round(x; digits=2)`, `range(0, 10; length=200)`, etc. In practice, a comma often works too, but the semicolon is the conventional style.

## Visualising the contact matrix

```@example contacts
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

for i in 1:n, j in 1:n
    text!(ax, j, i; text=string(round(M[i,j]; digits=1)),
          align=(:center, :center), fontsize=11,
          color=M[i,j] > maximum(M)/2 ? :white : :black)
end

fig
```

Contact matrices have a characteristic structure reflecting social mixing patterns:

Key patterns:
- **School-age children** (5–18) have the highest within-group contact rates
- **Working-age adults** (18–40, 40–65) show moderate assortative mixing
- Off-diagonal contacts between children and adults reflect household mixing

## From contact matrices to epidemiological analyses

Contact matrices can be used in a wide range of epidemiological analyses:

- Calculating the basic reproduction number while accounting for heterogeneous contacts between age groups
- Computing the final size of an epidemic
- Assessing the impact of interventions by comparing pre- and post-intervention contact matrices
- In mathematical models of transmission, to account for group-specific contact patterns

All of these applications require additional calculations. In particular, when simulating an epidemic, we often want to ensure that the model's basic reproduction number R₀ is consistent with known values for the pathogen. Rather than using raw contact numbers, we normalise the contact matrix by scaling it so that its largest eigenvalue is 1. This preserves the relative contact patterns while making it straightforward to set R₀ to any desired value.

## Contact matrices in mathematical models

To understand how contact matrices enter epidemic models, consider the standard SIR model where individuals are Susceptible (S), Infectious (I), or Recovered (R):

$$\frac{dS}{dt} = -\beta S \frac{I}{N}, \quad \frac{dI}{dt} = \beta S \frac{I}{N} - \gamma I, \quad \frac{dR}{dt} = \gamma I$$

To add age structure, we introduce separate equations for each age group $i$ and replace the single transmission term with one that accounts for heterogeneous contacts via the contact matrix $C$:

$$\frac{dS_i}{dt} = -\beta S_i \sum_j C_{ij} \frac{I_j}{N_j}$$

$$\frac{dI_i}{dt} = \beta S_i \sum_j C_{ij} \frac{I_j}{N_j} - \gamma I_i$$

$$\frac{dR_i}{dt} = \gamma I_i$$

Susceptible individuals in age group $i$ become infected at a rate that depends on their contacts with infectious individuals across *all* age groups, weighted by the contact matrix.

### Normalising the contact matrix

When simulating an epidemic, we typically want to set a specific value of R₀. The contact matrix $C$ captures the *pattern* of contacts, but we need to scale it so the model produces the desired R₀.

The approach is to normalise the matrix so its largest eigenvalue (spectral radius) equals 1. This normalised matrix $\hat{C}$ preserves the relative contact patterns. We then define the transmission rate as $\beta = R_0 \cdot \gamma / \rho(C)$, where $\rho(C)$ is the spectral radius of the original contact matrix. This ensures that the model's basic reproduction number matches the target R₀.

```@example contacts
C = Matrix(cm)
eigenvalue = maximum(real.(eigvals(C)))
C_normalised = C ./ eigenvalue
println("Largest eigenvalue: $(round(eigenvalue; digits=2))")
println("Largest eigenvalue of normalised matrix: $(round(maximum(real.(eigvals(C_normalised))); digits=2))")
```

In `Epidemics.jl`, this normalisation happens internally when you pass the contact matrix to the model — you don't need to do it manually.

### Synthetic contact matrices

Empirical contact data from diary studies is not available for all countries. Prem et al. (2017, 2021) used POLYMOD data within a Bayesian hierarchical model to project contact matrices for 177 countries. These synthetic matrices provide estimates of age-specific contact patterns even where no survey has been conducted, and are commonly used when local empirical data is unavailable.

## Final size of an epidemic

The final size equation gives the proportion of each age group ultimately infected, assuming no intervention and no replenishment of susceptibles. With heterogeneous mixing, groups with more contacts have higher attack rates, whilst groups with fewer contacts are partially shielded.

Given a contact matrix and R₀, we compute the final size:

```@example contacts
pop = pop_age(uk_pop, [0, 5, 18, 40, 65])
demography = pop.population ./ sum(pop.population)
fs = final_size(1.5, cm; demography=demography)
```

```@example contacts
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
```

School-age children have the highest attack rate due to their high contact rates, whilst the elderly are partially shielded. The red dashed line shows the final size under the homogeneous mixing assumption — the difference demonstrates that ignoring contact structure can lead to misleading projections. In reality, age-structured mixing means some groups bear a disproportionate burden of infection.

## Sensitivity to R₀

Since R₀ is often uncertain or varies between settings, it is useful to explore how the final size changes across a range of values:

```@example contacts
R0_values = [1.1, 1.3, 1.5, 2.0, 2.5, 3.0]
groups = sort(unique(fs.group))
colors = Makie.wong_colors()

fig = Figure(size=(700, 400))
ax = Axis(fig[1, 1];
          xlabel="R₀", ylabel="Proportion infected",
          title="Final size by age group across R₀ values")

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

homo = [final_size(R0) for R0 in R0_values]
lines!(ax, R0_values, homo; linewidth=2, color=:black, linestyle=:dash,
       label="Homogeneous")

axislegend(ax; position=:lt)
fig
```

## Vaccination scenario

A key question in vaccination policy is *who* to vaccinate. Different strategies involve different trade-offs:

- **Vaccinating high-contact groups** (e.g. children) reduces transmission across the entire population through indirect (herd) protection
- **Vaccinating high-risk groups** (e.g. the elderly) directly protects those most likely to experience severe outcomes

The optimal strategy depends on whether the primary goal is to minimise total infections or to minimise severe outcomes. Here we model 50% coverage with a vaccine that reduces susceptibility by 80%:

```@example contacts
n_groups = nrow(fs)
susceptibility = repeat([1.0 0.2], n_groups, 1)
p_susceptibility = repeat([0.5 0.5], n_groups, 1)

fs_vacc = final_size(1.5, cm;
                     demography=demography,
                     susceptibility=susceptibility,
                     p_susceptibility=p_susceptibility)
```

```@example contacts
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
```

Vaccinating the high-contact group (children) has a larger impact on overall transmission, while vaccinating the elderly better protects that specific group. The optimal strategy depends on the objective.

## Contact groups beyond age

While we have focused on age-structured contact matrices, the same framework applies to any grouping of the population. For example:

- **Geographic areas**: a matrix representing contact rates between regions in a metapopulation model
- **Risk groups**: high-risk vs low-risk occupational categories
- **Social settings**: separate matrices for household, workplace, school, and community contacts, which can be combined or modified independently when modelling setting-specific interventions

The dimension of the contact matrix must match the number of groups in the model.

## Key points

- **Contact matrices** capture who-contacts-whom patterns by age group
- **POLYMOD** is the standard social contact survey; `SocialMixr.jl` provides easy access
- Matrices should be **symmetrised** using population data to enforce reciprocity
- **Final size** calculations reveal heterogeneous attack rates — age structure matters
- **Heatmaps** are the natural visualisation for contact matrices
- The **spectral radius** (largest eigenvalue) of the contact matrix determines R₀ — normalisation makes it easy to set a target R₀
- Synthetic contact matrices (Prem et al.) provide estimates for countries without survey data
- Contact matrices can represent any population grouping, not just age


---

*Adapted from the [Epiverse-TRACE tutorials](https://epiverse-trace.github.io/tutorials/), © Epiverse-TRACE contributors, licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).*
