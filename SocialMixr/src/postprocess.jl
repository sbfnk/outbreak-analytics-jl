"""
    pop_age(pop, age_limits; pop_age_column=:lower_age_limit, pop_column=:population)

Adjust population data to have age groups matching the given `age_limits`.
Linearly interpolates within bands if needed, then aggregates by summing.

# Arguments
- `pop` — a `DataFrame` with population counts by age group
- `age_limits` — lower age limits to map to
- `pop_age_column` — column name for lower age limits (default `:lower_age_limit`)
- `pop_column` — column name for population (default `:population`)

# Examples
```julia
pop_age(uk_pop, [0, 18, 65])
```
"""
function pop_age(
    pop::DataFrame,
    age_limits::AbstractVector{<:Integer};
    pop_age_column::Symbol = :lower_age_limit,
    pop_column::Symbol = :population,
)
    pop = copy(pop)
    sort!(pop, pop_age_column)

    age_limits = sort(age_limits)
    max_pop_age = maximum(pop[!, pop_age_column])

    # Check if interpolation is needed
    existing_ages = pop[!, pop_age_column]
    missing_ages = setdiff(
        [a for a in age_limits if a <= max_pop_age],
        existing_ages,
    )

    if !isempty(missing_ages)
        # Linearly interpolate within bands
        pop = _interpolate_pop(pop, age_limits, pop_age_column, pop_column)
    end

    # Filter to ages >= minimum requested
    pop = pop[pop[!, pop_age_column] .>= minimum(age_limits), :]

    # Map to requested age groups and aggregate
    pop[!, :_group] = reduce_agegroups(pop[!, pop_age_column], age_limits)
    # Remove rows with missing group
    pop = pop[.!ismissing.(pop._group), :]
    result = combine(groupby(pop, :_group), pop_column => sum => pop_column)
    rename!(result, :_group => pop_age_column)
    sort!(result, pop_age_column)

    result
end

function _interpolate_pop(pop, age_limits, pop_age_column, pop_column)
    max_pop_age = maximum(pop[!, pop_age_column])

    # Build expanded DataFrame with original band info
    ages = pop[!, pop_age_column]
    pops = pop[!, pop_column]
    n = length(ages)

    # Create upper age limits for original bands
    upper = vcat(ages[2:end], [max_pop_age + 1])

    # All ages we need (from age_limits that are within range, plus originals)
    target_ages = sort(unique(vcat(
        [a for a in age_limits if a <= max_pop_age],
        ages,
    )))

    new_ages = Int[]
    new_pops = Float64[]

    for a in target_ages
        # Find which original band this age falls in
        band_idx = searchsortedlast(ages, a)
        band_idx < 1 && continue

        band_lo = ages[band_idx]
        band_hi = upper[band_idx]
        band_pop = pops[band_idx]
        band_width = band_hi - band_lo

        # Find the upper limit of the sub-band starting at `a`
        next_targets = [t for t in target_ages if t > a && t <= band_hi]
        sub_hi = isempty(next_targets) ? band_hi : minimum(next_targets)
        sub_width = sub_hi - a

        # Linearly interpolate population for this sub-band
        if band_width > 0
            sub_pop = round(band_pop * sub_width / band_width)
        else
            sub_pop = Float64(band_pop)
        end

        push!(new_ages, a)
        push!(new_pops, sub_pop)
    end

    DataFrame(pop_age_column => new_ages, pop_column => new_pops)
end

"""
    symmetrise(result, survey_pop; kwargs...) → NamedTuple

Symmetrise a contact matrix so that `c_ij * N_i ≈ c_ji * N_j`.

# Arguments
- `result` — a `NamedTuple` as returned by `compute_matrix`
- `survey_pop` — a `DataFrame` with `:lower_age_limit` and `:population`
"""
function symmetrise(result::NamedTuple, survey_pop::DataFrame; kwargs...)
    cm = result.matrix
    population = _population_for_matrix(cm, survey_pop; kwargs...)
    cm_sym = make_symmetric(cm, population)

    (matrix = cm_sym, participants = result.participants)
end

"""
    split_matrix(result, survey_pop; kwargs...) → NamedTuple

Decompose a contact matrix into components:
- `mean_contacts` — mean contacts across the population
- `normalisation` — spectral radius / mean contacts
- `contacts` — age-specific relative contact rates
- `matrix` — assortativity matrix

# Model
``m_{ij} = c \\cdot q \\cdot d_i \\cdot a_{ij} \\cdot n_j``
"""
function split_matrix(result::NamedTuple, survey_pop::DataFrame; kwargs...)
    cm = result.matrix
    population = _population_for_matrix(cm, survey_pop; kwargs...)

    m = Float64.(Matrix(cm))

    num_contacts = vec(sum(m; dims = 2))
    age_proportions = population ./ sum(population)

    mean_contacts = sum(population .* num_contacts) / sum(population)

    # Spectral radius
    eigvals = eigen(m).values
    spectral_radius = maximum(real.(eigvals))

    normalisation = spectral_radius / mean_contacts

    # Assortativity matrix: a_ij = (1/d_i) * m_ij * (1/n_j)
    assort = similar(m)
    for i in axes(m, 1)
        for j in axes(m, 2)
            d = num_contacts[i]
            n = age_proportions[j]
            assort[i, j] = (d > 0 && n > 0) ? m[i, j] / d / n : 0.0
        end
    end

    # Normalise contacts by spectral radius
    contacts_normalised = num_contacts ./ spectral_radius

    labels = groupings(cm).labels[1]
    assort_cm = ContactMatrix(assort, labels)

    (
        matrix = assort_cm,
        participants = result.participants,
        mean_contacts = mean_contacts,
        normalisation = normalisation,
        contacts = contacts_normalised,
    )
end

"""
    per_capita(result, survey_pop; kwargs...) → NamedTuple

Divide the contact matrix by population sizes to get per-capita contact rates.

``c_{ij} = m_{ij} / N_j``
"""
function per_capita(result::NamedTuple, survey_pop::DataFrame; kwargs...)
    cm = result.matrix
    population = _population_for_matrix(cm, survey_pop; kwargs...)

    m = Float64.(Matrix(cm))

    # Divide each column by its population
    for j in axes(m, 2)
        population[j] > 0 && (m[:, j] ./= population[j])
    end

    labels = groupings(cm).labels[1]
    pc_cm = ContactMatrix(m, labels)

    (matrix = pc_cm, participants = result.participants)
end

"""
    _extract_age_limits(cm::ContactMatrix)

Extract age limits from a ContactMatrix's grouping labels, sorted numerically.
"""
function _extract_age_limits(cm::ContactMatrix)
    labels = groupings(cm).labels[1]
    limits = Int[]
    for label in labels
        m = match(r"^\[?(\d+)", label)
        m !== nothing && push!(limits, parse(Int, m.captures[1]))
    end
    sort(limits)
end

"""
    _population_for_matrix(cm, survey_pop; kwargs...)

Return a population vector ordered to match the ContactMatrix's internal label
ordering (which is alphabetical). This is critical because `Matrix(cm)` and
`make_symmetric` use the alphabetical label order.
"""
function _population_for_matrix(cm::ContactMatrix, survey_pop::DataFrame; kwargs...)
    age_limits = _extract_age_limits(cm)
    resolved_pop = pop_age(survey_pop, age_limits; kwargs...)

    # Build a lookup from numeric age limit → population
    pop_lookup = Dict{Int,Float64}()
    for row in eachrow(resolved_pop)
        pop_lookup[row.lower_age_limit] = Float64(row.population)
    end

    # Return population in the same order as the ContactMatrix's labels
    labels = groupings(cm).labels[1]
    population = Float64[]
    for label in labels
        m = match(r"^\[?(\d+)", label)
        if m !== nothing
            age = parse(Int, m.captures[1])
            push!(population, get(pop_lookup, age, 0.0))
        end
    end

    population
end
