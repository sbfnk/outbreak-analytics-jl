"""
    filter_survey(survey; countries=nothing, filter=nothing) → ContactSurvey

Filter a survey by country and/or contact-level columns.

# Arguments
- `countries` — vector of country names to keep (filters participants)
- `filter` — `Dict` mapping column names (Symbol or String) to values; contacts
  not matching all filters are removed

# Examples
```julia
survey |> s -> filter_survey(s; countries=["United Kingdom"])
survey |> s -> filter_survey(s; filter=Dict(:cnt_school => 1))
```
"""
function filter_survey(
    survey::ContactSurvey;
    countries::Union{Nothing,Vector{String}} = nothing,
    filter::Union{Nothing,Dict} = nothing,
)
    survey = copy_survey(survey)
    participants = survey.participants
    contacts = survey.contacts

    # Filter by country
    if countries !== nothing && hasproperty(participants, :country)
        DataFrames.filter!(row -> row.country in countries, participants)
        nrow(participants) > 0 || throw(
            ArgumentError("no participants left after selecting countries: $countries"),
        )
        valid_ids = Set(participants.part_id)
        DataFrames.filter!(row -> row.part_id in valid_ids, contacts)
    end

    # Apply contact-level filters
    if filter !== nothing
        for (col, val) in filter
            col_sym = Symbol(col)
            if hasproperty(contacts, col_sym)
                DataFrames.filter!(row -> !ismissing(row[col_sym]) && row[col_sym] == val, contacts)
            end
        end
    end

    ContactSurvey(participants, contacts, survey.reference)
end

"""
    compute_matrix(survey; counts=false, weight_threshold=nothing) → NamedTuple

Compute a contact matrix from a survey that has been processed by
`assign_age_groups` (and optionally `weigh`).

Returns a `NamedTuple` with:
- `matrix` — a `ContactMatrix` from ContactMatrices.jl
- `participants` — a `DataFrame` with participant counts per age group

# Arguments
- `counts` — if `true`, return total contact counts instead of means (default `false`)
- `weight_threshold` — if set, truncate standardised weights above this value
  and re-normalise (default `nothing`)
"""
function compute_matrix(
    survey::ContactSurvey;
    counts::Bool = false,
    weight_threshold::Union{Nothing,Real} = nothing,
)
    survey = copy_survey(survey)
    participants = survey.participants
    contacts = survey.contacts

    :age_group in propertynames(participants) || throw(
        ArgumentError("column :age_group not found; call assign_age_groups first"),
    )
    :contact_age_group in propertynames(contacts) || throw(
        ArgumentError("column :contact_age_group not found; call assign_age_groups first"),
    )

    # Initialise weight if not present
    if !hasproperty(participants, :weight)
        participants.weight = ones(Float64, nrow(participants))
    end

    # Post-stratification normalisation with optional threshold
    normalise_weights!(participants; by = :age_group, threshold = weight_threshold)

    # Merge participants and contacts
    merged = innerjoin(contacts, participants; on = :part_id, makeunique = true)

    # Get ordered age group labels
    age_groups = _ordered_age_groups(participants)

    n_groups = length(age_groups)
    matrix = zeros(Float64, n_groups, n_groups)

    # Build index mapping
    group_idx = Dict(g => i for (i, g) in enumerate(age_groups))

    # Accumulate weighted contacts
    for row in eachrow(merged)
        pg = get(group_idx, row.age_group, 0)
        cg = get(group_idx, row.contact_age_group, 0)
        (pg == 0 || cg == 0) && continue
        w = hasproperty(merged, :weight) ? row.weight : 1.0
        ismissing(w) && continue
        matrix[pg, cg] += w
    end

    if !counts
        # Normalise by sum of weights per participant age group
        norm_vector = zeros(Float64, n_groups)
        for row in eachrow(participants)
            g = get(group_idx, row.age_group, 0)
            g == 0 && continue
            w = row.weight
            ismissing(w) && continue
            norm_vector[g] += w
        end
        for i in 1:n_groups
            if norm_vector[i] > 0
                matrix[i, :] ./= norm_vector[i]
            else
                matrix[i, :] .= NaN
            end
        end
    end

    # Build ContactMatrix
    cm = ContactMatrix(matrix, age_groups)

    # Participant population summary
    part_pop = combine(groupby(participants, :age_group), nrow => :participants)
    part_pop.proportion = part_pop.participants ./ sum(part_pop.participants)
    sort!(part_pop, :age_group; by = _age_group_lower_limit)

    (matrix = cm, participants = part_pop)
end

"""
    _ordered_age_groups(participants)

Extract ordered age group labels from participant data. Labels are sorted by
their numeric lower limit.
"""
function _ordered_age_groups(participants)
    groups = unique(skipmissing(participants.age_group))
    # Sort by numeric lower limit
    sort(collect(groups); by = _age_group_lower_limit)
end

function _age_group_lower_limit(label::AbstractString)
    m = match(r"^\[?(\d+)", label)
    m === nothing ? 0 : parse(Int, m.captures[1])
end
