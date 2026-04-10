"""
    weigh(survey, by; target=nothing, groups=nothing, kwargs...) → ContactSurvey

Apply weights to survey participants. Weights are multiplicative: multiple
calls to `weigh` compose.

# Weighting modes

- **Direct** (`target=nothing, groups=nothing`): multiply weight by the numeric
  values in column `by`.
- **Grouped** (unnamed `target` vector + `groups`): map column values to groups,
  weight by `target[g] / n_in_group`.
- **Named** (named `target` dict): weight by `target[value] / n_with_value`.
- **Population** (`target` is a DataFrame): post-stratify against population
  data by age (requires `:part_age` column from `assign_age_groups`).

# Examples
```julia
# Day-of-week weighting (POLYMOD: 0=Sunday, 1=Monday, ..., 6=Saturday)
survey |> s -> weigh(s, "dayofweek"; target=[5, 2], groups=[1:5, [0, 6]])

# Age weighting against population
survey |> s -> weigh(s, "part_age"; target=pop_df)
```
"""
function weigh(
    survey::ContactSurvey,
    by::Union{String,Symbol};
    target = nothing,
    groups = nothing,
    kwargs...,
)
    survey = copy_survey(survey)
    participants = survey.participants
    by_sym = Symbol(by)

    by_sym in propertynames(participants) || throw(
        ArgumentError("column :$by_sym not found in participant data"),
    )

    # Initialise weight column if not present
    if !hasproperty(participants, :weight)
        participants.weight = ones(Float64, nrow(participants))
    end

    if target === nothing && groups === nothing
        _weigh_direct!(participants, by_sym)
    elseif target isa DataFrame
        _weigh_population!(participants, target; kwargs...)
    elseif target isa AbstractDict
        _weigh_named!(participants, by_sym, target)
    elseif groups !== nothing
        _weigh_grouped!(participants, by_sym, target, groups)
    else
        throw(ArgumentError(
            "cannot determine weighting method; provide `groups` with an unnamed `target`, a Dict `target`, or a DataFrame `target`",
        ))
    end

    ContactSurvey(participants, survey.contacts, survey.reference)
end

function _weigh_direct!(participants, by_sym)
    col = participants[!, by_sym]
    for i in 1:nrow(participants)
        v = col[i]
        ismissing(v) && continue
        participants.weight[i] *= v
    end
end

function _weigh_grouped!(participants, by_sym, target, groups)
    length(target) == length(groups) || throw(
        ArgumentError("target (length $(length(target))) and groups (length $(length(groups))) must have the same length"),
    )

    col = participants[!, by_sym]
    n = nrow(participants)

    # Assign group index to each participant
    group_idx = Vector{Union{Int,Missing}}(fill(missing, n))
    for i in 1:n
        v = col[i]
        ismissing(v) && continue
        for (g, grp) in enumerate(groups)
            if v in grp
                group_idx[i] = g
                break
            end
        end
    end

    # Count per group
    group_counts = zeros(Int, length(groups))
    for g in group_idx
        ismissing(g) || (group_counts[g] += 1)
    end

    n_total = n
    for i in 1:n
        g = group_idx[i]
        if ismissing(g)
            # Unmatched: use average weight
            participants.weight[i] *= sum(target) / n_total
        else
            gc = group_counts[g]
            gc > 0 || continue
            participants.weight[i] *= target[g] / gc
        end
    end
end

function _weigh_named!(participants, by_sym, target)
    col = participants[!, by_sym]
    val_counts = Dict{Any,Int}()
    for v in col
        ismissing(v) && continue
        val_counts[v] = get(val_counts, v, 0) + 1
    end

    for i in 1:nrow(participants)
        v = col[i]
        ismissing(v) && continue
        if haskey(target, v)
            participants.weight[i] *= target[v] / val_counts[v]
        end
    end
end

function _weigh_population!(participants, target; kwargs...)
    :part_age in propertynames(participants) || throw(
        ArgumentError("column :part_age not found; call assign_age_groups first"),
    )

    # Expand population to single-year ages
    pop_df = _prepare_pop_for_age_weighting(target, participants; kwargs...)

    # Count participants per age
    n_total = nrow(participants)
    age_counts = Dict{Int,Int}()
    for a in participants.part_age
        ismissing(a) && continue
        age_counts[a] = get(age_counts, a, 0) + 1
    end

    # Build population lookup
    pop_lookup = Dict{Int,Float64}()
    total_pop = sum(pop_df.population)
    for row in eachrow(pop_df)
        pop_lookup[row.lower_age_limit] = row.population
    end

    for i in 1:nrow(participants)
        a = participants.part_age[i]
        ismissing(a) && continue
        pop_prop = get(pop_lookup, a, 0.0) / total_pop
        age_prop = get(age_counts, a, 1) / n_total
        age_prop > 0 && (participants.weight[i] *= pop_prop / age_prop)
    end
end

function _prepare_pop_for_age_weighting(target, participants; kwargs...)
    # Get age range from participants
    ages = collect(skipmissing(participants.part_age))
    isempty(ages) && return DataFrame(lower_age_limit = Int[], population = Float64[])
    lo, hi = extrema(ages)
    pop_age(target, collect(lo:(hi + 1)); kwargs...)
end

"""
    normalise_weights!(participants; by=:age_group, threshold=nothing)

Post-stratification normalisation: within each group defined by `by`, rescale
weights so they sum to the number of participants in that group. Optionally
truncate weights above `threshold` and re-normalise.
"""
function normalise_weights!(participants::DataFrame; by = :age_group, threshold = nothing)
    by_sym = Symbol(by)
    !hasproperty(participants, by_sym) && return participants

    for gdf in groupby(participants, by_sym)
        n = nrow(gdf)
        ws = sum(gdf.weight)
        ws > 0 || continue
        gdf.weight .= gdf.weight ./ ws .* n
    end

    if threshold !== nothing
        for i in 1:nrow(participants)
            if participants.weight[i] > threshold
                participants.weight[i] = threshold
            end
        end
        # Re-normalise after truncation
        for gdf in groupby(participants, by_sym)
            n = nrow(gdf)
            ws = sum(gdf.weight)
            ws > 0 || continue
            gdf.weight .= gdf.weight ./ ws .* n
        end
    end

    participants
end
