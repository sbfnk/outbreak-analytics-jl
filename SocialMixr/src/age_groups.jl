"""
    reduce_agegroups(ages, limits)

Map individual ages to their age group lower limit using the provided limits.

Ages below the minimum limit or `missing` are returned as `missing`.

# Examples
```julia
reduce_agegroups([3, 7, 12, 20], [0, 5, 15])  # → [0, 5, 5, 15]
```
"""
function reduce_agegroups(ages, limits)
    sorted_limits = sort(limits)
    map(ages) do age
        ismissing(age) && return missing
        idx = searchsortedlast(sorted_limits, age)
        idx == 0 ? missing : sorted_limits[idx]
    end
end

"""
    limits_to_agegroups(limits)

Convert age limits to age group labels using bracket notation.

# Examples
```julia
limits_to_agegroups([0, 1, 5, 15])  # → ["[0,1)", "[1,5)", "[5,15)", "15+"]
```
"""
function limits_to_agegroups(limits)
    sorted = sort(limits)
    n = length(sorted)
    labels = String[]
    for i in 1:(n - 1)
        push!(labels, "[$(sorted[i]),$(sorted[i+1]))")
    end
    push!(labels, "$(sorted[n])+")
    labels
end

"""
    assign_age_groups(survey; kwargs...) → ContactSurvey

Process age data in a survey: impute ages from ranges, handle missing values,
and assign age groups. Returns a new `ContactSurvey` with additional columns.

# Keyword Arguments
- `age_limits` — lower limits of age groups (default: inferred from data)
- `estimated_participant_age` — how to impute participant ages from ranges:
  `"mean"` (default), `"sample"`, or `"missing"`
- `estimated_contact_age` — how to impute contact ages from ranges:
  `"mean"` (default), `"sample"`, or `"missing"`
- `missing_participant_age` — what to do with missing participant ages:
  `"remove"` (default) or `"keep"`
- `missing_contact_age` — what to do with missing contact ages:
  `"remove"` (default), `"ignore"`, or `"keep"`
"""
function assign_age_groups(
    survey::ContactSurvey;
    age_limits::Union{Nothing,Vector{<:Integer}} = nothing,
    estimated_participant_age::String = "mean",
    estimated_contact_age::String = "mean",
    missing_participant_age::String = "remove",
    missing_contact_age::String = "remove",
)
    estimated_participant_age in ("mean", "sample", "missing") || throw(
        ArgumentError("estimated_participant_age must be \"mean\", \"sample\", or \"missing\""),
    )
    estimated_contact_age in ("mean", "sample", "missing") || throw(
        ArgumentError("estimated_contact_age must be \"mean\", \"sample\", or \"missing\""),
    )
    missing_participant_age in ("remove", "keep") || throw(
        ArgumentError("missing_participant_age must be \"remove\" or \"keep\""),
    )
    missing_contact_age in ("remove", "ignore", "keep") || throw(
        ArgumentError("missing_contact_age must be \"remove\", \"ignore\", or \"keep\""),
    )

    survey = copy_survey(survey)
    participants = survey.participants
    contacts = survey.contacts

    # Add part_age from part_age_exact if not present
    _add_age!(participants, :part_age, :part_age_exact)
    _add_age!(contacts, :cnt_age, :cnt_age_exact)

    # Impute ages from ranges
    _impute_ages!(participants, :part_age, :part_age_exact, :part_age_est_min, :part_age_est_max, estimated_participant_age)
    _impute_ages!(contacts, :cnt_age, :cnt_age_exact, :cnt_age_est_min, :cnt_age_est_max, estimated_contact_age)

    # Infer age limits if not provided
    if age_limits === nothing
        age_limits = _get_age_limits(participants, contacts)
    end
    age_limits = sort(age_limits)

    # Handle missing participant ages
    if missing_participant_age == "remove"
        filter!(row -> !ismissing(row.part_age) && row.part_age >= minimum(age_limits), participants)
    end

    # Drop contacts with ages below minimum age limit
    filter!(row -> ismissing(row.cnt_age) || row.cnt_age >= minimum(age_limits), contacts)

    # Handle missing contact ages
    if missing_contact_age == "remove"
        # Remove participants that have any contact with missing age
        missing_age_ids = unique(contacts[ismissing.(contacts.cnt_age), :part_id])
        if !isempty(missing_age_ids)
            filter!(row -> !(row.part_id in missing_age_ids), participants)
        end
        filter!(row -> !ismissing(row.cnt_age), contacts)
    elseif missing_contact_age == "ignore"
        filter!(row -> !ismissing(row.cnt_age), contacts)
    end
    # "keep" → leave as is

    # Compute max age for age breaks
    max_age = _max_participant_age(participants)
    age_breaks = vcat(age_limits, max(max_age, maximum(age_limits) + 1))

    # Create age group labels
    labels = limits_to_agegroups(age_limits)

    # Assign participant age groups
    participants.lower_age_limit = reduce_agegroups(participants.part_age, age_limits)
    participants.age_group = map(participants.part_age) do age
        ismissing(age) && return missing
        idx = searchsortedlast(age_breaks, age)
        (idx < 1 || idx > length(labels)) ? missing : labels[idx]
    end

    # Assign contact age groups
    # Adjust age breaks for contacts whose age may exceed participant max
    cnt_max = _max_contact_age(contacts)
    contact_breaks = copy(age_breaks)
    if cnt_max > contact_breaks[end]
        contact_breaks[end] = cnt_max
    end

    contacts.contact_age_group = map(contacts.cnt_age) do age
        ismissing(age) && return missing
        idx = searchsortedlast(contact_breaks, age)
        (idx < 1 || idx > length(labels)) ? missing : labels[idx]
    end

    # Keep only contacts whose part_id is in the (possibly filtered) participants
    valid_ids = Set(participants.part_id)
    filter!(row -> row.part_id in valid_ids, contacts)

    ContactSurvey(participants, contacts, survey.reference)
end

# --- Internal helpers ---

function _add_age!(df, age_col, exact_col)
    if hasproperty(df, exact_col)
        if !hasproperty(df, age_col)
            df[!, age_col] = Vector{Union{Missing,Int}}(
                _to_int_or_missing.(df[!, exact_col]),
            )
        else
            # Ensure column can hold Union{Missing,Int}
            _ensure_int_or_missing!(df, age_col)
            # Fill missing ages from exact column
            for i in 1:nrow(df)
                if ismissing(df[i, age_col]) && !ismissing(df[i, exact_col])
                    df[i, age_col] = _to_int_or_missing(df[i, exact_col])
                end
            end
        end
    elseif !hasproperty(df, age_col)
        df[!, age_col] = Vector{Union{Missing,Int}}(fill(missing, nrow(df)))
    end
end

function _ensure_int_or_missing!(df, col)
    T = eltype(df[!, col])
    if !(nonmissingtype(T) <: Integer)
        df[!, col] = Vector{Union{Missing,Int}}(df[!, col])
    end
end

_to_int_or_missing(x) = ismissing(x) ? missing : Int(floor(x))

function _impute_ages!(df, age_col, exact_col, min_col, max_col, method)
    method == "missing" && return
    !hasproperty(df, min_col) && return
    !hasproperty(df, max_col) && return

    # Ensure column can hold Union{Missing,Int}
    _ensure_int_or_missing!(df, age_col)

    for i in 1:nrow(df)
        # Only impute if exact age is missing (or age_col is missing)
        age_known = hasproperty(df, exact_col) && !ismissing(df[i, exact_col])
        age_known && continue

        has_age = !ismissing(df[i, age_col])
        has_age && continue

        lo = df[i, min_col]
        hi = df[i, max_col]
        (ismissing(lo) || ismissing(hi)) && continue

        if method == "mean"
            df[i, age_col] = Int(floor((lo + hi) / 2))
        elseif method == "sample"
            df[i, age_col] = rand(lo:hi)
        end
    end
end

function _get_age_limits(participants, contacts)
    ages = Int[]
    if hasproperty(participants, :part_age)
        for a in participants.part_age
            ismissing(a) || push!(ages, a)
        end
    end
    if hasproperty(contacts, :cnt_age)
        for a in contacts.cnt_age
            ismissing(a) || push!(ages, a)
        end
    end
    union([0], sort(unique(ages)))
end

function _max_participant_age(participants)
    max_age = 0
    if hasproperty(participants, :part_age_exact)
        for a in participants.part_age_exact
            ismissing(a) || (max_age = max(max_age, Int(floor(a))))
        end
    end
    if hasproperty(participants, :part_age_est_max)
        for a in participants.part_age_est_max
            ismissing(a) || (max_age = max(max_age, Int(floor(a))))
        end
    end
    if hasproperty(participants, :part_age)
        for a in participants.part_age
            ismissing(a) || (max_age = max(max_age, a))
        end
    end
    max_age + 1
end

function _max_contact_age(contacts)
    max_age = 0
    if hasproperty(contacts, :cnt_age)
        for a in contacts.cnt_age
            ismissing(a) || (max_age = max(max_age, a))
        end
    end
    max_age + 1
end
