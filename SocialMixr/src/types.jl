"""
    ContactSurvey

A social contact survey containing participant and contact data.

# Fields
- `participants::DataFrame` — one row per participant; must have `:part_id`
  column and typically `:part_age_exact`, `:country`, `:dayofweek`
- `contacts::DataFrame` — one row per reported contact; must have `:part_id`
  and typically `:cnt_age_exact`, `:cnt_age_est_min`, `:cnt_age_est_max`
- `reference::Dict{String,String}` — metadata (title, author, doi, year)
"""
struct ContactSurvey
    participants::DataFrame
    contacts::DataFrame
    reference::Dict{String,String}

    function ContactSurvey(
        participants::DataFrame,
        contacts::DataFrame,
        reference::Dict{String,String} = Dict{String,String}(),
    )
        :part_id in propertynames(participants) || throw(
            ArgumentError("participants must have a :part_id column"),
        )
        :part_id in propertynames(contacts) || throw(
            ArgumentError("contacts must have a :part_id column"),
        )
        new(participants, contacts, reference)
    end
end

function Base.show(io::IO, s::ContactSurvey)
    np = nrow(s.participants)
    nc = nrow(s.contacts)
    print(io, "ContactSurvey($np participants, $nc contacts)")
    if haskey(s.reference, "title")
        print(io, "\n  ", s.reference["title"])
    end
end

"""
    copy_survey(survey::ContactSurvey)

Return a deep copy of a `ContactSurvey`, copying both DataFrames.
"""
function copy_survey(survey::ContactSurvey)
    ContactSurvey(copy(survey.participants), copy(survey.contacts), copy(survey.reference))
end
