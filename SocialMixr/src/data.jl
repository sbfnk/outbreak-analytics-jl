const DATA_DIR = joinpath(@__DIR__, "..", "data")

"""
    polymod()

Load the bundled POLYMOD survey data as a `ContactSurvey`.

The dataset contains contact data from Italy, Germany, and other European
countries, collected in 2005-2006 (Mossong et al., 2008).
"""
function polymod()
    participants = CSV.read(
        joinpath(DATA_DIR, "polymod_participants.csv"),
        DataFrame;
        missingstring = "NA",
    )
    contacts = CSV.read(
        joinpath(DATA_DIR, "polymod_contacts.csv"),
        DataFrame;
        missingstring = "NA",
    )
    ref = Dict{String,String}(
        "title" => "Social Contacts and Mixing Patterns Relevant to the Spread of Infectious Diseases",
        "author" => "Mossong et al.",
        "doi" => "10.1371/journal.pmed.0050074",
        "year" => "2008",
    )
    ContactSurvey(participants, contacts, ref)
end

"""
    polymod_population(; countries = nothing)

Load bundled WPP population data for POLYMOD countries.

Returns a `DataFrame` with columns `:lower_age_limit` and `:population`,
optionally filtered by country.
"""
function polymod_population(; countries::Union{Nothing,Vector{String}} = nothing)
    pop = CSV.read(
        joinpath(DATA_DIR, "wpp_polymod_by_country.csv"),
        DataFrame;
        missingstring = "NA",
    )
    # Normalise column names: lower.age.limit → lower_age_limit
    rename!(pop, Symbol("lower.age.limit") => :lower_age_limit)
    if countries !== nothing
        pop = filter(row -> row.country in countries, pop)
        nrow(pop) > 0 || throw(
            ArgumentError("no population data found for countries: $countries"),
        )
    end
    # Aggregate across countries
    result = combine(groupby(pop, :lower_age_limit), :population => sum => :population)
    sort!(result, :lower_age_limit)
    result
end

"""
    load_population(path::AbstractString)

Load population data from a CSV file with columns `lower.age.limit` and
`population`.
"""
function load_population(path::AbstractString)
    pop = CSV.read(path, DataFrame; missingstring = "NA")
    if hasproperty(pop, Symbol("lower.age.limit"))
        rename!(pop, Symbol("lower.age.limit") => :lower_age_limit)
    end
    :lower_age_limit in propertynames(pop) || throw(
        ArgumentError("population data must have a :lower_age_limit column"),
    )
    :population in propertynames(pop) || throw(
        ArgumentError("population data must have a :population column"),
    )
    sort!(pop, :lower_age_limit)
    pop
end
