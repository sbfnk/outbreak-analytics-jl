"""
    epiparameter(; disease=nothing, epi_name=nothing, author=nothing) → Vector{EpiParam}

Query the epidemiological parameter database.

Filters are combined with AND logic. String matching is case-insensitive and
uses substring matching (e.g. `disease="covid"` matches "COVID-19").

# Examples
```julia
# All incubation periods for COVID-19
epiparameter(disease="COVID-19", epi_name="incubation period")

# All Ebola parameters
epiparameter(disease="Ebola")

# Parameters by a specific author
epiparameter(author="Lessler")
```
"""
function epiparameter(; disease=nothing, epi_name=nothing, author=nothing)
    db = _get_database()
    results = filter(db) do p
        _matches(p, disease, epi_name, author)
    end
    return results
end

function _matches(p::EpiParam, disease, epi_name, author)
    if !isnothing(disease)
        _icontains(p.disease, disease) || return false
    end
    if !isnothing(epi_name)
        _icontains(p.epi_name, epi_name) || return false
    end
    if !isnothing(author)
        _author_matches(p.citation, author) || return false
    end
    return true
end

"""Case-insensitive substring match."""
_icontains(haystack::String, needle) = occursin(lowercase(string(needle)), lowercase(haystack))

function _author_matches(citation::Dict, author)
    authors = get(citation, "author", Any[])
    for a in authors
        family = something(get(a, "family", nothing), "")
        given = something(get(a, "given", nothing), "")
        if _icontains(family, author) || _icontains(given, author)
            return true
        end
    end
    return false
end

"""
    list_diseases() → Vector{String}

List all unique disease names in the database.
"""
function list_diseases()
    db = _get_database()
    sort(unique(p.disease for p in db))
end

"""
    list_parameters() → Vector{String}

List all unique parameter type names in the database.
"""
function list_parameters()
    db = _get_database()
    sort(unique(p.epi_name for p in db))
end
