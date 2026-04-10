"""
    EpiParam

An epidemiological parameter entry from the database.

# Fields
- `disease::String` — disease name (e.g. "COVID-19")
- `pathogen::Union{String, Nothing}` — pathogen name
- `epi_name::String` — parameter type (e.g. "incubation period")
- `distribution::Union{UnivariateDistribution, Nothing}` — fitted distribution,
  or `nothing` if only summary statistics are available
- `offset::Float64` — shift applied to the distribution (0 if unshifted)
- `summary_stats::Dict{String, Any}` — summary statistics (mean, median, sd, quantiles, etc.)
- `citation::Dict{String, Any}` — citation information (author, title, year, doi)
- `metadata::Dict{String, Any}` — units, sample_size, region, etc.
- `method_assessment::Dict{String, Any}` — truncation, censoring, discretisation flags
- `notes::Union{String, Nothing}` — additional notes
"""
struct EpiParam
    disease::String
    pathogen::Union{String, Nothing}
    epi_name::String
    distribution::Union{UnivariateDistribution, Nothing}
    offset::Float64
    summary_stats::Dict{String, Any}
    citation::Dict{String, Any}
    metadata::Dict{String, Any}
    method_assessment::Dict{String, Any}
    notes::Union{String, Nothing}
end

function Base.show(io::IO, p::EpiParam)
    print(io, "EpiParam: ", p.disease, " — ", p.epi_name)
    if !isnothing(p.distribution)
        print(io, " <", p.distribution, ">")
        if p.offset != 0
            print(io, " + ", p.offset)
        end
    end
    year = get(p.citation, "year", nothing)
    author = _first_author(p.citation)
    if !isnothing(author) && !isnothing(year)
        print(io, " (", author, " et al. ", year, ")")
    end
end

function _first_author(citation::Dict)
    authors = get(citation, "author", nothing)
    isnothing(authors) && return nothing
    isempty(authors) && return nothing
    first = authors[1]
    return get(first, "family", nothing)
end
