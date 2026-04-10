const _DATABASE = Ref{Vector{EpiParam}}()
const _DATABASE_LOADED = Ref(false)

"""Load the bundled parameter database (lazy, cached)."""
function _get_database()
    if !_DATABASE_LOADED[]
        path = joinpath(@__DIR__, "..", "data", "parameters.json")
        raw = JSON3.read(read(path, String))
        _DATABASE[] = [_parse_entry(entry) for entry in raw]
        _DATABASE_LOADED[] = true
    end
    return _DATABASE[]
end

"""Parse a single JSON entry into an EpiParam."""
function _parse_entry(entry)
    disease = String(entry[:disease])
    pathogen = _get_nullable_string(entry, :pathogen)
    epi_name = String(entry[:epi_name])

    pd = entry[:probability_distribution]
    dist, offset = _parse_distribution(pd)

    summary_stats = _to_dict(get(entry, :summary_statistics, nothing))
    citation = _to_dict(get(entry, :citation, nothing))
    metadata = _to_dict(get(entry, :metadata, nothing))
    method_assessment = _to_dict(get(entry, :method_assessment, nothing))
    notes = _get_nullable_string(entry, :notes)

    EpiParam(disease, pathogen, epi_name, dist, offset,
             summary_stats, citation, metadata, method_assessment, notes)
end

"""Parse the probability_distribution field into a Distributions.jl object."""
function _parse_distribution(pd)
    dist_name = get(pd, :prob_distribution, nothing)
    offset = Float64(get(pd, :offset, 0))
    params = get(pd, :parameters, nothing)

    if isnothing(dist_name) || isnothing(params)
        return (nothing, offset)
    end

    # Check if parameters are empty
    param_str = string(params)
    if param_str == "{}" || param_str == "{ }"
        return (nothing, offset)
    end

    dist = _build_distribution(String(dist_name), params)
    return (dist, offset)
end

"""Build a Distributions.jl distribution from R-style name and parameters."""
function _build_distribution(name::String, params)
    try
        if name == "gamma"
            shape = Float64(params[:shape])
            scale = Float64(params[:scale])
            return Gamma(shape, scale)
        elseif name == "lnorm"
            meanlog = get(params, :meanlog, nothing)
            sdlog = get(params, :sdlog, nothing)
            if !isnothing(meanlog) && !isnothing(sdlog)
                return LogNormal(Float64(meanlog), Float64(sdlog))
            end
            return nothing  # only precision reported, can't reconstruct
        elseif name == "weibull"
            shape = Float64(params[:shape])
            scale = Float64(params[:scale])
            return Weibull(shape, scale)
        elseif name == "nbinom"
            μ = Float64(params[:mean])
            k = Float64(params[:dispersion])
            # Distributions.jl NegativeBinomial(r, p) where r=dispersion, p=r/(r+μ)
            p = k / (k + μ)
            return NegativeBinomial(k, p)
        elseif name == "norm"
            μ = Float64(params[:mean])
            σ = Float64(params[:sd])
            return Normal(μ, σ)
        elseif name == "geom"
            p = Float64(params[:prob])
            return Geometric(p)
        elseif name == "pois"
            λ = Float64(params[:mean])
            return Poisson(λ)
        else
            return nothing
        end
    catch
        return nothing
    end
end

"""Convert a JSON3 object to a plain Dict{String, Any}."""
function _to_dict(obj)
    isnothing(obj) && return Dict{String, Any}()
    d = Dict{String, Any}()
    for (k, v) in pairs(obj)
        d[String(k)] = _convert_value(v)
    end
    return d
end

function _convert_value(v)
    if v isa JSON3.Object
        return _to_dict(v)
    elseif v isa JSON3.Array
        return Any[_convert_value(x) for x in v]
    elseif v isa AbstractString
        return String(v)
    else
        return v
    end
end

function _get_nullable_string(entry, key)
    v = get(entry, key, nothing)
    isnothing(v) && return nothing
    return String(v)
end
