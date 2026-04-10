"""
    Groupings

Describes the grouping dimensions of a `ContactMatrix`.

Each grouping (e.g. age, sex) has a name (`Symbol`) and a vector of sorted
unique labels (`Vector{String}`).  A single unnamed grouping is stored with
name `:group`.

# Fields
- `names::Vector{Symbol}` — dimension names, e.g. `[:age]` or `[:age, :sex]`
- `labels::Vector{Vector{String}}` — sorted unique labels per dimension
"""
struct Groupings
    names::Vector{Symbol}
    labels::Vector{Vector{String}}

    function Groupings(names::Vector{Symbol}, labels::Vector{Vector{String}})
        length(names) == length(labels) ||
            throw(ArgumentError("`names` and `labels` must have the same length"))
        length(names) > 0 ||
            throw(ArgumentError("at least one grouping is required"))
        for (n, l) in zip(names, labels)
            length(l) > 0 ||
                throw(ArgumentError("labels for grouping :$n must not be empty"))
            l == sort(unique(l)) ||
                throw(ArgumentError("labels for grouping :$n must be sorted and unique"))
        end
        new(names, labels)
    end
end

"""
    Groupings(labels::AbstractVector{<:AbstractString})

Convenience constructor for a single unnamed grouping.
"""
function Groupings(labels::AbstractVector{<:AbstractString})
    sorted = sort(unique(String.(labels)))
    Groupings([:group], [sorted])
end

"""
    Groupings(; kwargs...)

Named multi-grouping via keyword arguments.

# Examples
```julia
Groupings(age = ["[0,5)", "[5,10)"], sex = ["female", "male"])
```
"""
function Groupings(; kwargs...)
    ns = Symbol[]
    ls = Vector{String}[]
    for (k, v) in kwargs
        push!(ns, k)
        push!(ls, sort(unique(String.(v))))
    end
    Groupings(ns, ls)
end

"""
    ndimgroups(g::Groupings)

Number of grouping dimensions.
"""
ndimgroups(g::Groupings) = length(g.names)

"""
    grouplength(g::Groupings, i::Int)

Number of labels in grouping dimension `i`.
"""
grouplength(g::Groupings, i::Int) = length(g.labels[i])

"""
    grouplengths(g::Groupings)

Vector of label counts for each grouping dimension.
"""
grouplengths(g::Groupings) = [length(l) for l in g.labels]

function Base.:(==)(a::Groupings, b::Groupings)
    a.names == b.names && a.labels == b.labels
end

function Base.show(io::IO, g::Groupings)
    n = ndimgroups(g)
    if n == 1
        print(io, "Groupings(", g.names[1], ": ", join(g.labels[1], ", "), ")")
    else
        parts = [string(nm, ": ", join(lb, ", ")) for (nm, lb) in zip(g.names, g.labels)]
        print(io, "Groupings(", join(parts, "; "), ")")
    end
end
