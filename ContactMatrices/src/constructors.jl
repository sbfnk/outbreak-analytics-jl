# ---- Tidy-format constructors ------------------------------------------------

"""
    ContactMatrix(; of, with, value, fill=0, setting=:all)

Construct a `ContactMatrix` from tidy-format vectors (matching the R
`new_contactmatrix` interface with "of"/"with" terminology).

`of` and `with` can be:
- `AbstractVector{<:AbstractString}` — single unnamed grouping
- `NamedTuple` of `AbstractVector{<:AbstractString}` — multi-grouping

Entries not specified in the triplet are filled with `fill` (default `0`).

# Examples
```julia
# Single grouping
cm = ContactMatrix(
    of    = ["[0,5)", "[5,10)", "[5,10)"],
    with  = ["[0,5)", "[10,15)", "[15,20)"],
    value = [0.32, 0.46, 0.72]
)

# Multi grouping
cm = ContactMatrix(
    of    = (age = ["young", "young", "old"], sex = ["male", "female", "female"]),
    with  = (age = ["old", "old", "young"], sex = ["female", "female", "female"]),
    value = [1.0, 2.0, 2.0]
)
```
"""
function ContactMatrix(;
    of::Union{AbstractVector{<:AbstractString}, NamedTuple},
    with::Union{AbstractVector{<:AbstractString}, NamedTuple},
    value::AbstractVector{<:Real},
    fill::Real = 0,
    setting::Symbol = :all,
)
    of_vecs, with_vecs, gnames = _parse_of_with(of, with)
    _validate_tidy_inputs(of_vecs, with_vecs, value)
    grp, data = _build_array(of_vecs, with_vecs, value, gnames, fill)
    T = eltype(data)
    N = ndims(data)
    ContactMatrix{T, N}(data, grp, setting)
end

# ---- From-matrix constructor -------------------------------------------------

"""
    ContactMatrix(matrix::AbstractMatrix{<:Real}, labels::AbstractVector{<:AbstractString}; setting=:all)

Construct a 2D `ContactMatrix` from a square matrix and a vector of group
labels (rows = "of", columns = "with").

# Examples
```julia
cm = ContactMatrix(
    [0.32 0.0; 0.0 0.46],
    ["[0,5)", "[5,10)"]
)
```
"""
function ContactMatrix(
    matrix::AbstractMatrix{<:Real},
    labels::AbstractVector{<:AbstractString};
    setting::Symbol = :all,
)
    m, n = size(matrix)
    m == n || throw(ArgumentError("matrix must be square, got $(m)×$(n)"))
    m == length(labels) || throw(ArgumentError(
        "labels length $(length(labels)) does not match matrix size $(m)"
    ))
    sorted = sort(String.(labels))
    sorted == unique(sorted) || throw(ArgumentError("labels must be unique"))

    # If labels aren't already sorted, permute the matrix to match sorted order
    perm = sortperm(String.(labels))
    data = Matrix{eltype(matrix)}(matrix[perm, perm])

    grp = Groupings(sorted)
    ContactMatrix{eltype(data), 2}(data, grp, setting)
end

# ---- Internal helpers --------------------------------------------------------

function _parse_of_with(
    of::AbstractVector{<:AbstractString},
    with::AbstractVector{<:AbstractString},
)
    return [String.(of)], [String.(with)], [:group]
end

function _parse_of_with(of::NamedTuple, with::NamedTuple)
    keys(of) == keys(with) || throw(ArgumentError(
        "grouping names in `of` and `with` must match"
    ))
    length(keys(of)) > 0 || throw(ArgumentError("at least one grouping required"))
    gnames = collect(Symbol, keys(of))
    of_vecs = [String.(of[k]) for k in keys(of)]
    with_vecs = [String.(with[k]) for k in keys(with)]
    return of_vecs, with_vecs, gnames
end

function _validate_tidy_inputs(of_vecs, with_vecs, value)
    n = length(value)
    for (i, (ov, wv)) in enumerate(zip(of_vecs, with_vecs))
        length(ov) == n || throw(ArgumentError(
            "grouping dimension $i: `of` length $(length(ov)) ≠ value length $n"
        ))
        length(wv) == n || throw(ArgumentError(
            "grouping dimension $i: `with` length $(length(wv)) ≠ value length $n"
        ))
    end
end

function _build_array(of_vecs, with_vecs, value, gnames, fill)
    ngrp = length(gnames)

    # Compute sorted unique labels per grouping (union of of/with)
    all_labels = Vector{Vector{String}}(undef, ngrp)
    for i in 1:ngrp
        all_labels[i] = sort(unique(vcat(of_vecs[i], with_vecs[i])))
    end

    grp = Groupings(gnames, all_labels)

    # Build dims: [of_1, ..., of_n, with_1, ..., with_n]
    dims = Tuple(vcat(grouplengths(grp), grouplengths(grp)))

    T = promote_type(eltype(value), typeof(fill))
    data = Array{T}(undef, dims...)
    Base.fill!(data, T(fill))

    # Create label → index mappings
    label_to_idx = [Dict(l => j for (j, l) in enumerate(all_labels[i])) for i in 1:ngrp]

    # Fill values
    for k in eachindex(value)
        idx = ntuple(2ngrp) do d
            if d <= ngrp
                label_to_idx[d][of_vecs[d][k]]
            else
                label_to_idx[d - ngrp][with_vecs[d - ngrp][k]]
            end
        end
        data[idx...] = value[k]
    end

    return grp, data
end
