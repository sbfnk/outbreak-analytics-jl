"""
    ContactMatrix{T, N} <: AbstractArray{T, N}

An N-dimensional contact matrix with grouping metadata and a setting label.

For `n` grouping dimensions the array has `2n` dimensions: the first `n`
correspond to the "of" groups (row-side / denominator population) and the
last `n` to the "with" groups (column-side).

# Fields
- `data::Array{T, N}` — the underlying array of contact rates
- `groupings::Groupings` — group labels and dimension names
- `setting::Symbol` — e.g. `:all`, `:home`, `:work`, `:school`, `:other`

# Examples
```julia
cm = ContactMatrix(
    of    = ["[0,5)", "[5,10)"],
    with  = ["[0,5)", "[10,15)"],
    value = [0.32, 0.46]
)
```
"""
struct ContactMatrix{T<:Real, N} <: AbstractArray{T, N}
    data::Array{T, N}
    groupings::Groupings
    setting::Symbol

    function ContactMatrix{T, N}(
        data::Array{T, N}, groupings::Groupings, setting::Symbol
    ) where {T<:Real, N}
        n = ndimgroups(groupings)
        2n == N || throw(ArgumentError(
            "expected $(2n) dimensions for $n groupings, got $N"
        ))
        expected = vcat(grouplengths(groupings), grouplengths(groupings))
        size(data) == Tuple(expected) || throw(ArgumentError(
            "data dimensions $(size(data)) do not match grouping sizes $(Tuple(expected))"
        ))
        new{T, N}(data, groupings, setting)
    end
end

# --- AbstractArray interface ---------------------------------------------------

Base.size(cm::ContactMatrix) = size(cm.data)
Base.getindex(cm::ContactMatrix, I...) = getindex(cm.data, I...)
Base.IndexStyle(::Type{<:ContactMatrix}) = IndexLinear()

"""
    setting(cm::ContactMatrix)

Return the setting symbol of the contact matrix.
"""
setting(cm::ContactMatrix) = cm.setting

"""
    groupings(cm::ContactMatrix)

Return the `Groupings` of the contact matrix.
"""
groupings(cm::ContactMatrix) = cm.groupings
