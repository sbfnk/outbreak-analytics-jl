# ---- Arithmetic --------------------------------------------------------------

"""
    +(cm1::ContactMatrix, cm2::ContactMatrix)

Element-wise addition of two contact matrices. Groupings must match.
The resulting setting is `:combined`.
"""
function Base.:+(cm1::ContactMatrix{T1, N}, cm2::ContactMatrix{T2, N}) where {T1, T2, N}
    cm1.groupings == cm2.groupings || throw(ArgumentError(
        "cannot add contact matrices with different groupings"
    ))
    T = promote_type(T1, T2)
    data = cm1.data .+ cm2.data
    ContactMatrix{T, N}(Array{T}(data), cm1.groupings, :combined)
end

"""
    *(scalar::Real, cm::ContactMatrix)
    *(cm::ContactMatrix, scalar::Real)

Scale contact rates by a scalar.
"""
function Base.:*(a::Real, cm::ContactMatrix{T, N}) where {T, N}
    S = promote_type(typeof(a), T)
    ContactMatrix{S, N}(Array{S}(a .* cm.data), cm.groupings, cm.setting)
end

function Base.:*(cm::ContactMatrix, a::Real)
    a * cm
end

# ---- Matrix extraction -------------------------------------------------------

"""
    Matrix(cm::ContactMatrix{T, 2}) where T

Extract the underlying 2D matrix. Errors if the contact matrix has more than
one grouping dimension.
"""
function Base.Matrix(cm::ContactMatrix{T, 2}) where T
    copy(cm.data)
end

function Base.Matrix(cm::ContactMatrix)
    throw(ArgumentError(
        "cannot convert a $(ndimgroups(cm.groupings))-grouping ContactMatrix to Matrix; " *
        "only single-grouping (2D) matrices are supported"
    ))
end

# ---- Group reduction ---------------------------------------------------------

"""
    reduce_groups(cm::ContactMatrix{T, 2}, mapping::Dict{String, String};
                  population::Union{Nothing, AbstractVector{<:Real}} = nothing) where T

Aggregate a 2D `ContactMatrix` to coarser groups defined by `mapping`.

`mapping` maps each current label to a new (coarser) label. If `population` is
provided, the reduction is population-weighted; otherwise it is a simple sum.

# Examples
```julia
cm = ContactMatrix([1.0 2.0; 3.0 4.0], ["[0,5)", "[5,10)"])
coarse = reduce_groups(cm, Dict("[0,5)" => "[0,10)", "[5,10)" => "[0,10)"))
```
"""
function reduce_groups(
    cm::ContactMatrix{T, 2},
    mapping::Dict{String, String};
    population::Union{Nothing, AbstractVector{<:Real}} = nothing,
) where T
    labels = cm.groupings.labels[1]
    n = length(labels)

    # Validate mapping covers all labels
    for l in labels
        haskey(mapping, l) || throw(ArgumentError("mapping missing label: \"$l\""))
    end

    new_labels = sort(unique(values(mapping)))
    new_n = length(new_labels)
    new_idx = Dict(l => findfirst(==(l), new_labels) for l in new_labels)

    S = promote_type(T, Float64)
    new_data = zeros(S, new_n, new_n)

    if isnothing(population)
        # Simple sum
        for i in 1:n, j in 1:n
            ni = new_idx[mapping[labels[i]]]
            nj = new_idx[mapping[labels[j]]]
            new_data[ni, nj] += cm.data[i, j]
        end
    else
        length(population) == n || throw(ArgumentError(
            "population length $(length(population)) does not match matrix size $n"
        ))
        pop = Vector{Float64}(population)
        # Population-weighted: total contacts then divide by new group population
        new_pop = zeros(Float64, new_n)
        for i in 1:n
            new_pop[new_idx[mapping[labels[i]]]] += pop[i]
        end
        for i in 1:n, j in 1:n
            ni = new_idx[mapping[labels[i]]]
            nj = new_idx[mapping[labels[j]]]
            new_data[ni, nj] += cm.data[i, j] * pop[i]
        end
        for i in 1:new_n
            new_pop[i] > 0 && (new_data[i, :] ./= new_pop[i])
        end
    end

    grp = Groupings(new_labels)
    ContactMatrix{S, 2}(new_data, grp, cm.setting)
end
