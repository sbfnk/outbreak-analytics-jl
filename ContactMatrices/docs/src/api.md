# API Reference

## Types

```@docs
ContactMatrix
Groupings
```

## Constructors

```@docs
ContactMatrix(::AbstractMatrix{<:Real}, ::AbstractVector{<:AbstractString})
```

## Functions

```@docs
setting
groupings
ndimgroups
make_symmetric
reduce_groups
```

## Arithmetic

```@docs
Base.:+(::ContactMatrix, ::ContactMatrix)
Base.:*(::Real, ::ContactMatrix)
Base.Matrix(::ContactMatrix{T, 2}) where T
```
