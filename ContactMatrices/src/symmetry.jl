"""
    make_symmetric(cm::ContactMatrix{T, 2}, population::AbstractVector{<:Real}) where T

Return a new `ContactMatrix` with symmetrised contact rates.

Symmetry is defined as `c_ij * N_i == c_ji * N_j`.  The symmetrised matrix is
computed as:

    c_sym[i,j] = (N_i * c[i,j] + N_j * c[j,i]) / (2 * N_i)

Only 2D (single-grouping) matrices are supported, matching the R package.

# Arguments
- `cm` — a 2D `ContactMatrix`
- `population` — population sizes for each group, in the same order as the
  grouping labels

# Examples
```julia
cm = ContactMatrix([4.0 1.0; 2.0 3.0], ["young", "old"])
pop = [1000.0, 500.0]
cm_sym = make_symmetric(cm, pop)
```
"""
function make_symmetric(cm::ContactMatrix{T, 2}, population::AbstractVector{<:Real}) where T
    n = size(cm, 1)
    length(population) == n || throw(ArgumentError(
        "population length $(length(population)) does not match matrix size $n"
    ))

    pop = Vector{Float64}(population)
    c = Matrix{Float64}(cm.data)

    popx = pop .* c
    sym = (popx .+ popx') ./ 2.0 ./ pop

    ContactMatrix{Float64, 2}(sym, cm.groupings, cm.setting)
end
