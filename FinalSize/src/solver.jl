"""
    final_size(R0::Real) → Float64

Solve the homogeneous final size equation `φ = 1 - exp(-R₀ φ)`.

Returns the proportion of the population infected at the end of an SIR
epidemic with basic reproduction number `R0` and homogeneous mixing.
"""
function final_size(R0::Real)
    R0 < 0 && throw(DomainError(R0, "R₀ must be non-negative"))
    R0 ≤ 1 && return 0.0

    # Fixed-point iteration: φ_{n+1} = 1 - exp(-R₀ φ_n)
    φ = 1 - 1 / R0  # good starting point
    for _ in 1:100
        φ_new = 1 - exp(-R0 * φ)
        abs(φ_new - φ) < 1e-12 && return φ_new
        φ = φ_new
    end
    return φ
end

"""
    final_size(R0::Real, cm::ContactMatrix;
               demography::AbstractVector,
               susceptibility = ones(size(Matrix(cm), 1), 1),
               p_susceptibility = ones(size(Matrix(cm), 1), 1)) → DataFrame

Solve the heterogeneous final size equation using a contact matrix.

Returns a `DataFrame` with columns `:group`, `:susc_group`, `:susceptibility`,
and `:p_infected`.

# Arguments
- `R0`: basic reproduction number (scales the contact matrix)
- `cm`: a `ContactMatrix` (must be 2D / single grouping)
- `demography`: population proportions per group (normalised internally)
- `susceptibility`: `n × k` matrix of relative susceptibilities per group and
  susceptibility class (default: all 1)
- `p_susceptibility`: `n × k` matrix of proportions in each susceptibility
  class per group (rows must sum to 1; default: all 1)
"""
function final_size(R0::Real, cm::ContactMatrix;
                    demography::AbstractVector,
                    susceptibility = ones(size(Matrix(cm), 1), 1),
                    p_susceptibility = ones(size(Matrix(cm), 1), 1))
    C = Matrix(cm)
    labels = groupings(cm).labels[1]
    _final_size_heterogeneous(R0, C, labels, demography,
                              susceptibility, p_susceptibility)
end

"""
    final_size(R0::Real, contact_matrix::AbstractMatrix;
               demography::AbstractVector,
               susceptibility = ones(size(contact_matrix, 1), 1),
               p_susceptibility = ones(size(contact_matrix, 1), 1)) → DataFrame

Solve the heterogeneous final size equation using a raw contact matrix.

See the `ContactMatrix` method for details on arguments and return value.
Group labels in the output are `"1"`, `"2"`, etc.
"""
function final_size(R0::Real, contact_matrix::AbstractMatrix;
                    demography::AbstractVector,
                    susceptibility = ones(size(contact_matrix, 1), 1),
                    p_susceptibility = ones(size(contact_matrix, 1), 1))
    C = Matrix(contact_matrix)
    labels = string.(1:size(C, 1))
    _final_size_heterogeneous(R0, C, labels, demography,
                              susceptibility, p_susceptibility)
end

"""
Internal solver for the heterogeneous final size equation.

Iterates `π_ik = 1 - exp(-s_ik Σ_j R_ij Σ_l p_jl π_jl)` to convergence,
where `R_ij = R₀ × C_ij_normalised × d_j`.
"""
function _final_size_heterogeneous(
    R0::Real,
    C::Matrix{<:Real},
    labels::Vector{String},
    demography::AbstractVector,
    susceptibility::AbstractMatrix,
    p_susceptibility::AbstractMatrix,
)
    R0 < 0 && throw(DomainError(R0, "R₀ must be non-negative"))
    n = size(C, 1)
    size(C, 2) == n ||
        throw(DimensionMismatch("contact matrix must be square"))
    length(demography) == n ||
        throw(DimensionMismatch("demography length must match matrix dimensions"))
    size(susceptibility, 1) == n ||
        throw(DimensionMismatch("susceptibility rows must match number of groups"))
    size(p_susceptibility) == size(susceptibility) ||
        throw(DimensionMismatch("p_susceptibility must have same size as susceptibility"))

    # Normalise demography to proportions
    d = Float64.(demography) ./ sum(demography)
    k = size(susceptibility, 2)
    s = Float64.(susceptibility)
    p = Float64.(p_susceptibility)

    # Normalise contact matrix by its spectral radius so R0 scales it
    ρ = maximum(real.(eigvals(C .* d')))
    if ρ < 1e-15 || R0 ≤ 0
        # Zero contacts or R₀ — no epidemic
        return _build_result(labels, k, s, p, zeros(n, k))
    end
    R = R0 .* C ./ ρ

    # R_ij includes demography: effective next-gen component
    # The equation is: π_ik = 1 - exp(-s_ik Σ_j R_ij Σ_l p_jl π_jl)
    # where R_ij = R0 * C_ij / ρ * d_j
    R_scaled = R .* d'

    # Iterate to convergence
    π = fill(0.5, n, k)
    for _ in 1:200
        # Σ_l p_jl π_jl for each group j
        weighted_pi = vec(sum(p .* π; dims=2))  # length n

        # Σ_j R_ij_scaled * weighted_pi_j for each group i
        force = R_scaled * weighted_pi  # length n

        π_new = 1 .- exp.(.-s .* force)  # n × k

        if maximum(abs.(π_new .- π)) < 1e-12
            return _build_result(labels, k, s, p, π_new)
        end
        π = π_new
    end

    return _build_result(labels, k, s, p, π)
end

"""Build the output DataFrame from converged final size proportions."""
function _build_result(labels, k, susceptibility, p_susceptibility, π)
    rows = NamedTuple{(:group, :susc_group, :susceptibility, :p_infected), Tuple{String, Int, Float64, Float64}}[]
    n = length(labels)
    for i in 1:n
        for j in 1:k
            push!(rows, (group=labels[i], susc_group=j,
                         susceptibility=susceptibility[i, j],
                         p_infected=π[i, j]))
        end
    end
    DataFrame(rows)
end
