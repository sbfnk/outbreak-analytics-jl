"""
    simulate(model::SEIR, cm::ContactMatrix;
             demography, initial_infected, tspan,
             interventions=Intervention[], vaccination=Vaccination[],
             saveat=1.0) → DataFrame

Simulate an age-structured SEIR model using a contact matrix.

Returns a `DataFrame` with columns `:time`, `:group`, `:S`, `:E`, `:I`, `:R`
(absolute population counts).

# Arguments
- `model::SEIR` — model parameters (beta, sigma, gamma)
- `cm::ContactMatrix` — contact matrix (must be 2D / single grouping)
- `demography` — population size per age group
- `initial_infected` — initial number infected per age group
- `tspan` — `(t_start, t_end)` simulation time range
- `interventions` — vector of `Intervention`s (default: none)
- `vaccination` — vector of `Vaccination`s (default: none)
- `saveat` — output time step (default: 1.0)
"""
function simulate(model::SEIR, cm::ContactMatrix;
                  demography::AbstractVector{<:Real},
                  initial_infected::AbstractVector{<:Real},
                  tspan::Tuple{<:Real, <:Real},
                  interventions::Vector{Intervention} = Intervention[],
                  vaccination::Vector{Vaccination} = Vaccination[],
                  saveat::Real = 1.0)
    C = Matrix(cm)
    labels = groupings(cm).labels[1]
    _simulate_internal(model, C, labels, demography, initial_infected,
                       tspan, interventions, vaccination, saveat)
end

"""
    simulate(model::SEIR, contact_matrix::AbstractMatrix;
             demography, initial_infected, tspan, ...) → DataFrame

Simulate using a raw contact matrix. Group labels default to `"1"`, `"2"`, etc.
"""
function simulate(model::SEIR, contact_matrix::AbstractMatrix{<:Real};
                  demography::AbstractVector{<:Real},
                  initial_infected::AbstractVector{<:Real},
                  tspan::Tuple{<:Real, <:Real},
                  interventions::Vector{Intervention} = Intervention[],
                  vaccination::Vector{Vaccination} = Vaccination[],
                  saveat::Real = 1.0)
    C = Matrix(contact_matrix)
    labels = string.(1:size(C, 1))
    _simulate_internal(model, C, labels, demography, initial_infected,
                       tspan, interventions, vaccination, saveat)
end

function _simulate_internal(model, C, labels, demography, initial_infected,
                            tspan, interventions, vaccination, saveat)
    N = Float64.(demography)
    n = length(N)

    size(C, 1) == n || throw(DimensionMismatch(
        "contact matrix size $(size(C,1)) does not match demography length $n"))
    size(C, 2) == n || throw(DimensionMismatch(
        "contact matrix must be square"))
    length(initial_infected) == n || throw(DimensionMismatch(
        "initial_infected length must match number of age groups"))

    I0 = Float64.(initial_infected)
    u0 = zeros(4n)
    u0[1:n]      .= N .- I0    # S
    u0[2n+1:3n]  .= I0         # I

    p = (beta=model.beta, sigma=model.sigma, gamma=model.gamma,
         C=C, N=N, interventions=interventions, vaccinations=vaccination, n=n)

    prob = ODEProblem(seir_ode!, u0, Float64.(tspan), p)
    sol = solve(prob, Tsit5(); saveat=Float64(saveat))

    _solution_to_dataframe(sol, n, labels)
end
