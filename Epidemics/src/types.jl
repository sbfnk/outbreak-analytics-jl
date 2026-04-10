"""
    SEIR(; beta, sigma, gamma)

Parameters for an SEIR compartmental model.

# Fields
- `beta::Float64` — transmission rate (per contact per unit time)
- `sigma::Float64` — rate of progression from E to I (1 / latent period)
- `gamma::Float64` — recovery rate (1 / infectious period)
"""
struct SEIR
    beta::Float64
    sigma::Float64
    gamma::Float64

    function SEIR(; beta, sigma, gamma)
        beta >= 0 || throw(DomainError(beta, "beta must be non-negative"))
        sigma > 0 || throw(DomainError(sigma, "sigma must be positive"))
        gamma > 0 || throw(DomainError(gamma, "gamma must be positive"))
        new(Float64(beta), Float64(sigma), Float64(gamma))
    end
end

"""
    Intervention(; time_begin, time_end, reduction)

A time-limited reduction in transmission rate (e.g. social distancing).

During `[time_begin, time_end]`, the effective transmission rate is multiplied
by `(1 - reduction)`.
"""
struct Intervention
    time_begin::Float64
    time_end::Float64
    reduction::Float64

    function Intervention(; time_begin, time_end, reduction)
        0 <= reduction <= 1 || throw(DomainError(reduction, "reduction must be in [0, 1]"))
        time_begin < time_end || throw(ArgumentError("time_begin must be < time_end"))
        new(Float64(time_begin), Float64(time_end), Float64(reduction))
    end
end

"""
    Vaccination(; time_begin, time_end=Inf, rate, groups=nothing)

A vaccination campaign moving people from S to R at a constant rate.

# Fields
- `time_begin::Float64` — start time
- `time_end::Float64` — end time (`Inf` for ongoing)
- `rate::Float64` — vaccination rate (proportion of remaining S per unit time)
- `groups::Union{Nothing, Vector{Bool}}` — which age groups to target (`nothing` = all)
"""
struct Vaccination
    time_begin::Float64
    time_end::Float64
    rate::Float64
    groups::Union{Nothing, Vector{Bool}}

    function Vaccination(; time_begin, time_end=Inf, rate, groups=nothing)
        rate > 0 || throw(DomainError(rate, "rate must be positive"))
        time_begin < time_end || throw(ArgumentError("time_begin must be < time_end"))
        new(Float64(time_begin), Float64(time_end), Float64(rate), groups)
    end
end
