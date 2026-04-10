"""
    cfr_static(data::DataFrame;
               delay_density=nothing,
               confidence_level=0.95) → CfrEstimate

Estimate the static (overall) case fatality ratio.

Without `delay_density`, returns the naive CFR (total deaths / total cases).
With a delay distribution, applies the Nishiura et al. (2009) correction for
right-censoring.

# Arguments
- `data`: DataFrame with columns `:date`, `:cases`, `:deaths`
- `delay_density`: a `UnivariateDistribution` for the onset-to-death delay,
  or `nothing` for naive CFR
- `confidence_level`: confidence level for the interval (default 0.95)
"""
function cfr_static(data::DataFrame;
                    delay_density=nothing,
                    confidence_level=0.95)
    _validate_data(data)

    total_deaths = sum(data.deaths)
    total_cases = sum(data.cases)

    if isnothing(delay_density)
        # Naive CFR
        return _binomial_cfr(total_deaths, total_cases, confidence_level)
    end

    # Delay-adjusted: u_t correction
    corrected = estimate_outcomes(data, delay_density)
    known = sum(corrected.known_outcomes)

    if known < 1e-10
        return CfrEstimate(NaN, NaN, NaN)
    end

    # Adjusted denominator: proportion with known outcomes × total cases
    u = min(known / total_cases, 1.0)
    adjusted_cases = u * total_cases

    _binomial_cfr(total_deaths, adjusted_cases, confidence_level)
end

"""
Compute CFR with Clopper-Pearson exact binomial confidence interval.

Uses the Beta distribution to compute quantiles:
- Lower: Beta(deaths, n - deaths + 1) at α/2
- Upper: Beta(deaths + 1, n - deaths) at 1 - α/2
"""
function _binomial_cfr(deaths, n, confidence_level)
    if n < 1e-10
        return CfrEstimate(NaN, NaN, NaN)
    end

    p = deaths / n
    # Clamp to [0, 1] — can exceed 1 with delay adjustment if deaths > known
    p = clamp(p, 0.0, 1.0)

    α = 1 - confidence_level
    deaths_int = round(Int, min(deaths, n))
    n_int = max(round(Int, n), 1)

    if deaths_int ≤ 0
        lower = 0.0
    else
        lower = quantile(Beta(deaths_int, n_int - deaths_int + 1), α / 2)
    end

    if deaths_int ≥ n_int
        upper = 1.0
    else
        upper = quantile(Beta(deaths_int + 1, n_int - deaths_int), 1 - α / 2)
    end

    CfrEstimate(p, lower, upper)
end
