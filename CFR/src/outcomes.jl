"""
    estimate_outcomes(data::DataFrame, delay_density) → DataFrame

Convolve case incidence with the delay distribution CDF to estimate the
number of cases with known outcomes at each time point.

Returns a copy of `data` with an additional `:known_outcomes` column.

# Arguments
- `data`: DataFrame with columns `:date`, `:cases`, `:deaths`
- `delay_density`: a `UnivariateDistribution` representing the delay from
  onset to outcome (e.g. `Gamma(2.4, 3.33)`)
"""
function estimate_outcomes(data::DataFrame, delay_density)
    _validate_data(data)

    n = nrow(data)
    cases = data.cases
    known = zeros(Float64, n)

    # Compute the CDF values for delays 0, 1, 2, ..., n-1
    # F(t) = probability that outcome is known by delay t
    pmf = [cdf(delay_density, t) - cdf(delay_density, t - 1) for t in 0:(n - 1)]

    # Convolve: known_outcomes_t = Σ_{s≤t} cases_s × P(delay = t-s)
    for t in 1:n
        for s in 1:t
            delay = t - s
            known[t] += cases[s] * pmf[delay + 1]
        end
    end

    result = copy(data)
    result.known_outcomes = known
    return result
end

"""Validate that the input DataFrame has the required columns."""
function _validate_data(data::DataFrame)
    for col in [:date, :cases, :deaths]
        hasproperty(data, col) ||
            throw(ArgumentError("data must have a :$col column"))
    end
end
