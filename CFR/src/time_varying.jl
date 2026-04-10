"""
    cfr_time_varying(data::DataFrame;
                     delay_density=nothing,
                     burn_in=7,
                     smoothing_window=nothing,
                     confidence_level=0.95) → DataFrame

Estimate the time-varying case fatality ratio.

Returns a DataFrame with columns `:date`, `:cfr`, `:lower`, `:upper`,
`:cases`, `:deaths`, and `:known_outcomes` (if delay-adjusted).

# Arguments
- `data`: DataFrame with columns `:date`, `:cases`, `:deaths`
- `delay_density`: a `UnivariateDistribution` for the onset-to-death delay,
  or `nothing` for naive time-varying CFR
- `burn_in`: number of initial time points to discard (default 7)
- `smoothing_window`: if set, apply a rolling sum over this many days to
  cases and deaths before computing the CFR
- `confidence_level`: confidence level for intervals (default 0.95)
"""
function cfr_time_varying(data::DataFrame;
                          delay_density=nothing,
                          burn_in=7,
                          smoothing_window=nothing,
                          confidence_level=0.95)
    _validate_data(data)
    n = nrow(data)

    cases = Float64.(data.cases)
    deaths = Float64.(data.deaths)

    # Optional smoothing
    if !isnothing(smoothing_window)
        cases = _rolling_sum(cases, smoothing_window)
        deaths = _rolling_sum(deaths, smoothing_window)
    end

    if isnothing(delay_density)
        # Naive time-varying: cumulative deaths / cumulative cases
        cum_cases = cumsum(cases)
        cum_deaths = cumsum(deaths)

        cfr_vals = Float64[]
        lower_vals = Float64[]
        upper_vals = Float64[]

        for t in 1:n
            est = _binomial_cfr(cum_deaths[t], cum_cases[t], confidence_level)
            push!(cfr_vals, est.estimate)
            push!(lower_vals, est.lower)
            push!(upper_vals, est.upper)
        end

        result = DataFrame(
            date=data.date,
            cfr=cfr_vals,
            lower=lower_vals,
            upper=upper_vals,
        )
    else
        corrected = estimate_outcomes(data, delay_density)
        known = corrected.known_outcomes

        cum_cases = cumsum(cases)
        cum_deaths = cumsum(deaths)
        cum_known = cumsum(known)

        cfr_vals = Float64[]
        lower_vals = Float64[]
        upper_vals = Float64[]

        for t in 1:n
            if cum_known[t] < 1e-10
                push!(cfr_vals, NaN)
                push!(lower_vals, NaN)
                push!(upper_vals, NaN)
                continue
            end

            u_t = min(cum_known[t] / cum_cases[t], 1.0)
            adjusted = u_t * cum_cases[t]
            est = _binomial_cfr(cum_deaths[t], adjusted, confidence_level)
            push!(cfr_vals, est.estimate)
            push!(lower_vals, est.lower)
            push!(upper_vals, est.upper)
        end

        result = DataFrame(
            date=data.date,
            cfr=cfr_vals,
            lower=lower_vals,
            upper=upper_vals,
            known_outcomes=corrected.known_outcomes,
        )
    end

    # Apply burn-in
    if burn_in > 0 && burn_in < n
        result = result[(burn_in + 1):end, :]
    end

    return result
end

"""Compute a rolling sum with the given window size."""
function _rolling_sum(x::AbstractVector, window::Int)
    n = length(x)
    result = similar(x, Float64)
    for i in 1:n
        lo = max(1, i - window + 1)
        result[i] = sum(@view x[lo:i])
    end
    return result
end
