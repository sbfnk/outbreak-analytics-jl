# Create a Short-Term Forecast

Beyond estimating ``R_t``, we often need to **forecast** — predict future case
counts to inform resource planning and public health response.

EpiNow2.jl produces forecasts automatically by projecting the estimated ``R_t``
trajectory forward.

Given case data of an epidemic, we can create estimates of current and future case numbers by accounting for delays in reporting and under-reporting. To make predictions about the future, we need to make an assumption about how observations up to the present relate to what we expect to happen. The simplest approach is to assume "no change" — the reproduction number remains the same as its most recent estimate. This is the default approach in EpiNow2.

Internally, `epinow()` wraps two functions: `estimate_infections()` (estimates cases by date of infection from the observed data) and `forecast_infections()` (simulates future infections using the estimated ``R_t``).

!!! warning "Computation time"
    MCMC sampling takes several minutes per model fit.

## Basic forecast

The `forecast_opts(horizon=7)` argument specifies a 7-day forecast beyond the last observed data point. The forecast is generated as part of the MCMC sampling — for each posterior sample of ``R_t``, future infections are simulated, naturally propagating uncertainty.

```@example forecast
using EpiNow2
using EpiParameters
using Distributions
using DataFrames
using Dates
using CairoMakie

cases = example_confirmed()

gt = discretise(LogNormal(1.2, 0.5); max=14)
delay = discretise(LogNormal(1.4, 0.4); max=14)

result = epinow(cases;
    generation_time = gt_opts(gt),
    delays = delay_opts(delay),
    forecast = forecast_opts(horizon=7),
    inference = inference_opts(samples=1000, warmup=250, chains=3,
                               progress=false),
    verbose = false)
nothing # hide
```

### Visualising the forecast

The reports DataFrame includes both fitted values (for observed dates) and forecasts (for future dates). The estimates are split into three categories:

- **Estimate** (blue): uses all available data
- **Estimate based on partial data**: higher uncertainty because recent infections haven't been fully observed
- **Forecast** (red): projects into the future assuming ``R_t`` remains constant

The grey dashed line separates the observation period from the forecast.

```@example forecast
rep = result.estimates.reports
n_obs = nrow(cases)
n_total = nrow(rep)

fig = Figure(size=(750, 400))
ax = Axis(fig[1, 1]; xlabel="Date", ylabel="Cases",
          title="Case forecast (7-day horizon)")

barplot!(ax, 1:n_obs, cases.confirm; color=(:grey, 0.5), label="Observed")

band!(ax, 1:n_obs, rep.lower_90[1:n_obs], rep.upper_90[1:n_obs];
      color=(:steelblue, 0.15))
lines!(ax, 1:n_obs, rep.median[1:n_obs]; linewidth=2, color=:steelblue,
       label="Fitted")

if n_total > n_obs
    forecast_idx = (n_obs+1):n_total
    band!(ax, forecast_idx, rep.lower_90[forecast_idx], rep.upper_90[forecast_idx];
          color=(:firebrick, 0.15))
    band!(ax, forecast_idx, rep.lower_50[forecast_idx], rep.upper_50[forecast_idx];
          color=(:firebrick, 0.3))
    lines!(ax, forecast_idx, rep.median[forecast_idx]; linewidth=2,
           color=:firebrick, label="Forecast")
    vlines!(ax, [n_obs + 0.5]; color=:grey, linestyle=:dash)
end

ax.xticks = (1:10:n_total, string.(rep.date[1:10:n_total]))
ax.xticklabelrotation = π/4
axislegend(ax; position=:lt)
fig
```

The forecast uncertainty (red band) is wider than the fitted uncertainty — this correctly reflects that predicting the future is harder than explaining the past. This growing uncertainty is a fundamental feature of forecasting. For planning purposes, the width of the credible interval is as informative as the central estimate — it tells us the range of plausible outcomes.

## Accounting for incomplete observation

In many outbreaks, not all cases are reported. Surveillance captures only a fraction of true infections due to limited testing capacity, asymptomatic infections going undetected, and reporting delays. The observation model in EpiNow2 accounts for this by scaling estimated infections by a reporting fraction.

The `scale` parameter takes a `Normal` distribution to express uncertainty about the reporting fraction itself. For example, if we believe roughly 40% of cases are reported:

```@example forecast
result_scaled = epinow(cases;
    generation_time = gt_opts(gt),
    delays = delay_opts(delay),
    obs = obs_opts(scale=Normal(0.4, 0.01)),
    forecast = forecast_opts(horizon=7),
    inference = inference_opts(samples=1000, warmup=250, chains=3,
                               progress=false),
    verbose = false)
nothing # hide
```

!!! details "Coming from R?"
    In R's EpiNow2, the observation scale is specified as `obs_opts(scale = list(mean = 0.4, sd = 0.01))`. In Julia, we pass a `Normal(0.4, 0.01)` distribution object directly — Julia's type system makes this more explicit about what the scale parameter represents (a distribution, not just a list of numbers).

```@example forecast
inf1 = result.estimates.infections
inf2 = result_scaled.estimates.infections

fig = Figure(size=(700, 350))
ax = Axis(fig[1, 1]; xlabel="Date", ylabel="Estimated infections",
          title="Impact of observation scaling on infection estimates")

lines!(ax, 1:nrow(inf1), inf1.median; linewidth=2, color=:steelblue,
       label="100% observed")
lines!(ax, 1:nrow(inf2), inf2.median; linewidth=2, color=:firebrick,
       label="40% observed")

ax.xticks = (1:10:nrow(inf1), string.(inf1.date[1:10:end]))
ax.xticklabelrotation = π/4
axislegend(ax; position=:lt)
fig
```

With observation scaling, estimated **infections** are higher (by roughly 1/0.4 = 2.5×), but the **``R_t`` estimate** remains similar — scaling affects the level, not the trend of transmission.

This has important public health implications: if only 40% of cases are observed, the true burden is substantially higher than reported. This matters for hospital capacity planning and resource allocation. However, trend-based metrics like ``R_t`` and growth rate are robust to constant under-reporting.

!!! note "Discussion"
    Compare different observation percentages (e.g. 20%, 40%, 80%). How do the estimated infection numbers differ? What are the public health implications of assuming a higher or lower reporting fraction? Note that while the absolute number of estimated infections changes substantially, the reproduction number and growth rate estimates remain largely unchanged.

## Longer forecast horizon

For planning purposes, we might want a longer forecast. The reproduction number projected into the future can also be controlled — by default, the most recent ``R_t`` estimate (based on partial data, hence more uncertain) is used. An alternative is to use a less recent but more certain estimate, trading off recency for precision.

For a 14-day forecast:

```@example forecast
result_14 = epinow(cases;
    generation_time = gt_opts(gt),
    delays = delay_opts(delay),
    forecast = forecast_opts(horizon=14),
    inference = inference_opts(samples=1000, warmup=250, chains=3))
```

Longer horizons produce wider credible intervals — forecast uncertainty
grows with time.

## Forecasting secondary observations

EpiNow2.jl can also estimate and forecast **secondary observations** (e.g. deaths or hospitalisations) from primary observations (e.g. cases). This is useful for planning health system capacity.

The approach requires:
- A time series with both primary and secondary observations
- The delay distribution between primary and secondary events (e.g. case-to-death delay)
- A specification of the relationship type: **incidence** (secondary observations arise from previous primary observations, e.g. deaths from cases) or **prevalence** (secondary observations depend on both current primary and past secondary observations, e.g. hospital bed occupancy)

!!! warning "Time-scale caution"
    In the early stages of an outbreak, there can be substantial changes in testing and reporting. If testing practices change from one month to another, there will be bias in the model fit. Be cautious of the time-scale of data used for fitting and forecasting.

The secondary forecasting pipeline has two steps:

1. **Estimate the relationship** between primary and secondary observations using `estimate_secondary()` with historical data
2. **Forecast** secondary observations using `forecast_secondary()` with future primary case data (either observed or themselves forecast)

This functionality enables a complete pipeline from case data through to forecasts of downstream outcomes like deaths and hospital admissions.

## Challenge: Ebola outbreak analysis

Using the Ebola 1976 data and delay distributions from `EpiParameters.jl`:

1. Estimate whether cases are increasing or decreasing at the end of the available data
2. Account for a capacity to observe 80% of cases
3. Create a 14-day forecast of case numbers

You can use the following approach:

```julia
using CSV

ebola = CSV.read(joinpath(@__DIR__, "ebola1976.csv"), DataFrame)
rename!(ebola, :cases => :confirm)

# Get Ebola delay distributions from the literature
ebola_delays = epiparameter(disease="Ebola", epi_name="onset to death")
ebola_gt = discretise(
    filter(p -> !isnothing(p.distribution), ebola_delays)[1].distribution;
    max=20
)

ebola_result = epinow(ebola;
    generation_time = gt_opts(ebola_gt),
    forecast = forecast_opts(horizon=14),
    obs = obs_opts(scale=Normal(0.8, 0.01)),
    inference = inference_opts(samples=1000, warmup=250, chains=3))
```

Before running the model, consider:
- What assumptions are we making by projecting the current ``R_t`` forward?
- In what situations might this assumption break down? (e.g. if interventions are about to be introduced, or if a new variant is emerging)
- How would you interpret the width of the forecast credible intervals for public health planning?

### Controlling the forecast Rt

By default, the most recent ``R_t`` estimate is projected forward. Since this estimate is based on partial data, it has considerable uncertainty. The reproduction number used for forecasting can be controlled:

- **Latest estimate** (default): uses the most recent ``R_t``, which has wider uncertainty
- **Stable estimate**: uses a less recent but more certain ``R_t`` estimate — this trades off recency for precision

The choice matters: using the latest estimate captures the most current transmission dynamics but with more uncertainty, whilst using an earlier estimate provides tighter forecasts that may not reflect recent changes.

## Key points

- EpiNow2.jl automatically produces **short-term forecasts** by projecting ``R_t`` forward
- Forecasts are generated as part of the MCMC posterior — uncertainty is properly propagated
- **Observation scaling** (`obs_opts(scale=...)`) accounts for under-reporting
- Forecast uncertainty **grows with horizon** — always report credible intervals
- The **forecast horizon** is controlled via `forecast_opts(horizon=N)`
- **Secondary observations** (deaths, hospitalisations) can be forecast from primary case data
- Consider different **observation fractions** and their impact on absolute estimates vs trend metrics
- The Rt projection used for forecasting can be controlled to trade off recency vs precision

---

*These tutorials are adapted from the [Epiverse-TRACE tutorials](https://epiverse-trace.github.io/tutorials/), © Epiverse-TRACE contributors, licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).*
