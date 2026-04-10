### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4d03
begin
    import Pkg
    Pkg.develop(path=joinpath(@__DIR__, "..", "EpiParameters"))
    Pkg.develop(path=expanduser("~/code/EpiNow2.jl"))
    using EpiNow2
    using EpiParameters
    using Distributions
    using DataFrames
    using Dates
    using CSV
    using CairoMakie
end

# ╔═╡ 1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5d03
md"""
# Create a Short-Term Forecast

Beyond estimating $R_t$, we often need to **forecast** — predict future case
counts to inform resource planning and public health response. Given case data,
we can create estimates of current and future case numbers by accounting for
delays in reporting and under-reporting.

To make predictions about the future, we need an assumption about how the
present relates to what we expect. The simplest assumption is "no change" — the
reproduction number remains the same as its most recent estimate. This is the
approach EpiNow2 uses by default.

Under the hood, `epinow()` is a wrapper for two underlying functions:
`estimate_infections()` (which estimates cases by date of infection from
reported data) and `forecast_infections()` (which simulates future infections
using the fitted model). When you call `epinow()`, both steps are performed
automatically.

This tutorial covers:
1. Generating a 7-day forecast from case data
2. Adjusting for **incomplete observation** (reporting fraction)
3. Customising the **forecast horizon**
4. Interpreting forecast uncertainty

!!! warning "Computation time"
    MCMC sampling takes several minutes per model fit.
"""

# ╔═╡ 2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6d03
md"""
## Basic forecast

We start with the synthetic example data and produce a default 7-day forecast.
The `forecast_opts(horizon=7)` argument specifies that we want to forecast 7
days beyond the last observed data point. The forecast is generated as part of
the MCMC sampling — for each posterior sample of $R_t$, future infections are
simulated, naturally propagating uncertainty from the estimation into the
forecast.
"""

# ╔═╡ 3a4b5c6d-7e8f-9a0b-1c2d-3e4f5a6b7d03
cases = example_confirmed()

# ╔═╡ 4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8d03
begin
    gt = discretise(LogNormal(1.2, 0.5); max=14)
    delay = discretise(LogNormal(1.4, 0.4); max=14)

    result = epinow(cases;
        generation_time = gt_opts(gt),
        delays = delay_opts(delay),
        forecast = forecast_opts(horizon=7),
        inference = inference_opts(samples=1000, warmup=200, chains=2,
                                   progress=false),
        verbose = false)
end

# ╔═╡ 5a6b7c8d-9e0f-1a2b-3c4d-5e6f7a8b9d03
md"""
### Visualising the forecast

The reports DataFrame includes both fitted values (for observed dates) and
forecasts (for future dates). The estimates in the output are split into three
categories: **Estimate** uses all available data and has the tightest credible
intervals. **Estimate based on partial data** has higher uncertainty because
recent estimates are informed by fewer observations — right-truncation of the
delay distribution means recent data are incomplete. **Forecast** projects into
the future using the estimated $R_t$ trajectory. The grey dashed line in the
plot below separates the observation period from the forecast period.
"""

# ╔═╡ 6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0d03
let
    rep = result.estimates.reports
    n_obs = nrow(cases)
    n_total = nrow(rep)

    fig = Figure(size=(750, 400))
    ax = Axis(fig[1, 1]; xlabel="Date", ylabel="Cases",
              title="Case forecast (7-day horizon)")

    # Observed data
    barplot!(ax, 1:n_obs, cases.confirm; color=(:grey, 0.5), label="Observed")

    # Fitted (observation period)
    band!(ax, 1:n_obs, rep.lower_90[1:n_obs], rep.upper_90[1:n_obs];
          color=(:steelblue, 0.15))
    lines!(ax, 1:n_obs, rep.median[1:n_obs]; linewidth=2, color=:steelblue,
           label="Fitted")

    # Forecast (future)
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
end

# ╔═╡ 7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1d03
md"""
Note the forecast uncertainty (red band) is wider than the fitted uncertainty —
this correctly reflects that predicting the future is harder than explaining
the past. This growing uncertainty is a fundamental feature of forecasting — it
reflects genuine epistemic uncertainty about future transmission. Shorter
forecast horizons produce tighter intervals. For planning purposes, the width
of the credible interval is as informative as the central estimate: a wide
interval signals that decision-makers should prepare for a range of scenarios.
"""

# ╔═╡ 8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2d03
md"""
## Accounting for incomplete observation

In reality, 100% of cases are rarely reported — surveillance captures only a
fraction of true infections. Common reasons include: limited testing capacity,
asymptomatic infections going undetected, and reporting delays that cause cases
to be missed entirely. The observation model in EpiNow2 accounts for this by
scaling the estimated infections by a reporting fraction.

We can specify an **observation scaling** factor using `obs_opts(scale=...)`.
The `scale` parameter takes a `Normal` distribution to express uncertainty
about the reporting fraction itself. For example, if we believe roughly 40% of
cases are reported, with some uncertainty:
"""

# ╔═╡ 9a0b1c2d-3e4f-5a6b-7c8d-9e0f1a2b3d03
result_scaled = epinow(cases;
    generation_time = gt_opts(gt),
    delays = delay_opts(delay),
    obs = obs_opts(scale=Normal(0.4, 0.01)),
    forecast = forecast_opts(horizon=7),
    inference = inference_opts(samples=1000, warmup=200, chains=2,
                               progress=false),
    verbose = false)

# ╔═╡ aa1b2c3d-4e5f-6a7b-8c9d-ae1f2a3b4d03
md"""
With observation scaling, the estimated **infections** will be higher than
reported cases (by roughly 1/0.4 = 2.5×), but the **$R_t$ estimate** remains
similar — scaling affects the level, not the trend of transmission.
"""

# ╔═╡ ba2b3c4d-5e6f-7a8b-9c0d-be2f3a4b5d03
let
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
end

# ╔═╡ f1a2b3c4-d5e6-7f8a-9b0c-1d2e3f4a5b6c
md"""
The public health implications of under-reporting are significant: if only 40%
of cases are observed, the true burden is roughly 2.5× higher than reported.
This has direct consequences for hospital capacity planning, resource
allocation, and understanding the true scope of an outbreak. However,
trend-based metrics like $R_t$ and growth rate are robust to constant
under-reporting — a reassuring property that means we can track whether an
epidemic is growing or shrinking even without perfect surveillance.
"""

# ╔═╡ ca3b4c5d-6e7f-8a9b-0c1d-ce3f4a5b6d03
md"""
## Longer forecast horizon

For planning purposes, we might want a 14-day forecast. The forecast of $R_t$
into the future can be controlled. By default, the most recent $R_t$ estimate
(which is based on partial data and therefore more uncertain) is projected
forward. An alternative is to use a less recent but more certain $R_t$ estimate
as the basis for projection. This trades off recency for precision — a choice
that depends on how rapidly you expect the epidemiological situation to change.

```julia
result_14 = epinow(cases;
    generation_time = gt_opts(gt),
    delays = delay_opts(delay),
    forecast = forecast_opts(horizon=14),
    inference = inference_opts(samples=1000, warmup=200, chains=2),
    verbose = false)
```

Longer horizons naturally produce wider credible intervals — forecast
uncertainty grows with time.
"""

# ╔═╡ da4b5c6d-7e8f-9a0b-1c2d-de4f5a6b7d03
md"""
## Challenge: Ebola forecast

Before attempting the code below, think about: What assumptions are we making
by projecting the current $R_t$ forward? In what situations might this
assumption break down? For example, if interventions are about to be introduced
or lifted, or if a new variant is emerging, the "no change" assumption would be
misleading. Forecasts are conditional on the model assumptions — they tell us
what *would* happen if current trends continue, not what *will* happen.

Using the Ebola 1976 data and delay distributions from `EpiParameters.jl`:

1. Find the Ebola onset-to-death delay and use it as a proxy for the
   generation time (or find a serial interval estimate)
2. Estimate $R_t$ — is the epidemic growing at the end of the dataset?
3. Produce a 14-day forecast

```julia
ebola = CSV.read("ebola1976.csv", DataFrame)
rename!(ebola, :cases => :confirm)

ebola_delays = epiparameter(disease="Ebola", epi_name="onset to death")
ebola_gt = discretise(
    filter(p -> !isnothing(p.distribution), ebola_delays)[1].distribution;
    max=20
)

epinow(ebola;
    generation_time = gt_opts(ebola_gt),
    forecast = forecast_opts(horizon=14),
    obs = obs_opts(scale=Normal(0.8, 0.01)),
    inference = inference_opts(samples=1000, warmup=200, chains=2))
```
"""

# ╔═╡ ea5b6c7d-8e9f-0a1b-2c3d-ee5f6a7b8d03
md"""
## Key points

- EpiNow2.jl automatically produces **short-term forecasts** by projecting
  the estimated $R_t$ trajectory forward
- Forecasts are generated as part of the MCMC posterior — uncertainty is
  properly propagated
- **Observation scaling** (`obs_opts(scale=...)`) accounts for under-reporting
  without affecting $R_t$ estimates
- Forecast uncertainty **grows with horizon** — always report credible intervals
- The **forecast horizon** is controlled via `forecast_opts(horizon=N)`
- Combine with `EpiParameters.jl` for a complete literature-to-forecast pipeline
"""

# ╔═╡ Cell order:
# ╟─1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5d03
# ╠═0a1b2c3d-4e5f-6a7b-8c9d-0e1f2a3b4d03
# ╟─2a3b4c5d-6e7f-8a9b-0c1d-2e3f4a5b6d03
# ╠═3a4b5c6d-7e8f-9a0b-1c2d-3e4f5a6b7d03
# ╠═4a5b6c7d-8e9f-0a1b-2c3d-4e5f6a7b8d03
# ╟─5a6b7c8d-9e0f-1a2b-3c4d-5e6f7a8b9d03
# ╠═6a7b8c9d-0e1f-2a3b-4c5d-6e7f8a9b0d03
# ╟─7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1d03
# ╟─8a9b0c1d-2e3f-4a5b-6c7d-8e9f0a1b2d03
# ╠═9a0b1c2d-3e4f-5a6b-7c8d-9e0f1a2b3d03
# ╟─aa1b2c3d-4e5f-6a7b-8c9d-ae1f2a3b4d03
# ╠═ba2b3c4d-5e6f-7a8b-9c0d-be2f3a4b5d03
# ╟─f1a2b3c4-d5e6-7f8a-9b0c-1d2e3f4a5b6c
# ╟─ca3b4c5d-6e7f-8a9b-0c1d-ce3f4a5b6d03
# ╟─da4b5c6d-7e8f-9a0b-1c2d-de4f5a6b7d03
# ╟─ea5b6c7d-8e9f-0a1b-2c3d-ee5f6a7b8d03
