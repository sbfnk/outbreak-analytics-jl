# Aggregate and Visualise

Before diving into transmission modelling, we need to turn raw case data into
**incidence curves** — counts of cases aggregated by time period. In R, the
`incidence2` package provides this; in Julia, we use `DataFrames.jl` groupby
operations directly.

## Learning objectives

- Load case data from CSV
- Aggregate cases by day, week, or custom periods
- Plot epidemic curves with CairoMakie

## Loading data

```@example aggregate
using DataFrames
using Dates
using CSV
using CairoMakie

ebola = CSV.read(joinpath(@__DIR__, "ebola1976.csv"), DataFrame)
first(ebola, 5)
```

## Daily epidemic curve

The simplest visualisation — cases per day:

```@example aggregate
fig = Figure(size=(700, 350))
ax = Axis(fig[1, 1]; xlabel="Date", ylabel="Daily cases",
          title="Ebola 1976 — daily incidence")
barplot!(ax, 1:nrow(ebola), ebola.cases; color=(:steelblue, 0.7))
ax.xticks = (1:14:nrow(ebola), string.(ebola.date[1:14:end]))
ax.xticklabelrotation = π/4
fig
```

## Weekly aggregation

Daily counts can be noisy. Aggregate to weekly:

```@example aggregate
ebola.week = Dates.value.(ebola.date .- ebola.date[1]) .÷ 7

weekly = combine(groupby(ebola, :week),
                 :cases => sum => :cases,
                 :deaths => sum => :deaths,
                 :date => first => :week_start)

weekly
```

```@example aggregate
fig = Figure(size=(700, 350))
ax = Axis(fig[1, 1]; xlabel="Week starting", ylabel="Count",
          title="Ebola 1976 — weekly incidence")
barplot!(ax, weekly.week, weekly.cases; width=0.4, offset=-0.2,
         color=(:steelblue, 0.7), label="Cases")
barplot!(ax, weekly.week, weekly.deaths; width=0.4, offset=0.2,
         color=(:firebrick, 0.7), label="Deaths")
ax.xticks = (weekly.week[1:2:end], string.(weekly.week_start[1:2:end]))
ax.xticklabelrotation = π/4
axislegend(ax; position=:rt)
fig
```

## Cumulative incidence

```@example aggregate
ebola.cumulative_cases = cumsum(ebola.cases)
ebola.cumulative_deaths = cumsum(ebola.deaths)

fig = Figure(size=(700, 350))
ax = Axis(fig[1, 1]; xlabel="Date", ylabel="Cumulative count",
          title="Ebola 1976 — cumulative incidence")
lines!(ax, 1:nrow(ebola), ebola.cumulative_cases; linewidth=2,
       color=:steelblue, label="Cases")
lines!(ax, 1:nrow(ebola), ebola.cumulative_deaths; linewidth=2,
       color=:firebrick, label="Deaths")
ax.xticks = (1:14:nrow(ebola), string.(ebola.date[1:14:end]))
ax.xticklabelrotation = π/4
axislegend(ax; position=:lt)
fig
```

## The Julia approach: DataFrames.jl

In R, `incidence2` provides specialised incidence objects. In Julia, we use
`DataFrames.jl` directly — the same `groupby` + `combine` pattern handles
any aggregation:

```julia
# By week
combine(groupby(data, :week), :cases => sum)

# By month
data.month = Dates.month.(data.date)
combine(groupby(data, :month), :cases => sum)

# By age group (if available)
combine(groupby(data, [:week, :age_group]), :cases => sum)
```

No dedicated package needed — `DataFrames.jl` is the Swiss army knife.

## Key points

- **Epidemic curves** are the foundation of outbreak analytics — always plot
  your data before modelling
- Use `DataFrames.groupby` + `combine` for aggregation by any time period
- **Weekly aggregation** smooths daily noise while preserving epidemic shape
- **Cumulative curves** show total burden and whether the outbreak is slowing
- Julia's `DataFrames.jl` replaces R's `incidence2` — same `groupby` pattern,
  no extra package
