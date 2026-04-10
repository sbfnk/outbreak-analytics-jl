# Validate Case Data

After cleaning, we need to **validate** that our data meets expectations:
correct types, required columns present, and values within plausible ranges.
This tutorial shows how to build lightweight validation checks in Julia.

## Learning objectives

- Check that required columns exist with correct types
- Validate value ranges and date sequences
- Tag columns with epidemiological meaning
- Build reusable validation functions

## A validated linelist

After cleaning (previous tutorial), we expect our data to have specific
columns with specific types:

```@example validate
using DataFrames
using Dates

# Cleaned data (output of previous tutorial)
data = DataFrame(
    case_id = 1:8,
    date_onset = Date.(["2024-01-15", "2024-01-15", "2024-01-17",
                         "2024-01-18", "2024-01-19", "2024-01-20",
                         "2024-01-21", "2024-01-23"]),
    age = [25, 34, missing, 45, 12, 67, 28, 55],
    gender = ["male", "female", "male", "female", "male",
              "female", missing, "male"],
    status = ["confirmed", "confirmed", "suspected", "confirmed",
              "confirmed", "probable", "confirmed", "confirmed"],
)
data
```

## Column validation

Check that required columns exist and have the right types:

```@example validate
function validate_columns(df; required=Dict{Symbol, Type}())
    issues = String[]
    for (col, expected_type) in required
        if !(String(col) in names(df))
            push!(issues, "Missing column: $col")
        else
            actual = nonmissingtype(eltype(df[!, col]))
            if !(actual <: expected_type)
                push!(issues, "$col: expected $expected_type, got $actual")
            end
        end
    end
    if isempty(issues)
        println("✓ All columns valid")
    else
        for issue in issues
            println("✗ $issue")
        end
    end
    isempty(issues)
end

validate_columns(data;
    required = Dict(
        :case_id => Integer,
        :date_onset => Date,
        :age => Real,
        :gender => AbstractString,
        :status => AbstractString,
    ))
```

## Value range checks

```@example validate
function validate_values(df)
    issues = String[]

    # Age range
    if "age" in names(df)
        ages = skipmissing(df.age)
        if any(a -> a < 0 || a > 120, ages)
            push!(issues, "Age values outside [0, 120]")
        end
    end

    # Date range
    if "date_onset" in names(df)
        dates = skipmissing(df.date_onset)
        if any(d -> d > today(), dates)
            push!(issues, "Onset dates in the future")
        end
    end

    # Status values
    if "status" in names(df)
        valid_status = ["confirmed", "probable", "suspected"]
        invalid = setdiff(unique(skipmissing(df.status)), valid_status)
        if !isempty(invalid)
            push!(issues, "Invalid status values: $invalid")
        end
    end

    if isempty(issues)
        println("✓ All values valid")
    else
        for issue in issues
            println("✗ $issue")
        end
    end
    isempty(issues)
end

validate_values(data)
```

## Tagging columns with epidemiological meaning

In R, the `linelist` package tags columns with semantic roles (`id`,
`date_onset`, `gender`). In Julia, we can use a lightweight wrapper:

```@example validate
struct Linelist
    data::DataFrame
    tags::Dict{Symbol, Symbol}
end

function Linelist(df::DataFrame; tags...)
    tag_dict = Dict{Symbol, Symbol}(k => v for (k, v) in tags)
    # Validate all tagged columns exist
    for (role, col) in tag_dict
        @assert String(col) in names(df) "Tagged column :$col not found in DataFrame"
    end
    Linelist(df, tag_dict)
end

Base.getindex(ll::Linelist, tag::Symbol) = ll.data[!, ll.tags[tag]]
DataFrames.nrow(ll::Linelist) = nrow(ll.data)
```

```@example validate
ll = Linelist(data;
    id = :case_id,
    date_onset = :date_onset,
    age = :age,
    gender = :gender,
)

println("Cases: $(nrow(ll))")
println("Onset dates: $(extrema(skipmissing(ll[:date_onset])))")
println("Age range: $(extrema(skipmissing(ll[:age])))")
```

## Completeness summary

```@example validate
function completeness(df)
    result = DataFrame(
        column = names(df),
        n_total = nrow(df),
        n_present = [count(!ismissing, df[!, c]) for c in names(df)],
    )
    result.pct_complete = round.(result.n_present ./ result.n_total .* 100; digits=1)
    sort(result, :pct_complete)
end

completeness(data)
```

## Key points

- **Validate columns** exist with expected types before analysis
- **Check value ranges** — ages, dates, categorical levels
- **Tag columns** with epidemiological roles for safer downstream analysis
- **Report completeness** — know your missing data rates
- Julia's type system and `missing` make validation straightforward with
  plain DataFrames.jl — no dedicated validation package needed
