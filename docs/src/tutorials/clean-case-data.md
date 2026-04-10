# Clean Case Data

Real outbreak data is messy: inconsistent column names, mixed date formats,
duplicates, missing values encoded as strings, and data entry errors.
In Julia, `DataFrames.jl` provides all the tools needed for data cleaning.

## Learning objectives

- Standardise column names
- Handle missing values and duplicates
- Parse and validate dates
- Recode categorical variables

## Example: messy outbreak data

Let's create a realistic messy dataset and clean it step by step:

```@example cleandata
using DataFrames
using Dates
using CairoMakie

# Simulate messy data (mimicking real surveillance exports)
messy = DataFrame(
    "Case ID" => [1, 2, 3, 4, 5, 5, 6, 7, 8, 9],
    "Date_Onset" => ["2024-01-15", "15/01/2024", "2024-01-17",
                      "", "2024-01-19", "2024-01-19", "2024-01-20",
                      "2024-01-21", "Jan 22, 2024", "2024-01-23"],
    "age" => [25, 34, missing, 45, 12, 12, 67, 28, 41, 55],
    "Gender" => ["M", "female", "m", "F", "Male", "Male", "f", "", "FEMALE", "M"],
    "STATUS" => ["confirmed", "confirmed", "suspected", "confirmed",
                 "confirmed", "confirmed", "probable", "confirmed",
                 "confirmed", "confirmed"],
    "Lab" => fill("Central", 10),  # constant column
)
messy
```

## Step 1: Standardise column names

Convert to lowercase with underscores:

```@example cleandata
function standardise_names!(df)
    new_names = [Symbol(lowercase(replace(string(n), r"[\s\-]+" => "_")))
                 for n in names(df)]
    rename!(df, names(df) .=> new_names)
end

standardise_names!(messy)
names(messy)
```

## Step 2: Remove constant columns

Columns with a single unique value carry no information:

```@example cleandata
constant_cols = [n for n in names(messy) if length(unique(skipmissing(messy[!, n]))) <= 1]
println("Constant columns: $constant_cols")
select!(messy, Not(constant_cols))
names(messy)
```

## Step 3: Remove duplicates

```@example cleandata
println("Rows before: $(nrow(messy))")
unique!(messy)
println("Rows after: $(nrow(messy))")
```

## Step 4: Handle missing values

Replace empty strings with `missing`:

```@example cleandata
for col in names(messy)
    if eltype(messy[!, col]) <: Union{Missing, AbstractString}
        messy[!, col] = [x === missing || x == "" ? missing : x
                         for x in messy[!, col]]
    end
end

println("Missing values per column:")
for col in names(messy)
    n = count(ismissing, messy[!, col])
    n > 0 && println("  $col: $n")
end
```

## Step 5: Parse dates

Real data has mixed date formats. Parse with fallbacks:

```@example cleandata
function parse_date_flexible(s)
    ismissing(s) && return missing
    for fmt in [dateformat"yyyy-mm-dd", dateformat"dd/mm/yyyy", dateformat"u dd, yyyy"]
        try
            return Date(s, fmt)
        catch
            continue
        end
    end
    return missing  # unparseable
end

messy.date_onset = parse_date_flexible.(messy.date_onset)
messy
```

## Step 6: Standardise categorical variables

Recode gender to a consistent format:

```@example cleandata
function standardise_gender(g)
    ismissing(g) && return missing
    g_lower = lowercase(strip(g))
    g_lower in ["m", "male"] && return "male"
    g_lower in ["f", "female"] && return "female"
    return missing
end

messy.gender = standardise_gender.(messy.gender)
messy
```

## Step 7: Validate date sequences

Check that onset dates are reasonable:

```@example cleandata
valid_dates = dropmissing(messy, :date_onset)
date_range = extrema(valid_dates.date_onset)
println("Date range: $(date_range[1]) to $(date_range[2])")

# Flag any dates outside expected range
expected_start = Date(2024, 1, 1)
expected_end = Date(2024, 12, 31)
outside = filter(r -> !ismissing(r.date_onset) &&
                      (r.date_onset < expected_start || r.date_onset > expected_end),
                 messy)
println("Dates outside expected range: $(nrow(outside))")
```

## Building a cleaning pipeline

Combine all steps into a reusable function:

```julia
function clean_linelist(raw::DataFrame)
    df = copy(raw)
    standardise_names!(df)

    # Remove constant columns
    constant = [n for n in names(df) if length(unique(skipmissing(df[!, n]))) <= 1]
    select!(df, Not(constant))

    # Remove duplicates
    unique!(df)

    # Missing values
    for col in names(df)
        if eltype(df[!, col]) <: Union{Missing, AbstractString}
            df[!, col] = [x === missing || x == "" ? missing : x for x in df[!, col]]
        end
    end

    # Parse dates (adapt column names as needed)
    if "date_onset" in names(df)
        df.date_onset = parse_date_flexible.(df.date_onset)
    end

    df
end
```

## Key points

- **Standardise column names** to lowercase with underscores
- **Remove constant columns** and **duplicates** early
- Replace empty strings with `missing` — Julia's `missing` propagates correctly
- **Parse dates flexibly** with multiple format fallbacks
- **Recode categorical variables** to consistent values
- Build a reusable **cleaning pipeline** for your project
- Unlike R's `cleanepi`, Julia uses DataFrames.jl directly — no extra package needed
