# Read Case Data

The first step in outbreak analytics is getting data into your analysis
environment. In Julia, `CSV.jl` and `DataFrames.jl` handle the vast majority
of tabular data formats.

## Learning objectives

- Read case data from CSV files
- Handle common file format issues (delimiters, missing values, dates)
- Inspect and summarise the data

## Reading a CSV file

```@example readdata
using CSV
using DataFrames
using Dates

ebola = CSV.read(joinpath(@__DIR__, "ebola1976.csv"), DataFrame)
first(ebola, 5)
```

```@example readdata
println("Rows: $(nrow(ebola)), Columns: $(ncol(ebola))")
println("Column types: $(eltype.(eachcol(ebola)))")
```

## Handling common issues

### Missing value strings

Real-world data often uses `""`, `"NA"`, `"N/A"`, or `"unknown"` for missing
values. Specify these at read time:

```@example readdata
# CSV.read("file.csv", DataFrame; missingstring=["", "NA", "N/A", "unknown"])
nothing # hide
```

### Date parsing

CSV.jl auto-detects ISO 8601 dates (`YYYY-MM-DD`). For other formats, specify
the `dateformat`:

```@example readdata
# For DD/MM/YYYY format:
# CSV.read("file.csv", DataFrame; dateformat="dd/mm/yyyy")

# For mixed formats, read as String then parse:
dates_str = ["2024-01-15", "2024-02-20", "2024-03-10"]
dates = Date.(dates_str)
```

### Delimiters

Tab-separated, semicolon-separated, or other delimiters:

```@example readdata
# CSV.read("file.tsv", DataFrame; delim='\t')
# CSV.read("file.csv", DataFrame; delim=';')
nothing # hide
```

## Inspecting the data

```@example readdata
describe(ebola)
```

Check for missing values:

```@example readdata
for col in names(ebola)
    n_missing = count(ismissing, ebola[!, col])
    n_missing > 0 && println("$col: $n_missing missing values")
end
println("No missing values found") # hide
```

## Reading other formats

| Format | Julia package |
|--------|--------------|
| CSV / TSV | `CSV.jl` |
| Excel (.xlsx) | `XLSX.jl` |
| JSON | `JSON3.jl` |
| Compressed (.gz) | `CSV.jl` handles transparently |
| SQLite | `SQLite.jl` + `DBInterface.jl` |

## Key points

- `CSV.read(path, DataFrame)` is the workhorse for tabular data
- Handle missing values with `missingstring=` at read time
- Specify `dateformat=` for non-ISO date formats
- Use `describe(df)` and column-level checks to inspect data quality
- Julia's `missing` type is first-class — no special handling needed
