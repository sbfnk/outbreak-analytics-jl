"""
    CfrEstimate

A point estimate of the case fatality ratio with confidence interval.

# Fields
- `estimate::Float64` — point estimate
- `lower::Float64` — lower bound of the confidence interval
- `upper::Float64` — upper bound of the confidence interval
"""
struct CfrEstimate
    estimate::Float64
    lower::Float64
    upper::Float64
end

function Base.show(io::IO, x::CfrEstimate)
    print(io, "CFR estimate: ",
          round(x.estimate; digits=4),
          " (", round(x.lower; digits=4),
          ", ", round(x.upper; digits=4), ")")
end
