module CFR

using DataFrames
using Distributions

export CfrEstimate, cfr_static, cfr_time_varying, estimate_outcomes

include("types.jl")
include("outcomes.jl")
include("static.jl")
include("time_varying.jl")

end
