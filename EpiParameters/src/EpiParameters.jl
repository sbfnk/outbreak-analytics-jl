module EpiParameters

using Distributions
using JSON3

export EpiParam, epiparameter, list_diseases, list_parameters

include("types.jl")
include("database.jl")
include("query.jl")

end
