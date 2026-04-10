module SocialMixr

using CSV
using DataFrames
using LinearAlgebra

using ContactMatrices

export ContactSurvey
export polymod, polymod_population, load_population
export reduce_agegroups, limits_to_agegroups, assign_age_groups
export weigh
export filter_survey, compute_matrix
export symmetrise, split_matrix, per_capita, pop_age

include("types.jl")
include("data.jl")
include("age_groups.jl")
include("weights.jl")
include("compute_matrix.jl")
include("postprocess.jl")

end # module SocialMixr
