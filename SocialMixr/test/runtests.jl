using Test
using SocialMixr
using DataFrames
using ContactMatrices

@testset "SocialMixr" begin
    include("test_age_groups.jl")
    include("test_compute_matrix.jl")
    include("test_weights.jl")
    include("test_postprocess.jl")
end
