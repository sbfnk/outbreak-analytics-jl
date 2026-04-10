using Test

@testset "ContactMatrices" begin
    include("test_constructors.jl")
    include("test_symmetry.jl")
    include("test_operations.jl")
end
