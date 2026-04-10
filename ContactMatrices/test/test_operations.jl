using ContactMatrices
using Test

@testset "Addition" begin
    cm1 = ContactMatrix([1.0 2.0; 3.0 4.0], ["a", "b"]; setting = :home)
    cm2 = ContactMatrix([0.5 0.5; 0.5 0.5], ["a", "b"]; setting = :work)
    cm_sum = cm1 + cm2

    @test cm_sum[1, 1] ≈ 1.5
    @test cm_sum[1, 2] ≈ 2.5
    @test cm_sum[2, 1] ≈ 3.5
    @test cm_sum[2, 2] ≈ 4.5
    @test setting(cm_sum) == :combined
end

@testset "Addition mismatched groupings" begin
    cm1 = ContactMatrix([1.0 2.0; 3.0 4.0], ["a", "b"])
    cm2 = ContactMatrix([1.0 2.0; 3.0 4.0], ["x", "y"])
    @test_throws ArgumentError cm1 + cm2
end

@testset "Scalar multiplication" begin
    cm = ContactMatrix([1.0 2.0; 3.0 4.0], ["a", "b"])
    cm2 = 2.0 * cm
    @test cm2[1, 1] ≈ 2.0
    @test cm2[2, 2] ≈ 8.0
    @test setting(cm2) == :all

    cm3 = cm * 3.0
    @test cm3[1, 1] ≈ 3.0
end

@testset "Matrix extraction" begin
    mat = [1.0 2.0; 3.0 4.0]
    cm = ContactMatrix(mat, ["a", "b"])
    extracted = Matrix(cm)
    @test extracted == mat
    @test extracted !== cm.data  # should be a copy
end

@testset "reduce_groups simple sum" begin
    cm = ContactMatrix(
        [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0],
        ["[0,5)", "[10,15)", "[5,10)"]
    )
    mapping = Dict(
        "[0,5)" => "[0,10)",
        "[5,10)" => "[0,10)",
        "[10,15)" => "[10,15)"
    )
    coarse = reduce_groups(cm, mapping)
    @test size(coarse) == (2, 2)
    labels = groupings(coarse).labels[1]
    @test labels == ["[0,10)", "[10,15)"]

    # [0,10),[0,10) = [0,5),[0,5) + [0,5),[5,10) + [5,10),[0,5) + [5,10),[5,10)
    # Note: labels are sorted, so [0,5) is idx 1, [10,15) is idx 2, [5,10) is idx 3
    @test coarse[1, 1] ≈ 1.0 + 3.0 + 7.0 + 9.0
end

@testset "reduce_groups population-weighted" begin
    # Labels sorted: old=1, young=2
    cm = ContactMatrix([4.0 3.0; 1.0 2.0], ["old", "young"])
    mapping = Dict("old" => "all", "young" => "all")
    pop = [500.0, 1000.0]  # old=500, young=1000
    coarse = reduce_groups(cm, mapping; population = pop)
    @test size(coarse) == (1, 1)
    # Total contacts = sum of pop[i] * c[i,j] for all i,j, divided by total pop
    total_contacts = 500.0 * 4.0 + 500.0 * 3.0 + 1000.0 * 1.0 + 1000.0 * 2.0
    total_pop = 1500.0
    @test coarse[1, 1] ≈ total_contacts / total_pop
end

@testset "AbstractArray interface" begin
    cm = ContactMatrix([1.0 2.0; 3.0 4.0], ["a", "b"])
    @test length(cm) == 4
    @test eltype(cm) == Float64
    @test ndims(cm) == 2

    # Iteration
    vals = vec(collect(cm))
    @test sort(vals) == [1.0, 2.0, 3.0, 4.0]

    # Broadcasting
    doubled = cm .* 2
    @test doubled == [2.0 4.0; 6.0 8.0]
end
