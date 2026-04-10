using ContactMatrices
using Test

@testset "Tidy-format single grouping" begin
    cm = ContactMatrix(
        of    = ["[0,5)", "[5,10)", "[5,10)"],
        with  = ["[0,5)", "[10,15)", "[15,20)"],
        value = [0.32, 0.46, 0.72]
    )
    @test size(cm) == (4, 4)
    # String-sorted labels: [0,5)=1, [10,15)=2, [15,20)=3, [5,10)=4
    @test groupings(cm).labels[1] == ["[0,5)", "[10,15)", "[15,20)", "[5,10)"]
    @test cm[1, 1] == 0.32   # [0,5) of, [0,5) with
    @test cm[4, 2] == 0.46   # [5,10) of, [10,15) with
    @test cm[4, 3] == 0.72   # [5,10) of, [15,20) with
    @test cm[3, 1] == 0.0    # unfilled
    @test setting(cm) == :all
end

@testset "Tidy-format multi grouping" begin
    cm = ContactMatrix(
        of    = (age = ["young", "young", "old"], sex = ["male", "female", "female"]),
        with  = (age = ["old", "old", "young"], sex = ["female", "female", "female"]),
        value = [1.0, 2.0, 2.0]
    )
    @test ndims(cm) == 4  # 2 groupings → 4D
    @test ndimgroups(groupings(cm)) == 2
    @test groupings(cm).names == [:age, :sex]
    @test setting(cm) == :all
end

@testset "From-matrix constructor" begin
    mat = [0.32 0.0; 0.0 0.46]
    cm = ContactMatrix(mat, ["[0,5)", "[5,10)"])
    @test size(cm) == (2, 2)
    @test cm[1, 1] == 0.32
    @test cm[2, 2] == 0.46
    @test cm[1, 2] == 0.0

    # Setting keyword
    cm2 = ContactMatrix(mat, ["[0,5)", "[5,10)"]; setting = :home)
    @test setting(cm2) == :home
end

@testset "From-matrix sorts labels" begin
    mat = [4.0 2.0; 1.0 3.0]
    cm = ContactMatrix(mat, ["beta", "alpha"])
    @test groupings(cm).labels[1] == ["alpha", "beta"]
    # After sorting: alpha→row1, beta→row2
    # Original: beta(row1)→alpha(row2) mapping
    @test cm[1, 1] == 3.0  # alpha,alpha was mat[2,2]
    @test cm[2, 2] == 4.0  # beta,beta was mat[1,1]
    @test cm[1, 2] == 1.0  # alpha,beta was mat[2,1]
    @test cm[2, 1] == 2.0  # beta,alpha was mat[1,2]
end

@testset "Constructor validation" begin
    @test_throws ArgumentError ContactMatrix(
        of = ["a", "b"], with = ["a"], value = [1.0, 2.0]
    )
    @test_throws ArgumentError ContactMatrix(
        [1.0 2.0; 3.0 4.0], ["a", "a"]  # duplicate labels
    )
    @test_throws ArgumentError ContactMatrix(
        [1.0 2.0 3.0; 4.0 5.0 6.0], ["a", "b"]  # non-square
    )
end

@testset "Fill value" begin
    cm = ContactMatrix(
        of = ["a"], with = ["a"], value = [5.0], fill = -1.0
    )
    @test cm[1, 1] == 5.0
    # Only one element, so no fill visible; test with larger
    cm2 = ContactMatrix(
        of = ["a", "b"], with = ["a", "b"], value = [1.0, 2.0], fill = 99.0
    )
    @test cm2[1, 1] == 1.0
    @test cm2[2, 2] == 2.0
    @test cm2[1, 2] == 99.0  # filled
    @test cm2[2, 1] == 99.0  # filled
end
