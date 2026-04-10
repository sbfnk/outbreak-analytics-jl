using ContactMatrices
using Test

@testset "make_symmetric" begin
    # Labels sorted: old=1, young=2
    # So we specify the matrix in sorted order directly
    cm = ContactMatrix([3.0 2.0; 1.0 4.0], ["old", "young"])
    pop = [500.0, 1000.0]  # old=500, young=1000
    cm_sym = make_symmetric(cm, pop)

    # Verify symmetry: c_ij * N_i == c_ji * N_j
    for i in 1:2, j in 1:2
        @test cm_sym[i, j] * pop[i] ≈ cm_sym[j, i] * pop[j]
    end

    # Check formula: c_sym[i,j] = (N_i * c[i,j] + N_j * c[j,i]) / (2 * N_i)
    @test cm_sym[1, 1] ≈ 3.0   # diagonal unchanged
    @test cm_sym[2, 2] ≈ 4.0
    @test cm_sym[1, 2] ≈ (500.0 * 2.0 + 1000.0 * 1.0) / (2.0 * 500.0)
    @test cm_sym[2, 1] ≈ (1000.0 * 1.0 + 500.0 * 2.0) / (2.0 * 1000.0)

    @test setting(cm_sym) == :all
    @test groupings(cm_sym) == groupings(cm)
end

@testset "make_symmetric population validation" begin
    cm = ContactMatrix([1.0 2.0; 3.0 4.0], ["a", "b"])
    @test_throws ArgumentError make_symmetric(cm, [100.0])  # wrong length
end
