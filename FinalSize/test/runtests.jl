using Test
using FinalSize
using ContactMatrices
using DataFrames
using LinearAlgebra

@testset "FinalSize" begin

    @testset "Homogeneous final size" begin
        # R₀ = 2 → φ ≈ 0.7968
        φ = final_size(2.0)
        @test φ ≈ 0.7968 atol=1e-4
        @test φ ≈ 1 - exp(-2.0 * φ) atol=1e-12  # self-consistency

        # R₀ = 1.5 → φ ≈ 0.5828
        φ15 = final_size(1.5)
        @test φ15 ≈ 0.5828 atol=1e-4
        @test φ15 ≈ 1 - exp(-1.5 * φ15) atol=1e-12

        # R₀ ≤ 1 → no epidemic
        @test final_size(0.5) == 0.0
        @test final_size(1.0) == 0.0
        @test final_size(0.0) == 0.0

        # R₀ < 0 → error
        @test_throws DomainError final_size(-1.0)
    end

    @testset "Heterogeneous — raw matrix" begin
        # 2 groups with symmetric contact matrix
        C = [8.0 2.0; 2.0 4.0]
        d = [0.6, 0.4]

        result = final_size(1.5, C; demography=d)
        @test result isa DataFrame
        @test names(result) == ["group", "susc_group", "susceptibility", "p_infected"]
        @test nrow(result) == 2
        @test result.group == ["1", "2"]

        # All infected proportions should be between 0 and 1
        @test all(0 .< result.p_infected .< 1)

        # Group 1 (more contacts) should have higher attack rate
        @test result.p_infected[1] > result.p_infected[2]
    end

    @testset "Heterogeneous — single group reduces to homogeneous" begin
        # Single group with contact rate C, demography [1.0]
        # Should give same result as homogeneous with same R₀
        C = reshape([1.0], 1, 1)
        result = final_size(2.0, C; demography=[1.0])
        @test result.p_infected[1] ≈ final_size(2.0) atol=1e-6
    end

    @testset "ContactMatrix integration" begin
        labels = ["[0,20)", "20+"]
        C = [8.0 3.0; 3.0 5.0]
        cm = ContactMatrix(C, labels)
        d = [0.3, 0.7]

        result = final_size(1.5, cm; demography=d)
        @test result isa DataFrame
        @test result.group == sort(["[0,20)", "20+"])
        @test all(0 .< result.p_infected .< 1)
    end

    @testset "Susceptibility groups" begin
        C = [8.0 2.0; 2.0 4.0]
        d = [0.6, 0.4]

        # Two susceptibility classes: fully susceptible and vaccinated (50% reduction)
        s = [1.0 0.5; 1.0 0.5]
        p = [0.5 0.5; 0.5 0.5]

        result = final_size(1.5, C; demography=d,
                            susceptibility=s, p_susceptibility=p)
        @test nrow(result) == 4  # 2 groups × 2 susc classes
        @test result.susc_group == [1, 2, 1, 2]

        # Vaccinated class should have lower attack rate
        g1_full = result.p_infected[result.group .== "1" .&& result.susc_group .== 1][1]
        g1_vacc = result.p_infected[result.group .== "1" .&& result.susc_group .== 2][1]
        @test g1_full > g1_vacc
    end

    @testset "Edge cases" begin
        C = [5.0 1.0; 1.0 3.0]
        d = [0.5, 0.5]

        # R₀ = 0 → no epidemic
        result = final_size(0.0, C; demography=d)
        @test all(result.p_infected .== 0.0)

        # Demography normalisation
        result1 = final_size(1.5, C; demography=[50, 50])
        result2 = final_size(1.5, C; demography=[0.5, 0.5])
        @test result1.p_infected ≈ result2.p_infected atol=1e-10
    end

    @testset "Dimension mismatch errors" begin
        C = [5.0 1.0; 1.0 3.0]
        @test_throws DimensionMismatch final_size(1.5, C; demography=[0.5])
        @test_throws DimensionMismatch final_size(1.5, C; demography=[0.5, 0.5],
                                                  susceptibility=ones(3, 1))
    end
end
