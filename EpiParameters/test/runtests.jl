using Test
using EpiParameters
using Distributions

@testset "EpiParameters" begin

    @testset "Database loads" begin
        diseases = list_diseases()
        @test length(diseases) > 20
        @test "COVID-19" in diseases
        @test "Ebola Virus Disease" in diseases

        params = list_parameters()
        @test "incubation period" in params
        @test "serial interval" in params
        @test "onset to death" in params
    end

    @testset "Query by disease" begin
        results = epiparameter(disease="COVID-19")
        @test length(results) > 0
        @test all(r -> occursin("COVID", r.disease), results)
    end

    @testset "Query by parameter type" begin
        results = epiparameter(epi_name="incubation period")
        @test length(results) > 10
        @test all(r -> r.epi_name == "incubation period", results)
    end

    @testset "Query combined filters" begin
        results = epiparameter(disease="COVID-19", epi_name="incubation period")
        @test length(results) > 0
        @test all(r -> occursin("COVID", r.disease), results)
        @test all(r -> r.epi_name == "incubation period", results)
    end

    @testset "Case-insensitive matching" begin
        r1 = epiparameter(disease="covid")
        r2 = epiparameter(disease="COVID")
        @test length(r1) == length(r2)
        @test length(r1) > 0
    end

    @testset "Author filter" begin
        results = epiparameter(author="Lessler")
        @test length(results) > 0
    end

    @testset "Distribution objects" begin
        # Ebola onset-to-death has a known Gamma distribution
        results = epiparameter(disease="Ebola", epi_name="onset to death")
        with_dist = filter(r -> !isnothing(r.distribution), results)
        @test length(with_dist) > 0

        d = with_dist[1].distribution
        @test d isa Gamma
        # Gamma(2.4, 3.3333) — Barry et al. 2018
        @test shape(d) ≈ 2.4 atol=0.1
        @test scale(d) ≈ 3.33 atol=0.1
    end

    @testset "Lognormal distribution" begin
        results = epiparameter(disease="COVID-19", epi_name="incubation period")
        with_dist = filter(r -> r.distribution isa LogNormal, results)
        @test length(with_dist) > 0

        d = with_dist[1].distribution
        @test d isa LogNormal
        @test params(d)[1] > 0  # meanlog
        @test params(d)[2] > 0  # sdlog
    end

    @testset "Negative binomial (offspring distribution)" begin
        results = epiparameter(epi_name="offspring distribution")
        with_dist = filter(r -> r.distribution isa NegativeBinomial, results)
        @test length(with_dist) > 0

        d = with_dist[1].distribution
        @test d isa NegativeBinomial
        @test mean(d) > 0
    end

    @testset "Offset" begin
        # Some Mpox entries have offset=4
        results = epiparameter(disease="Mpox")
        with_offset = filter(r -> r.offset > 0, results)
        @test length(with_offset) > 0
    end

    @testset "Summary statistics" begin
        results = epiparameter(disease="Adenovirus")
        @test length(results) > 0
        p = results[1]
        @test haskey(p.summary_stats, "median")
        @test p.summary_stats["median"] > 0
    end

    @testset "Citation metadata" begin
        results = epiparameter(disease="Ebola", epi_name="onset to death")
        p = results[1]
        @test haskey(p.citation, "year")
        @test haskey(p.citation, "doi")
        @test haskey(p.metadata, "units")
    end

    @testset "Display" begin
        results = epiparameter(disease="Ebola", epi_name="onset to death")
        with_dist = filter(r -> !isnothing(r.distribution), results)
        s = sprint(show, with_dist[1])
        @test occursin("Ebola", s)
        @test occursin("onset to death", s)
    end

    @testset "Empty query returns all" begin
        all_results = epiparameter()
        @test length(all_results) > 100
    end

    @testset "No matches returns empty" begin
        results = epiparameter(disease="Nonexistent Disease XYZ")
        @test isempty(results)
    end
end
