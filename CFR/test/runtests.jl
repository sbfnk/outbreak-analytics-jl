using Test
using CFR
using DataFrames
using Distributions
using Dates
using Random

@testset "CFR" begin

    # Helper: simple epidemic data
    function make_epidemic(; n=50, true_cfr=0.1, seed=42)
        rng = Xoshiro(seed)
        dates = Date(2024, 1, 1) .+ Day.(0:(n - 1))
        # Bell-shaped epidemic curve
        peak = n ÷ 2
        cases = [max(1, round(Int, 100 * exp(-(t - peak)^2 / (2 * (n / 5)^2)))) for t in 1:n]
        deaths = [round(Int, c * true_cfr) for c in cases]
        DataFrame(date=dates, cases=cases, deaths=deaths)
    end

    @testset "CfrEstimate display" begin
        est = CfrEstimate(0.05, 0.03, 0.08)
        s = sprint(show, est)
        @test occursin("0.05", s)
        @test occursin("0.03", s)
        @test occursin("0.08", s)
    end

    @testset "Naive CFR" begin
        data = DataFrame(
            date=Date(2024, 1, 1) .+ Day.(0:4),
            cases=[10, 20, 30, 20, 10],
            deaths=[1, 2, 3, 2, 1],
        )
        result = cfr_static(data)
        @test result.estimate ≈ 9 / 90 atol=1e-10
        @test result.lower < result.estimate
        @test result.upper > result.estimate
        @test 0 < result.lower
        @test result.upper < 1
    end

    @testset "Delay with point mass at 0 equals naive" begin
        data = DataFrame(
            date=Date(2024, 1, 1) .+ Day.(0:9),
            cases=[10, 20, 30, 40, 50, 40, 30, 20, 10, 5],
            deaths=[1, 2, 3, 4, 5, 4, 3, 2, 1, 0],
        )
        naive = cfr_static(data)
        # A distribution concentrated very near 0
        adjusted = cfr_static(data; delay_density=Exponential(0.001))
        @test adjusted.estimate ≈ naive.estimate atol=0.01
    end

    @testset "Delay-adjusted ≥ naive" begin
        data = make_epidemic(n=50, true_cfr=0.1)
        delay = Gamma(2.4, 3.33)
        naive = cfr_static(data)
        adjusted = cfr_static(data; delay_density=delay)
        @test adjusted.estimate ≥ naive.estimate - 1e-10
    end

    @testset "estimate_outcomes" begin
        data = DataFrame(
            date=Date(2024, 1, 1) .+ Day.(0:4),
            cases=[100, 0, 0, 0, 0],
            deaths=[0, 0, 0, 0, 10],
        )
        delay = Gamma(2.0, 1.5)
        result = estimate_outcomes(data, delay)
        @test hasproperty(result, :known_outcomes)
        @test length(result.known_outcomes) == 5
        # All known outcomes should be non-negative
        @test all(result.known_outcomes .≥ 0)
        # Total known outcomes ≤ total cases
        @test sum(result.known_outcomes) ≤ sum(data.cases) + 1e-10
    end

    @testset "Time-varying CFR" begin
        data = make_epidemic(n=30, true_cfr=0.1)

        # Naive time-varying
        result = cfr_time_varying(data; burn_in=7)
        @test result isa DataFrame
        @test :date in propertynames(result)
        @test :cfr in propertynames(result)
        @test :lower in propertynames(result)
        @test :upper in propertynames(result)
        @test nrow(result) == 30 - 7

        # Delay-adjusted time-varying
        delay = Gamma(2.4, 3.33)
        result_adj = cfr_time_varying(data; delay_density=delay, burn_in=7)
        @test :known_outcomes in propertynames(result_adj)
        @test nrow(result_adj) == 30 - 7
    end

    @testset "Smoothing window" begin
        data = make_epidemic(n=30, true_cfr=0.1)
        result = cfr_time_varying(data; burn_in=7, smoothing_window=3)
        @test nrow(result) == 30 - 7
    end

    @testset "Edge cases" begin
        # Zero deaths
        data = DataFrame(
            date=Date(2024, 1, 1) .+ Day.(0:4),
            cases=[10, 20, 30, 20, 10],
            deaths=[0, 0, 0, 0, 0],
        )
        result = cfr_static(data)
        @test result.estimate == 0.0
        @test result.lower == 0.0

        # Zero cases
        data_zero = DataFrame(
            date=Date(2024, 1, 1) .+ Day.(0:4),
            cases=[0, 0, 0, 0, 0],
            deaths=[0, 0, 0, 0, 0],
        )
        result_zero = cfr_static(data_zero)
        @test isnan(result_zero.estimate)
    end

    @testset "Input validation" begin
        bad_data = DataFrame(x=[1, 2, 3])
        @test_throws ArgumentError cfr_static(bad_data)
        @test_throws ArgumentError cfr_time_varying(bad_data)
    end

    @testset "Confidence level" begin
        data = make_epidemic(n=50, true_cfr=0.1)
        ci90 = cfr_static(data; confidence_level=0.90)
        ci95 = cfr_static(data; confidence_level=0.95)
        # 90% CI should be narrower than 95% CI
        @test (ci95.upper - ci95.lower) ≥ (ci90.upper - ci90.lower)
    end
end
