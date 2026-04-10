@testset "Post-processing" begin
    survey = polymod()
    uk_pop = SocialMixr.load_population(
        joinpath(pkgdir(SocialMixr), "data", "wpp_uk.csv"),
    )

    result = survey |>
        s -> filter_survey(s; countries = ["United Kingdom"]) |>
        s -> assign_age_groups(s; age_limits = [0, 1, 5, 15]) |>
        compute_matrix

    @testset "pop_age" begin
        # Basic aggregation
        aggregated = pop_age(uk_pop, [0, 15, 60])
        @test nrow(aggregated) == 3
        @test aggregated.lower_age_limit == [0, 15, 60]
        @test all(aggregated.population .> 0)

        # Total population should be approximately preserved
        @test sum(aggregated.population) ≈ sum(uk_pop.population) atol = 1
    end

    @testset "pop_age interpolation" begin
        # Age limit 18 is not in the 5-year bands, needs interpolation
        interp = pop_age(uk_pop, [0, 18, 65])
        @test nrow(interp) == 3
        @test interp.lower_age_limit == [0, 18, 65]
        @test all(interp.population .> 0)
        # Total should approximately match
        @test sum(interp.population) ≈ sum(uk_pop.population) atol = sum(uk_pop.population) * 0.01
    end

    @testset "symmetrise" begin
        sym = symmetrise(result, uk_pop)
        cm = sym.matrix

        @test cm isa ContactMatrix
        m = Matrix(cm)

        # Symmetry check: c_ij * N_i ≈ c_ji * N_j
        pop = SocialMixr._population_for_matrix(cm, uk_pop)
        for i in axes(m, 1)
            for j in axes(m, 2)
                @test m[i, j] * pop[i] ≈ m[j, i] * pop[j] atol = 1e-6
            end
        end
    end

    @testset "split_matrix" begin
        sp = split_matrix(result, uk_pop)

        # Expected values from R vignette
        @test sp.mean_contacts ≈ 11.55 atol = 0.1
        @test sp.normalisation ≈ 1.039 atol = 0.01

        # Contacts should be a vector of length 4
        @test length(sp.contacts) == 4
        @test all(sp.contacts .> 0)

        # Assortativity matrix should exist
        @test sp.matrix isa ContactMatrix
    end

    @testset "per_capita" begin
        de_pop = polymod_population(; countries = ["Germany"])
        de_result = survey |>
            s -> filter_survey(s; countries = ["Germany"]) |>
            s -> assign_age_groups(s; age_limits = [0, 60]) |>
            compute_matrix

        # Symmetric + per_capita should give a fully symmetric matrix
        de_sym = symmetrise(de_result, de_pop)
        de_pc = per_capita(de_sym, de_pop)

        m = Matrix(de_pc.matrix)
        @test m[1, 2] ≈ m[2, 1] atol = 1e-15

        # Per capita rates should be small (contacts per person per individual)
        @test all(m .< 1e-3)
        @test all(m .> 0)
    end

    @testset "pipeline composition" begin
        # Full pipeline: filter → age groups → weigh → compute → symmetrise
        pop = polymod_population(; countries = ["United Kingdom"])
        full_result = survey |>
            s -> filter_survey(s; countries = ["United Kingdom"]) |>
            s -> assign_age_groups(s; age_limits = [0, 18, 60]) |>
            s -> weigh(s, "dayofweek"; target = [5, 2], groups = [1:5, [0, 6]]) |>
            compute_matrix |>
            r -> symmetrise(r, pop)

        @test full_result.matrix isa ContactMatrix
        m = Matrix(full_result.matrix)
        @test size(m) == (3, 3)
        @test all(m .>= 0)

        # Symmetry check
        pop_vec = SocialMixr._population_for_matrix(full_result.matrix, pop)
        for i in 1:3, j in 1:3
            @test m[i, j] * pop_vec[i] ≈ m[j, i] * pop_vec[j] atol = 1e-8
        end
    end
end
