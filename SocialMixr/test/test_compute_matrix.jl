@testset "Compute matrix" begin
    survey = polymod()

    @testset "basic UK matrix" begin
        result = survey |>
            s -> filter_survey(s; countries = ["United Kingdom"]) |>
            s -> assign_age_groups(s; age_limits = [0, 1, 5, 15]) |>
            compute_matrix

        cm = result.matrix
        @test cm isa ContactMatrix

        labels = groupings(cm).labels[1]
        @test length(labels) == 4
        @test Set(labels) == Set(["[0,1)", "[1,5)", "[5,15)", "15+"])

        # Matrix should have reasonable values (all positive means)
        m = Matrix(cm)
        @test all(m .>= 0)

        # Participant counts
        @test nrow(result.participants) == 4
        @test sum(result.participants.participants) == 1011
        @test all(result.participants.proportion .> 0)
        @test sum(result.participants.proportion) ≈ 1.0
    end

    @testset "counts mode" begin
        result_mean = survey |>
            s -> filter_survey(s; countries = ["United Kingdom"]) |>
            s -> assign_age_groups(s; age_limits = [0, 1, 5, 15]) |>
            compute_matrix

        result_count = survey |>
            s -> filter_survey(s; countries = ["United Kingdom"]) |>
            s -> assign_age_groups(s; age_limits = [0, 1, 5, 15]) |>
            s -> compute_matrix(s; counts = true)

        # Counts should be >= means (for groups with >1 participant)
        m_mean = Matrix(result_mean.matrix)
        m_count = Matrix(result_count.matrix)
        @test all(m_count .>= m_mean .- 1e-10)
    end

    @testset "filtering" begin
        # School contacts only
        school_result = survey |>
            s -> filter_survey(s; filter = Dict(:cnt_school => 1)) |>
            s -> assign_age_groups(s; age_limits = [0, 20, 60]) |>
            compute_matrix

        cm = school_result.matrix
        m = Matrix(cm)
        @test all(m .>= 0)

        # All-contacts matrix should have higher values
        all_result = survey |>
            s -> assign_age_groups(s; age_limits = [0, 20, 60]) |>
            compute_matrix

        m_all = Matrix(all_result.matrix)
        # School contacts should be less than or equal to all contacts
        # (using same label order)
        @test sum(m) <= sum(m_all) + 1e-10
    end

    @testset "country filtering" begin
        uk_result = survey |>
            s -> filter_survey(s; countries = ["United Kingdom"]) |>
            s -> assign_age_groups(s; age_limits = [0, 5, 15]) |>
            compute_matrix

        de_result = survey |>
            s -> filter_survey(s; countries = ["Germany"]) |>
            s -> assign_age_groups(s; age_limits = [0, 5, 15]) |>
            compute_matrix

        # Different countries should give different matrices
        @test Matrix(uk_result.matrix) != Matrix(de_result.matrix)
    end

    @testset "error on missing age groups" begin
        @test_throws ArgumentError compute_matrix(survey)
    end
end
