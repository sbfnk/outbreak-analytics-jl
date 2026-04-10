@testset "Weights" begin
    survey = polymod()

    @testset "day-of-week weighting" begin
        uk = filter_survey(survey; countries = ["United Kingdom"])
        grouped = assign_age_groups(uk; age_limits = [0, 18, 60])

        # Apply day-of-week weighting using grouped method
        # POLYMOD: 0=Sunday, 1=Monday, ..., 6=Saturday
        weighted = weigh(grouped, "dayofweek"; target = [5, 2], groups = [1:5, [0, 6]])

        @test hasproperty(weighted.participants, :weight)
        @test all(!ismissing, weighted.participants.weight)

        # Weights should differ between weekday and weekend participants
        weekday_w = [row.weight for row in eachrow(weighted.participants) if
                     !ismissing(row.dayofweek) && row.dayofweek in 1:5]
        weekend_w = [row.weight for row in eachrow(weighted.participants) if
                     !ismissing(row.dayofweek) && row.dayofweek in [0, 6]]

        @test length(weekday_w) > 0
        @test length(weekend_w) > 0
        # Weekday and weekend weights should be different
        @test weekday_w[1] != weekend_w[1]
    end

    @testset "direct weighting" begin
        parts = DataFrame(
            part_id = [1, 2, 3],
            part_age_exact = [5, 12, 25],
            hh_size = [2, 4, 1],
        )
        conts = DataFrame(
            part_id = [1, 2, 3],
            cnt_age_exact = [3, 10, 20],
        )
        survey_small = ContactSurvey(parts, conts)
        grouped = assign_age_groups(survey_small; age_limits = [0, 10, 20])

        weighted = weigh(grouped, "hh_size")
        @test weighted.participants.weight[1] == 2.0
        @test weighted.participants.weight[2] == 4.0
        @test weighted.participants.weight[3] == 1.0
    end

    @testset "weight composition" begin
        # Multiple weigh calls should compose (multiplicative)
        parts = DataFrame(
            part_id = [1, 2, 3, 4],
            part_age_exact = [5, 12, 25, 30],
            dayofweek = [1, 2, 6, 0],  # 1,2=weekday, 6,0=weekend
            hh_size = [2, 3, 1, 4],
        )
        conts = DataFrame(
            part_id = [1, 2, 3, 4],
            cnt_age_exact = [3, 10, 20, 15],
        )
        survey_small = ContactSurvey(parts, conts)
        grouped = assign_age_groups(survey_small; age_limits = [0, 10, 20])

        # First weigh by dayofweek, then by hh_size
        w1 = weigh(grouped, "dayofweek"; target = [5, 2], groups = [1:5, [0, 6]])
        w2 = weigh(w1, "hh_size")

        # Weights should be product of both factors
        @test all(!ismissing, w2.participants.weight)
        # hh_size weights are multiplicative on top of dayofweek weights
        for i in 1:nrow(w2.participants)
            expected = w1.participants.weight[i] * w2.participants.hh_size[i]
            @test w2.participants.weight[i] ≈ expected
        end
    end

    @testset "normalise_weights!" begin
        parts = DataFrame(
            age_group = ["A", "A", "B", "B"],
            weight = [2.0, 4.0, 1.0, 3.0],
        )
        SocialMixr.normalise_weights!(parts)
        # Within group A: sum should be 2 (2 participants)
        a_weights = parts.weight[parts.age_group .== "A"]
        @test sum(a_weights) ≈ 2.0
        # Within group B: sum should be 2 (2 participants)
        b_weights = parts.weight[parts.age_group .== "B"]
        @test sum(b_weights) ≈ 2.0
    end

    @testset "weight threshold" begin
        parts = DataFrame(
            age_group = ["A", "A", "A", "A"],
            weight = [0.1, 0.1, 0.1, 10.0],
        )
        SocialMixr.normalise_weights!(parts; threshold = 3.0)
        # Sum should still equal group size after normalisation
        @test sum(parts.weight) ≈ 4.0
        # The extreme weight should have been reduced compared to no threshold
        parts2 = DataFrame(
            age_group = ["A", "A", "A", "A"],
            weight = [0.1, 0.1, 0.1, 10.0],
        )
        SocialMixr.normalise_weights!(parts2)
        @test maximum(parts.weight) < maximum(parts2.weight)
    end
end
