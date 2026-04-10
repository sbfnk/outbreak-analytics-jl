@testset "Age groups" begin
    @testset "reduce_agegroups" begin
        @test reduce_agegroups([3, 7, 12, 20], [0, 5, 15]) == [0, 5, 5, 15]
        @test reduce_agegroups([0, 1, 4, 5, 14, 15, 99], [0, 5, 15]) == [0, 0, 0, 5, 5, 15, 15]
        # Missing values
        result = reduce_agegroups([missing, 5], [0, 5])
        @test ismissing(result[1])
        @test result[2] == 5
    end

    @testset "limits_to_agegroups" begin
        @test limits_to_agegroups([0, 1, 5, 15]) == ["[0,1)", "[1,5)", "[5,15)", "15+"]
        @test limits_to_agegroups([0, 5, 10]) == ["[0,5)", "[5,10)", "10+"]
        @test limits_to_agegroups([0]) == ["0+"]
    end

    @testset "assign_age_groups" begin
        survey = polymod()

        # Basic assignment
        grouped = survey |>
            s -> filter_survey(s; countries = ["United Kingdom"]) |>
            s -> assign_age_groups(s; age_limits = [0, 1, 5, 15])

        @test hasproperty(grouped.participants, :age_group)
        @test hasproperty(grouped.participants, :lower_age_limit)
        @test hasproperty(grouped.contacts, :contact_age_group)

        # Check age group labels
        age_groups = sort(unique(skipmissing(grouped.participants.age_group)))
        @test Set(age_groups) == Set(["[0,1)", "[1,5)", "[5,15)", "15+"])

        # All participants should have valid age groups
        @test !any(ismissing, grouped.participants.age_group)

        # Lower age limits should correspond to groups
        for row in eachrow(grouped.participants)
            if row.age_group == "[0,1)"
                @test row.lower_age_limit == 0
            elseif row.age_group == "[1,5)"
                @test row.lower_age_limit == 1
            elseif row.age_group == "[5,15)"
                @test row.lower_age_limit == 5
            elseif row.age_group == "15+"
                @test row.lower_age_limit == 15
            end
        end
    end

    @testset "assign_age_groups with inferred limits" begin
        # Small synthetic survey
        parts = DataFrame(
            part_id = [1, 2, 3],
            part_age_exact = [5, 12, 25],
        )
        conts = DataFrame(
            part_id = [1, 1, 2, 3],
            cnt_age_exact = [3, 10, 20, 30],
        )
        survey = ContactSurvey(parts, conts)
        grouped = assign_age_groups(survey)

        @test hasproperty(grouped.participants, :age_group)
        @test nrow(grouped.participants) == 3
    end

    @testset "age imputation from ranges" begin
        parts = DataFrame(
            part_id = [1, 2],
            part_age_exact = [missing, missing],
            part_age_est_min = [10, 20],
            part_age_est_max = [14, 30],
        )
        conts = DataFrame(
            part_id = [1, 2],
            cnt_age_exact = [5, 15],
        )
        survey = ContactSurvey(parts, conts)

        # Mean imputation
        grouped = assign_age_groups(survey; age_limits = [0, 10, 20],
            estimated_participant_age = "mean")
        @test grouped.participants.part_age[1] == 12  # floor((10+14)/2)
        @test grouped.participants.part_age[2] == 25  # floor((20+30)/2)
    end

    @testset "missing_contact_age handling" begin
        parts = DataFrame(
            part_id = [1, 2, 3],
            part_age_exact = [5, 12, 25],
        )
        conts = DataFrame(
            part_id = [1, 1, 2, 3],
            cnt_age_exact = [3, missing, 20, 30],
        )

        # "remove" removes participant 1 (has a contact with missing age)
        survey = ContactSurvey(parts, conts)
        grouped = assign_age_groups(survey; age_limits = [0, 10, 20],
            missing_contact_age = "remove")
        @test !(1 in grouped.participants.part_id)
        @test Set(grouped.participants.part_id) == Set([2, 3])

        # "ignore" keeps participant 1 but drops the missing contact
        survey = ContactSurvey(copy(parts), copy(conts))
        grouped = assign_age_groups(survey; age_limits = [0, 10, 20],
            missing_contact_age = "ignore")
        @test 1 in grouped.participants.part_id
        @test !any(ismissing, grouped.contacts.cnt_age)
    end
end
