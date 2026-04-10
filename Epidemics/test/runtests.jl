using Test
using Epidemics
using ContactMatrices
using DataFrames

@testset "Epidemics" begin

    @testset "Homogeneous SEIR" begin
        # Single group, should produce an epidemic curve
        model = SEIR(beta=0.3, sigma=1/5, gamma=1/10)
        C = reshape([10.0], 1, 1)
        result = simulate(model, C;
                          demography=[1e6],
                          initial_infected=[10.0],
                          tspan=(0.0, 300.0))

        @test result isa DataFrame
        @test :time in propertynames(result)
        @test :S in propertynames(result)
        @test :E in propertynames(result)
        @test :I in propertynames(result)
        @test :R in propertynames(result)

        # Epidemic should have run its course — most people recovered
        final = result[result.time .== maximum(result.time), :]
        @test final.R[1] > 0.5e6  # significant epidemic
        @test final.I[1] < 100    # epidemic is over
    end

    @testset "Conservation of population" begin
        model = SEIR(beta=0.3, sigma=1/5, gamma=1/10)
        C = [8.0 2.0; 2.0 4.0]
        N = [5e5, 5e5]

        result = simulate(model, C;
                          demography=N,
                          initial_infected=[10.0, 0.0],
                          tspan=(0.0, 200.0))

        # S + E + I + R should equal N for each group at each time
        for g in unique(result.group)
            rows = result[result.group .== g, :]
            totals = rows.S .+ rows.E .+ rows.I .+ rows.R
            expected = g == "1" ? N[1] : N[2]
            @test all(abs.(totals .- expected) .< 1.0)  # within 1 person
        end
    end

    @testset "No transmission (beta=0)" begin
        model = SEIR(beta=0.0, sigma=1/5, gamma=1/10)
        C = [10.0 5.0; 5.0 8.0]

        result = simulate(model, C;
                          demography=[1e6, 1e6],
                          initial_infected=[100.0, 0.0],
                          tspan=(0.0, 100.0))

        # Only the initial infected should move through E→I→R
        final = result[result.time .== maximum(result.time), :]
        g1 = final[final.group .== "1", :]
        g2 = final[final.group .== "2", :]

        @test g1.S[1] ≈ 1e6 - 100 atol=1.0
        @test g2.S[1] ≈ 1e6 atol=1.0
        @test g2.R[1] ≈ 0.0 atol=1.0
    end

    @testset "ContactMatrix dispatch" begin
        model = SEIR(beta=0.3, sigma=1/5, gamma=1/10)
        labels = ["young", "old"]
        C_raw = [8.0 2.0; 2.0 4.0]
        cm = ContactMatrix(C_raw, labels)
        N = [5e5, 5e5]

        result = simulate(model, cm;
                          demography=N,
                          initial_infected=[10.0, 0.0],
                          tspan=(0.0, 200.0))

        @test "old" in result.group
        @test "young" in result.group

        # Compare with raw matrix (labels will differ but dynamics should match)
        result_raw = simulate(model, C_raw;
                              demography=N,
                              initial_infected=[10.0, 0.0],
                              tspan=(0.0, 200.0))

        # Note: ContactMatrix sorts labels alphabetically, so "old" comes first
        # in the ContactMatrix version. The raw matrix keeps original order.
        # Dynamics should still match when we compare by the same group ordering.
        final_cm = sort(result[result.time .== maximum(result.time), :], :group)
        final_raw = sort(result_raw[result_raw.time .== maximum(result_raw.time), :], :group)
        @test final_cm.R ≈ final_raw.R atol=1.0
    end

    @testset "Intervention reduces peak" begin
        model = SEIR(beta=0.3, sigma=1/5, gamma=1/10)
        C = reshape([10.0], 1, 1)
        N = [1e6]
        I0 = [10.0]

        # No intervention
        base = simulate(model, C; demography=N, initial_infected=I0,
                        tspan=(0.0, 300.0))
        peak_base = maximum(base.I)

        # Strong intervention during growth phase
        iv = Intervention(time_begin=10.0, time_end=100.0, reduction=0.7)
        with_iv = simulate(model, C; demography=N, initial_infected=I0,
                           tspan=(0.0, 300.0), interventions=[iv])
        peak_iv = maximum(with_iv.I)

        @test peak_iv < peak_base
    end

    @testset "Vaccination reduces epidemic" begin
        model = SEIR(beta=0.3, sigma=1/5, gamma=1/10)
        C = reshape([10.0], 1, 1)
        N = [1e6]
        I0 = [10.0]

        # No vaccination
        base = simulate(model, C; demography=N, initial_infected=I0,
                        tspan=(0.0, 300.0))

        # Early vaccination campaign
        vacc = Vaccination(time_begin=0.0, time_end=50.0, rate=0.02)
        with_vacc = simulate(model, C; demography=N, initial_infected=I0,
                             tspan=(0.0, 300.0), vaccination=[vacc])

        final_base = base[base.time .== maximum(base.time), :].R[1]
        final_vacc = with_vacc[with_vacc.time .== maximum(with_vacc.time), :].R[1]

        # With vaccination, fewer people get infected (more go S→R directly)
        # but total R should be similar or higher. Peak I should be lower.
        peak_base = maximum(base.I)
        peak_vacc = maximum(with_vacc.I)
        @test peak_vacc < peak_base
    end

    @testset "Targeted vaccination" begin
        model = SEIR(beta=0.3, sigma=1/5, gamma=1/10)
        C = [8.0 2.0; 2.0 4.0]
        N = [5e5, 5e5]

        # Vaccinate only group 2
        vacc = Vaccination(time_begin=0.0, time_end=30.0, rate=0.05,
                           groups=[false, true])
        result = simulate(model, C; demography=N,
                          initial_infected=[10.0, 0.0],
                          tspan=(0.0, 200.0), vaccination=[vacc])

        # Group 2 should have fewer infections at peak
        g1_peak = maximum(result[result.group .== "1", :].I)
        g2_peak = maximum(result[result.group .== "2", :].I)
        @test g2_peak < g1_peak
    end

    @testset "Dimension mismatches" begin
        model = SEIR(beta=0.3, sigma=1/5, gamma=1/10)
        C = [8.0 2.0; 2.0 4.0]

        @test_throws DimensionMismatch simulate(model, C;
            demography=[1e6], initial_infected=[10.0], tspan=(0.0, 100.0))

        @test_throws DimensionMismatch simulate(model, C;
            demography=[5e5, 5e5], initial_infected=[10.0], tspan=(0.0, 100.0))
    end

    @testset "Parameter validation" begin
        @test_throws DomainError SEIR(beta=-1.0, sigma=1/5, gamma=1/10)
        @test_throws DomainError SEIR(beta=0.3, sigma=0.0, gamma=1/10)
        @test_throws ArgumentError Intervention(time_begin=10.0, time_end=5.0, reduction=0.5)
        @test_throws DomainError Intervention(time_begin=0.0, time_end=10.0, reduction=1.5)
        @test_throws DomainError Vaccination(time_begin=0.0, rate=-0.01)
    end
end
