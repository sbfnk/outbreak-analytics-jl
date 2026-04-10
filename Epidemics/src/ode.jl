"""
Age-structured SEIR ODE system.

State vector layout: [S_1..S_n, E_1..E_n, I_1..I_n, R_1..R_n]
"""
function seir_ode!(du, u, p, t)
    (; beta, sigma, gamma, C, N, interventions, vaccinations, n) = p

    S = @view u[1:n]
    E = @view u[n+1:2n]
    I = @view u[2n+1:3n]

    dS = @view du[1:n]
    dE = @view du[n+1:2n]
    dI = @view du[2n+1:3n]
    dR = @view du[3n+1:4n]

    # Effective beta with interventions
    beta_eff = beta
    for iv in interventions
        if iv.time_begin <= t <= iv.time_end
            beta_eff *= (1 - iv.reduction)
        end
    end

    # Force of infection and compartment flows
    for i in 1:n
        λ_i = zero(eltype(u))
        for j in 1:n
            λ_i += C[i, j] * I[j] / N[j]
        end
        λ_i *= beta_eff

        new_infections = λ_i * S[i]
        progression = sigma * E[i]
        recovery = gamma * I[i]

        dS[i] = -new_infections
        dE[i] = new_infections - progression
        dI[i] = progression - recovery
        dR[i] = recovery
    end

    # Vaccination: S → R
    for vacc in vaccinations
        if vacc.time_begin <= t <= vacc.time_end
            for i in 1:n
                if isnothing(vacc.groups) || vacc.groups[i]
                    vacc_flow = vacc.rate * S[i]
                    dS[i] -= vacc_flow
                    dR[i] += vacc_flow
                end
            end
        end
    end

    nothing
end
