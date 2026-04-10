"""Convert ODE solution to a tidy DataFrame."""
function _solution_to_dataframe(sol, n, labels)
    rows = NamedTuple{(:time, :group, :S, :E, :I, :R),
                       Tuple{Float64, String, Float64, Float64, Float64, Float64}}[]
    sizehint!(rows, length(sol.t) * n)

    for (t_idx, t) in enumerate(sol.t)
        u = sol.u[t_idx]
        for i in 1:n
            push!(rows, (
                time  = t,
                group = labels[i],
                S     = max(u[i], 0.0),
                E     = max(u[n + i], 0.0),
                I     = max(u[2n + i], 0.0),
                R     = max(u[3n + i], 0.0),
            ))
        end
    end

    DataFrame(rows)
end
