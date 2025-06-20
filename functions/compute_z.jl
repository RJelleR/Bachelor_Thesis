using .FractionalMoments
include("fractional_moments.jl")

function compute_z(returns, gamma; H=4, n_sim=100_000, window=520, target_length=nothing)
    n = length(returns)
    max_index = n - H

    if isnothing(target_length)
        range = window:max_index
    else
        range = (max_index-target_length+1):max_index
    end

    z = Float64[]
    for t in range
        r_window = returns[(t-window+1):t]
        garch_model = fit(GARCH{1,1}, r_window)
        omega, alpha, beta = coef(garch_model)

        R_H_2 = simulate_conditional_variance(omega, alpha, beta, H, n_sim)
        distribution_R = kde(R_H_2)
        pdf_R = FractionalMoments.KDEDistribution(distribution_R)
        moment = FractionalMoments.fractional_moment(pdf_R, gamma; use_abs=true)
        push!(z, moment)
    end
    return z
end