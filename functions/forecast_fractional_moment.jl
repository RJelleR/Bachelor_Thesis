#Not used in any code/application
function forecast_power_fractional_moments(returns)
    Random.seed!(42)
    training_end = 1000
    n = length(returns)

    window = 500
    gammas = [0.5, 1, 1.5, 2]
    forecast_results = Dict(gamma => Float64[] for gamma in gammas)
    realized_volatility = Float64[]

    for t in (training_end-window+1):n-5
        r_window = returns[(t-window+1):t]

        garch_model = fit(GARCH{1,1}, r_window)
        omega, alpha, beta = coef(garch_model)


        R_t_1 = simulate_conditional_return(omega, alpha, beta, 1, 100_000)

        for gamma in gammas
            forecast = mean(abs.(R_t_1) .^ gamma)
            push!(forecast_results[gamma], forecast)
        end

        push!(realized_volatility, mean(returns[t+1:t+5] .^ 2))
    end
    return forecast_results, realized_volatility
end
