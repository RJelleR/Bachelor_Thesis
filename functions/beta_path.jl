function beta_path(regressor, y, parameters)
    omega_hat = parameters[1]
    phi_hat = parameters[2]
    alpha_hat = parameters[3]

    n = length(regressor)
    beta = zeros(n)
    beta[1] = omega_hat / (1 - phi_hat)

    for t in 2:n
        beta[t] = omega_hat + phi_hat .* beta[t-1] + alpha_hat .* (y[t-1] .- beta[t-1] .* regressor[t-1]) .* regressor[t-1]
    end
    return beta
end