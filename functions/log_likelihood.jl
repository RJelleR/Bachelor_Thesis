
function log_likelihood(y, regressor, theta)
    omega = theta[1]
    phi = exp(theta[2]) / (1 + exp(theta[2]))
    alpha = exp(theta[3])
    sigma_2 = exp(theta[4])

    n = length(y)
    beta = zeros(n)
    beta[1] = omega / (1 - phi)

    for t in 2:n

        beta[t] = omega + phi .* beta[t-1] + alpha .* (y[t-1] - beta[t-1] .* regressor[t-1]) .* regressor[t-1]

    end
    l = .-0.5 * log(sigma_2) .- 0.5 * (y .- beta .* regressor) .^ 2 / sigma_2
    llik = mean(l)

    return llik
end