function Hess_fun_beta(par, y, x)
    n = length(x)
    omega = par[1]
    phi = par[2]
    alpha = par[3]
    sigma_2 = par[4]

    n = length(y)
    beta = zeros(eltype(par), n)
    beta[1] = omega / (1 - phi)

    for t in 2:n

        beta[t] = omega + phi .* beta[t-1] + alpha .* (y[t-1] - beta[t-1] .* x[t-1]) .* x[t-1]

    end
    l = .-0.5 * log(sigma_2) .- 0.5 * (y .- beta .* x) .^ 2 / sigma_2
    llik = mean(l)

    return llik
end