function simulate_conditional_return(omega, alpha, beta, H, n_sim=100_000)
    R_H = zeros(n_sim)
    sigma_2_0 = omega / (1 - alpha - beta)

    for i in 1:n_sim
        r = 0.0
        sigma_2 = sigma_2_0
        for t in 1:H

            r = rand(Normal(0, sqrt(sigma_2)))
            sigma_2 = omega + alpha * r^2 + beta * sigma_2
        end
        R_H[i] = r
    end
    return R_H
end