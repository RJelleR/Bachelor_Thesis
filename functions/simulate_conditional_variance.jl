using Distributions

function simulate_conditional_variance(omega, alpha, beta, H, n_sim=100_000)
    R_H_2 = zeros(n_sim)
    sigma_2_0 = omega / (1 - alpha - beta)

    for i in 1:n_sim
        sigma_2 = sigma_2_0
        variance_R_H = 0.0 
        for _ in 1:H
            sigma = sqrt(sigma_2)
            r = rand(Normal(0, sigma))
            variance_R_H += r^2
            sigma_2 = omega + alpha * r^2 + beta * sigma_2
        end
        R_H_2[i] = variance_R_H
    end
    return R_H_2
end