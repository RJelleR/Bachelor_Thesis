function llik_fun_GARCH_pq(par, x, p, q)
    n = length(x)
    ω = exp(par[1])
    α = exp.(par[2: 1 + q]) ./ (1 .+ exp.(par[2:(1 + q)]))
    β = exp.(par[(2 + q):(1 + q + p)]) ./ (1 .+ exp.(par[(2 + q): (1 + q + p)]))
    sig_2 = fill(var(x), n)
    for t in (max(p, q) + 1):n
        arch_sum = sum(α .* x[(t -1):-1:(t - q)] .^2)
        garch_sum = sum(β .* sig_2[(t - 1):-1:(t - p)])
        sig_2[t] = ω + arch_sum + garch_sum
    end

    l = .-(1/2) * log(2π) .- (1/2) * log.(sig_2) .- (1/2) * (x.^2) ./ sig_2
    llik = mean(l)

    return llik
end 