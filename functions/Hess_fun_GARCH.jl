function Hess_fun_GARCH(par, x)
    n = length(x)
    omega = par[1]
    alpha = par[2]
    beta = par[3]
    sig_2 = zeros(eltype(par), n)
    sig_2[1] = var(x)

    for t in 2:n
        sig_2[t] = omega + alpha .* x[t - 1].^2 + beta .* sig_2[t - 1]
    end

    l =  .-(1/2) * log(2π) .- (1/2) * log.(sig_2) .- (1/2) * (x.^2) ./ sig_2
    llik = sum(l)

    return llik
 
end