module FractionalMoments

using Distributions, StatsBase, QuadGK, Interpolations, KernelDensity, Statistics
import Statistics: mean, std
export KDEDistribution, fractional_moment, fractional_moment_CF

struct KDEDistribution <: Distributions.ContinuousUnivariateDistribution
    kde::KernelDensity.UnivariateKDE
end

Distributions.pdf(d::KDEDistribution, x::Real) = pdf(d.kde, x)
Statistics.mean(f::Function, d::KDEDistribution) = sum(f.(d.kde.x) .* d.kde.density) / sum(d.kde.density)
Statistics.mean(d::KDEDistribution) = mean(identity, d)

Statistics.std(d::KDEDistribution) = begin
    mu = mean(d)
    sqrt(sum(((d.kde.x .- mu).^2) .* d.kde.density) / sum(d.kde.density))
end

# Helper function to get good integration limits based on KDE quantiles and density cutoff
function integration_limits(dist::KDEDistribution, μ, σ; quantile_low=0.001, quantile_high=0.999, density_threshold=1e-5)
    xs = dist.kde.x
    densities = dist.kde.density

    # Quantile-based limits (cuts extreme outliers)
    lower_q = quantile(xs, 0.005)  # tighten to 0.005 from 0.001
    upper_q = quantile(xs, 0.995)
    density_threshold = 1e-4

    # Density threshold limits (ignore near-zero density tails)
    valid_idx = findall(densities .> density_threshold)
    lower_dens = minimum(xs[valid_idx])
    upper_dens = maximum(xs[valid_idx])

    # Combine limits with mean±std limits, picking the tighter bounds
    lower = maximum([lower_q, lower_dens, μ - 6σ, 0.0]) # also ensure ≥ 0 for fractional moments
    upper = minimum([upper_q, upper_dens, μ + 6σ])

    # Fallback if any limit is not finite
    if !isfinite(lower)
        lower = max(0.0, μ - 6σ)
    end
    if !isfinite(upper)
        upper = μ + 6σ
    end

    return lower, upper
end


function fractional_moment(dist::Distribution, α::Real; use_abs=false, complex=false, symmetric=false)
    f = x -> begin
        val = use_abs ? abs(x) : (complex ? (x + 0im) : x)
        val^α
    end

    if dist isa ContinuousUnivariateDistribution || dist isa KDEDistribution
        μ = mean(dist)
        σ = std(dist)

        if dist isa KDEDistribution
            lower, upper = integration_limits(dist, μ, σ)
        else
            # For other continuous distributions, keep previous heuristic
            lower = max(0.0, μ - 4σ)
            upper = μ + 6σ
            if !isfinite(lower) || lower < 0
                lower = symmetric ? 0.0 : -1000.0
            end
            if !isfinite(upper)
                upper = 1000.0
            end
        end

        integrand(x) = f(x) * pdf(dist, x)
        moment, _ = quadgk(integrand, lower, upper)
        return symmetric ? 2 * moment : moment

    elseif dist isa DiscreteUnivariateDistribution
        support_vals = 0:1000
        return sum(f(x) * pdf(dist, x) for x in support_vals)

    else
        throw(ArgumentError("Distribution type $(typeof(dist)) not supported"))
    end
end


function fractional_moment_CF(dist::Distribution, α::Real; use_abs=false, complex=false, symmetric=false)
    n = floor(Int, α)
    β = α - n

    f = x -> begin
        ε = 1e-10  # numerical stability
        val = use_abs ? abs(x) : (complex ? (x + 0im) : x)
        numer = val^(n + 1)
        denom = (1 - β) * val + β + ε
        numer / denom
    end

    if dist isa ContinuousUnivariateDistribution || dist isa KDEDistribution
        μ = mean(dist)
        σ = std(dist)

        if dist isa KDEDistribution
            lower, upper = integration_limits(dist, μ, σ)
        else
            lower = max(0.0, μ - 4σ)
            upper = μ + 6σ
            if !isfinite(lower) || lower < 0
                lower = symmetric ? 0.0 : -1000.0
            end
            if !isfinite(upper)
                upper = 1000.0
            end
        end

        integrand(x) = f(x) * pdf(dist, x)
        moment, _ = quadgk(integrand, lower, upper)
        return symmetric ? 2 * moment : moment

    elseif dist isa DiscreteUnivariateDistribution
        support_vals = 0:1000
        return sum(f(x) * pdf(dist, x) for x in support_vals)

    else
        throw(ArgumentError("Distribution type $(typeof(dist)) not supported"))
    end
end

end 
