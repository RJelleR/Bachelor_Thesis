module FractionalMoments
using Distributions, Statistics, QuadGK

function is_well_defined(dist::Distribution)
    lower = minimum(dist)
    return isfinite(lower) && lower ≥ 0
end

function is_symmetric(dist::Distribution)
    μ = mean(dist)
    return pdf(dist, μ + 1.0) == pdf(dist, μ - 1.0)
end

function fractional_moment(dist::Distribution, α::Real; use_abs=false, complex=false)
    if dist isa ContinuousUnivariateDistribution
        lower = minimum(dist)
        upper = maximum(dist)
        print(lower)
        if !isfinite(lower)
            lower = -1000
            print(lower)
            if lower < 0 && is_symmetric(dist)
                lower = 0
            end
        elseif lower < 0 && is_symmetric(dist)
            print("ok")
            lower = 0
        elseif lower < 0
            lower = -1000
        end

        if !isfinite(upper)
            upper = 1000
        end

        integrand(x) = begin
            val = use_abs ? abs(x) : (complex ? (x + 0im) : x)
            val^α * pdf(dist, x)
        end

        moment, _ = quadgk(integrand, lower, upper)
        return is_symmetric(dist) ? 2 * moment : moment

    elseif dist isa DiscreteUnivariateDistribution
        support_vals = 0:1000
        return sum(begin
            val = use_abs ? abs(x) : (complex ? (x + 0im) : x)
            val^α * pdf(dist, x)
        end for x in support_vals)
    else
        throw(ArgumentError("Distribution type not supported"))
    end
end

function fractional_moment_CF(dist::Distribution, α::Real; use_abs=false, complex=false)
    n = floor(α)
    β = α - n

    if dist isa ContinuousUnivariateDistribution
        lower = minimum(dist)
        upper = maximum(dist)

        if !isfinite(lower)
            lower = -1000
            if lower < 0 && is_symmetric(dist)
                lower = 0
            end
        elseif lower < 0 && is_symmetric(dist)
            lower = 0
        elseif lower < 0
            lower = -1000
        end

        if !isfinite(upper)
            upper = 1000
        end

        integrand(x) = begin
            ε = 1e-10
            val = use_abs ? abs(x) : (complex ? (x + 0im) : x)
            numer = val^(n + 1)
            denom = (1 - β) * val + β + ε
            numer / denom * pdf(dist, x)
        end

        moment, _ = quadgk(integrand, lower, upper)
        print(moment)
        return is_symmetric(dist) ? 2 * moment : moment
    elseif dist isa DiscreteUnivariateDistribution
        support_vals = 0:1000
        return sum(begin
            ε = 1e-10
            val = use_abs ? abs(x) : (complex ? (x + 0im) : x)
            numer = val^(n + 1)
            denom = (1 - β) * val + β + ε
            numer / denom * pdf(dist, x)
        end for x in support_vals)
    else
        throw(ArgumentError("Distribution type not supported"))
    end
end
end