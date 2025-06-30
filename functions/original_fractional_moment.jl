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

function original_fractional_moment(dist::Distribution, α::Real; use_abs=false, complex=false)
    if dist isa ContinuousUnivariateDistribution
        μ = mean(dist)
        σ = std(dist)

        if use_abs
            lower = μ - 10σ
            upper = μ + 10σ
        else
            if is_symmetric(dist)
                lower = 0
            else
                lower = μ - 10σ
            end
            upper = μ + 10σ
        end

        integrand(x) = begin
            val = use_abs ? abs(x) : (complex ? (x + 0im) : x)
            val^α * pdf(dist, x)
        end

        moment, _ = quadgk(integrand, lower, upper)
        return (is_symmetric(dist) && lower ≥ 0 && !use_abs) ? 2 * moment : moment

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

function original_fractional_moment_CF(dist::Distribution, α::Real; use_abs=false, complex=false)
    n = floor(α)
    β = α - n

    if dist isa ContinuousUnivariateDistribution
         μ = mean(dist)
         σ = std(dist)

         if use_abs
            lower = μ - 10σ
            upper = μ + 10σ
         else
             if is_symmetric(dist)
                lower = 0
             else
                lower = μ - 10σ
             end
             upper = μ + 10σ
         end

        integrand(x) = begin
            ε = 1e-10
            val = use_abs ? abs(x) : (complex ? (x + 0im) : x)
            numer = val^(n + 1)
            denom = (1 - β) * val + β + ε
            numer / denom * pdf(dist, x)
        end

        moment, _ = quadgk(integrand, lower, upper)
        return (is_symmetric(dist) && lower ≥ 0 && !use_abs) ? 2 * moment : moment
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