# Negative-binomial distribution that avoids `InexactError` at very large
# means, plus the mean/cluster-factor constructor.

@doc raw"
A negative binomial distribution parameterised by `(r, p)` that avoids
`InexactError` at very large means.

The package uses a mean/cluster-factor parameterisation when constructing this
distribution from an expected count (see
[`NegativeBinomialMeanClust`](@ref)).

# Examples
```jldoctest SafeNegativeBinomial; output = false
using ComposableTuringIDModels, Distributions
bigμ = exp(48.0)
σ² = bigμ + 0.05 * bigμ^2
p = bigμ / σ²
r = bigμ * p / (1 - p)
d = SafeNegativeBinomial(r, p)
logpdf(d, 100)
nothing
# output
```
"
struct SafeNegativeBinomial{T <: Real} <: SafeDiscreteUnivariateDistribution
    r::T
    p::T

    SafeNegativeBinomial{T}(r::T, p::T) where {T <: Real} = new{T}(r, p)
end

SafeNegativeBinomial(r::T, p::T) where {T <: Real} = SafeNegativeBinomial{T}(r, p)
SafeNegativeBinomial(r::Real, p::Real) = SafeNegativeBinomial(promote(r, p)...)

_negbin(d::SafeNegativeBinomial) = NegativeBinomial(d.r, d.p; check_args = false)

Base.minimum(d::SafeNegativeBinomial) = 0
Base.maximum(d::SafeNegativeBinomial) = Inf
Distributions.insupport(d::SafeNegativeBinomial, x::Integer) = x >= 0
Distributions.params(d::SafeNegativeBinomial) = _negbin(d) |> params
Distributions.partype(::SafeNegativeBinomial{T}) where {T} = T
Distributions.succprob(d::SafeNegativeBinomial) = _negbin(d).p
Distributions.failprob(d::SafeNegativeBinomial{T}) where {T} = one(T) - _negbin(d).p
Distributions.mean(d::SafeNegativeBinomial) = _negbin(d) |> mean
Distributions.var(d::SafeNegativeBinomial) = _negbin(d) |> var
Distributions.std(d::SafeNegativeBinomial) = _negbin(d) |> std
Distributions.skewness(d::SafeNegativeBinomial) = _negbin(d) |> skewness
Distributions.kurtosis(d::SafeNegativeBinomial) = _negbin(d) |> kurtosis
Distributions.mode(d::SafeNegativeBinomial) = _negbin(d) |> mode
Distributions.logpdf(d::SafeNegativeBinomial, k::Real) = logpdf(_negbin(d), k)
Distributions.cdf(d::SafeNegativeBinomial, x::Real) = cdf(_negbin(d), x)
Distributions.ccdf(d::SafeNegativeBinomial, x::Real) = ccdf(_negbin(d), x)
Distributions.logcdf(d::SafeNegativeBinomial, x::Real) = logcdf(_negbin(d), x)
Distributions.logccdf(d::SafeNegativeBinomial, x::Real) = logccdf(_negbin(d), x)
Distributions.quantile(d::SafeNegativeBinomial, q::Real) = quantile(_negbin(d), q)

function Base.rand(rng::AbstractRNG, d::SafeNegativeBinomial)
    if isone(d.p)
        return 0
    else
        return rand(rng,
            SafePoisson(rand(rng, Gamma(d.r, (1 - d.p) / d.p; check_args = false))))
    end
end

@doc raw"
Construct a [`SafeNegativeBinomial`](@ref) from a mean `μ` and cluster factor
`α` using the variance relationship ``\sigma^2 = \mu + \alpha\,\mu^2``.

# Arguments

  - `μ`: the mean of the distribution.
  - `α`: the cluster factor relating mean and variance.

# Examples
```@example NegativeBinomialMeanClust
using ComposableTuringIDModels
NegativeBinomialMeanClust(10.0, 0.1)
```
"
function NegativeBinomialMeanClust(μ, α)
    μ² = μ^2
    σ² = μ + α * μ²
    p = μ / σ²
    r = 1 / α
    return SafeNegativeBinomial(r, p)
end
