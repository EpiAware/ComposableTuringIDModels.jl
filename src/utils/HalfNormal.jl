# Half-normal prior distribution parameterised by its mean.

@doc raw"
A half-normal prior distribution parameterised by its mean `μ`.

```math
X \sim \text{HalfNormal}(\mu) \iff X = |Y|,\quad
Y \sim \mathrm{Normal}\!\left(0, \mu\sqrt{\pi/2}\right)
```

so that ``\mathbb{E}[X] = \mu``.

# Examples
```jldoctest HalfNormal; output = false
using EpiAwarePrototype, Distributions
hn = HalfNormal(1.0)
nothing
# output
```
"
struct HalfNormal{T <: Real} <: ContinuousUnivariateDistribution
    μ::T
end

function Base.rand(rng::AbstractRNG, d::HalfNormal{T}) where {T <: Real}
    abs(rand(rng, Normal(0, d.μ * sqrt(π / 2))))
end

function Distributions.logpdf(d::HalfNormal{T}, x::Real) where {T <: Real}
    x < 0 ? -Inf : logpdf(Normal(0, d.μ * sqrt(π / 2)), x) - log(2)
end

function Distributions.cdf(d::HalfNormal{T}, x::Real) where {T <: Real}
    x < 0 ? 0.0 :
    cdf(Normal(0, d.μ * sqrt(π / 2)), x) -
    cdf(Normal(0, d.μ * sqrt(π / 2)), -x)
end

function Distributions.quantile(d::HalfNormal{T}, q::Real) where {T <: Real}
    quantile(Normal(0, d.μ * sqrt(π / 2)), q + (1 - q) / 2)
end

Base.minimum(d::HalfNormal) = 0.0
Base.maximum(d::HalfNormal) = Inf
Distributions.insupport(d::HalfNormal, x::Real) = x >= 0
Statistics.mean(d::HalfNormal{T}) where {T <: Real} = d.μ
Statistics.var(d::HalfNormal{T}) where {T <: Real} = d.μ^2 * (π / 2 - 1)
