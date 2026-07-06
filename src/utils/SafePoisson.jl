# Poisson distribution that avoids `InexactError` at very large means.

@doc raw"
A Poisson distribution parameterised by its mean `λ` that avoids `InexactError`
for very large means.

# Examples
```jldoctest SafePoisson; output = false
using ComposableTuringIDModels, Distributions
d = SafePoisson(exp(48.0))
logpdf(d, 100)
nothing
# output
```
"
struct SafePoisson{T <: Real} <: SafeDiscreteUnivariateDistribution
    λ::T

    SafePoisson{T}(λ::Real) where {T <: Real} = new{T}(λ)
    SafePoisson(λ::Real) = SafePoisson{eltype(λ)}(λ)
end

SafePoisson() = SafePoisson{Float64}(1.0)

_poisson(d::SafePoisson) = Poisson(d.λ; check_args = false)

Distributions.params(d::SafePoisson) = _poisson(d) |> params
Distributions.partype(::SafePoisson{T}) where {T} = T
Distributions.rate(d::SafePoisson) = d.λ
Distributions.mean(d::SafePoisson) = d.λ
Distributions.mode(d::SafePoisson) = floor(d.λ)
Distributions.var(d::SafePoisson) = d.λ
Distributions.skewness(d::SafePoisson) = one(typeof(d.λ)) / sqrt(d.λ)
Distributions.kurtosis(d::SafePoisson) = one(typeof(d.λ)) / d.λ
Distributions.logpdf(d::SafePoisson, x::Real) = logpdf(_poisson(d), x)
Distributions.pdf(d::SafePoisson, x::Integer) = pdf(_poisson(d), x)
Distributions.cdf(d::SafePoisson, x::Integer) = cdf(_poisson(d), x)
Distributions.ccdf(d::SafePoisson, x::Integer) = ccdf(_poisson(d), x)
Distributions.quantile(d::SafePoisson, q::Real) = quantile(_poisson(d), q)

Base.minimum(d::SafePoisson) = 0
Base.maximum(d::SafePoisson) = Inf
Distributions.insupport(d::SafePoisson, x::Integer) = x >= 0

# Poisson sampling that stays exact for arbitrarily large means. Small means use
# inversion by exponential interarrivals; large means use the Ahrens-Dieter
# normal-based method. Both floor to an integer via `_safe_int_floor`, which
# promotes to `BigInt` rather than throwing `InexactError` when the mean exceeds
# the `Int` range (which `Distributions.Poisson`'s own sampler does not handle).
#
# Ahrens-Dieter (1982), "Computer Generation of Poisson Deviates from Modified
# Normal Distributions", ACM TOMS 8(2):163-179. Ported from the original
# `EpiAware` package, which adapted it from PoissonRandom.jl.

function _count_rand(rng::AbstractRNG, λ)
    n = 0
    c = randexp(rng)
    while c < λ
        n += 1
        c += randexp(rng)
    end
    return n
end

# log(1+x) - x, accurate over the reduced range used by `_log1pmx`.
function _log1pmx_kernel(x::Float64)
    r = x / (x + 2.0)
    t = r * r
    w = @evalpoly(t, 6.66666666666666667e-1, 4.00000000000000000e-1,
        2.85714285714285714e-1, 2.22222222222222222e-1, 1.81818181818181818e-1,
        1.53846153846153846e-1, 1.33333333333333333e-1, 1.17647058823529412e-1)
    hxsq = 0.5 * x * x
    return r * (hxsq + w * t) - hxsq
end

function _log1pmx(x::Float64)
    if !(-0.7 < x < 0.9)
        return log1p(x) - x
    elseif x > 0.315
        u = (x - 0.5) / 1.5
        return _log1pmx_kernel(u) - 9.45348918918356180e-2 - 0.5 * u
    elseif x > -0.227
        return _log1pmx_kernel(x)
    elseif x > -0.4
        u = (x + 0.25) / 0.75
        return _log1pmx_kernel(u) - 3.76820724517809274e-2 + 0.25 * u
    elseif x > -0.6
        u = (x + 0.5) * 2.0
        return _log1pmx_kernel(u) - 1.93147180559945309e-1 + 0.5 * u
    else
        u = (x + 0.625) / 0.375
        return _log1pmx_kernel(u) - 3.55829253011726237e-1 + 0.625 * u
    end
end

# Procedure F of the Ahrens-Dieter algorithm.
function _procf(λ, K::SafeInt, s::Float64)
    ω = 0.3989422804014327 / s
    b1 = 0.041666666666666664 / λ
    b2 = 0.3 * b1 * b1
    c3 = 0.14285714285714285 * b1 * b2
    c2 = b2 - 15.0 * c3
    c1 = b1 - 6.0 * b2 + 45.0 * c3
    c0 = 1.0 - b1 + 3.0 * b2 - 15.0 * c3
    if K < 10
        px = -float(λ)
        py = λ^K / factorial(K)
    else
        δ = 0.08333333333333333 / K
        δ -= 4.8 * δ^3
        V = (λ - K) / K
        px = K * _log1pmx(V) - δ
        py = 0.3989422804014327 / sqrt(K)
    end
    X = (K - λ + 0.5) / s
    X2 = X^2
    fx = -0.5 * X2
    fy = ω * (((c3 * X2 + c2) * X2 + c1) * X2 + c0)
    return px, py, fx, fy
end

function _ad_rand(rng::AbstractRNG, λ)
    s = sqrt(λ)
    d = 6.0 * λ^2
    L = _safe_int_floor(λ - 1.1484)
    G = λ + s * randn(rng)
    if G >= 0.0
        K = _safe_int_floor(G)
        K >= L && return K
        U = rand(rng)
        d * U >= (λ - K)^3 && return K
        px, py, fx, fy = _procf(λ, K, s)
        fy * (1 - U) <= py * exp(px - fx) && return K
    end
    while true
        E = randexp(rng)
        U = 2.0 * rand(rng) - 1.0
        T = 1.8 + copysign(E, U)
        T <= -0.6744 && continue
        K = _safe_int_floor(λ + s * T)
        px, py, fx, fy = _procf(λ, K, s)
        c = 0.1069 / λ
        @fastmath if c * abs(U) <= py * exp(px + E) - fy * exp(fx + E)
            return K
        end
    end
end

function Base.rand(rng::AbstractRNG, d::SafePoisson)
    d.λ < 6 ? _count_rand(rng, d.λ) : _ad_rand(rng, d.λ)
end
