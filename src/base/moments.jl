# Moment-matched distributions: one guarded route from a `(mean, sd)` pair onto
# a distribution family, shared by the observation-error models and by the
# stochastic renewal.

# The field type a family is stored under. `typeof(LogNormal)` is `UnionAll`,
# which says nothing about which family it is, so every dispatch on the family
# would resolve at run time and whatever the moment core returns would infer as
# `Any`. That propagates through everything downstream, a renewal recursion
# included. `Type{LogNormal}` has one instance, so the family stays in the type
# domain and the caller stays concretely typed.
#
# A component storing a family calls this from its inner constructor. That is
# also where an instance handed in place of the family is caught, rather than
# deep inside the moment solve.
_family_type(::Type{F}) where {F} = Type{F}

function _family_type(dist)
    return throw(
        ArgumentError(
            "a family is a distribution type, e.g. `LogNormal`, not an " *
                "instance of one ($(dist)). The component supplies the " *
                "moments the family is matched to."
        )
    )
end

# What an invalid moment pair routes to. Scoring and sampling want opposite
# things from it, and this is the one place that can tell them apart.
#
# Scoring is the sampler's hot path and a diverging proposal is routine, so
# `logpdf` is `-Inf`: silent, cheap, and the proposal is rejected. Sampling has
# no proposal to reject. A `rand` that returned a number would put an `Inf`
# into a simulated series and tell nobody, so it raises instead, naming the
# moments that have no distribution.
#
# The support matches the family it stands in for, so the bijector a linked
# `VarInfo` derives does not change with the validity of the moments.
struct _NoDraw{F, T <: Real} <: ContinuousUnivariateDistribution
    "The mean that has no distribution in this family."
    mean::T
    "The standard deviation that has no distribution in this family."
    sd::T
end

function _NoDraw{F}(mean, sd) where {F}
    m, s = promote(float(mean), float(sd))
    return _NoDraw{F, typeof(m)}(m, s)
end

Base.minimum(::_NoDraw{Normal}) = -Inf
Base.minimum(::_NoDraw) = 0.0
Base.maximum(::_NoDraw) = Inf
Distributions.insupport(::_NoDraw, ::Real) = true

# A bare `Float64` `-Inf` is deliberate: giving it the input's type makes AD
# return `NaN` where it should return the gradient it already had.
Distributions.logpdf(::_NoDraw, ::Real) = -Inf
Distributions.pdf(::_NoDraw, ::Real) = 0.0
Distributions.logcdf(::_NoDraw, ::Real) = -Inf
Distributions.cdf(::_NoDraw, ::Real) = 0.0
Distributions.logccdf(::_NoDraw, ::Real) = 0.0
Distributions.ccdf(::_NoDraw, ::Real) = 1.0

function _no_draw(d::_NoDraw{F}) where {F}
    return throw(
        ArgumentError(
            "no $(nameof(F)) has mean $(d.mean) and standard deviation " *
                "$(d.sd), so this model cannot produce a value here. " *
                "Scoring such a point rejects it, but simulating one has " *
                "nothing to fall back on."
        )
    )
end

Base.rand(::AbstractRNG, d::_NoDraw) = _no_draw(d)
Distributions.quantile(d::_NoDraw, ::Real) = _no_draw(d)

_moment_reject(::Type{F}, mean, sd) where {F} = _NoDraw{F}(mean, sd)

# Whether a `(mean, sd)` pair describes a member of the family.
#
# The two closed-form families answer for themselves, in the coordinates their
# own conversion uses. `LogNormal` squares `sd / mean`, which overflows to `Inf`
# once the ratio passes `sqrt(floatmax)` and then reaches the constructor
# looking valid. `Normal` takes the mean as its location, so it accepts any
# finite one. Every other family defers to the registered `valid_moments`
# predicate, with a finiteness check on top.
function _moment_valid(::Type{Normal}, mean, sd)
    return isfinite(mean) && isfinite(sd) && sd >= 0
end

function _moment_valid(::Type{LogNormal}, mean, sd)
    return isfinite(mean) && mean > 0 && isfinite(sd) && sd > 0 &&
        isfinite((sd / mean)^2)
end

function _moment_valid(family, mean, sd)
    isfinite(mean) && isfinite(sd) || return false
    vals = promote(float(mean), float(sd))
    return valid_moments(family, Val((:mean, :sd)), vals)
end

# The moment solve itself, on a pair already known valid. `LogNormal` inverts
#   mean = exp(μ + s² / 2),  var = mean² (exp(s²) - 1)
# to `s² = log1p((sd / mean)²)` and `μ = log(mean) - s² / 2`, which is what
# `reparameterise` registers for the family. The fast path is the same algebra
# without the wrapper, so the two agree to the last bit and the return type
# stays a plain `LogNormal`. `Normal` is native in these coordinates and has no
# registered reparameterisation.
_native_dist(::Type{Normal}, mean, sd) = Normal(mean, sd)

function _native_dist(::Type{LogNormal}, mean, sd)
    s² = log1p((sd / mean)^2)
    return LogNormal(log(mean) - s² / 2, sqrt(s²))
end

function _native_dist(family, mean, sd)
    return reparameterise(family; mean = mean, sd = sd, check_args = false)
end

@doc raw"
The distribution of family `family` with mean `mean` and standard deviation
`sd`, or a rejection distribution when that pair describes no member of the
family.

The guarded entry point every moment-parameterised component builds its per-time-point distribution through.
An invalid pair scores ``-\infty``, so a diverging proposal is rejected rather than raising mid-gradient.
Asking one for a value raises instead, because there is no proposal to reject and no honest number to return.

# Arguments

  - `family`: the distribution family (a type, e.g. `LogNormal`).
  - `mean`: the mean the distribution is matched to.
  - `sd`: the standard deviation it is matched to.
"
function _moment_dist(family, mean, sd)
    _moment_valid(family, mean, sd) || return _moment_reject(family, mean, sd)
    return _native_dist(family, mean, sd)
end

# The moment-matched quantile at a standard normal draw `z`, on a pair already
# known valid.
#
# The two closed forms are exact. A `Normal` is its own location-scale family
# and a `LogNormal` is one on the log scale. Everything else round-trips through
# the normal CDF, which has no resolution left in the tails. It saturates at
# exactly 1 by eight standard deviations out, and the draw comes back as the
# support's endpoint.
_moment_quantile(::Type{Normal}, mean, sd, z) = mean + sd * z

function _moment_quantile(::Type{LogNormal}, mean, sd, z)
    s² = log1p((sd / mean)^2)
    return exp(log(mean) - s² / 2 + z * sqrt(s²))
end

function _moment_quantile(family, mean, sd, z)
    return quantile(_native_dist(family, mean, sd), cdf(Normal(), z))
end

@doc raw"
A non-centred draw from the family `family` matched to `mean` and `sd`: the
quantile at the standard normal draw `z`.

```math
X = F^{-1}_{m, s}(\Phi(z)), \qquad z \sim \mathrm{Normal}(0, 1)
```

`NaN` for a `(mean, sd)` pair describing no member of the family, which a guarded observation model downstream turns into a ``-\infty`` score.
The guard keeps `quantile` off such a pair, which it raises on as `rand` does.

# Arguments

  - `family`: the distribution family (a type, e.g. `LogNormal`).
  - `mean`: the mean the draw is matched to.
  - `sd`: the standard deviation it is matched to.
  - `z`: the standard normal draw.
"
function _moment_draw(family, mean, sd, z)
    _moment_valid(family, mean, sd) ||
        return convert(typeof(float(mean + sd * z)), NaN)
    return _moment_quantile(family, mean, sd, z)
end

# --- the soft cap -----------------------------------------------------------

# `softplus(x, k) = log1p(exp(k x)) / k`, through `log1pexp` because the naive
# form overflows to `Inf` once `k x` passes about 709 and a diverging proposal
# reaches that. `log1pexp` returns its argument there instead.
_softplus(x, k) = log1pexp(k * x) / k

@doc raw"
A smooth upper limit on `x`,

```math
\mathrm{soft\_upper}(x, u, k) = u - \mathrm{softplus}(u - x, k),
```

which tracks `x` well below the limit `u` and saturates at it above. `k` sets
how sharply the two meet. An infinite `u` is no limit at all and returns `x`
unchanged, which is how a component offers the exact, uncapped quantity.

# Arguments

  - `x`: the quantity being limited.
  - `u`: the upper limit.
  - `k`: the sharpness of the limit.
"
_soft_upper(x, u, k) = isfinite(u) ? u - _softplus(u - x, k) : x
