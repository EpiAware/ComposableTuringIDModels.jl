# Moment-matched distributions: one guarded route from a `(mean, sd)` pair onto
# a distribution family, shared by the observation-error models
# (`ObservationError`) and by the stochastic renewal (`InfectionNoise`,
# `StochasticRenewal`).
#
# Three things are wanted from that route and none of them is free.
#
# 1. GUARDS. A sampler exploring an unconstrained space reaches arguments that
#    throw. `check_args = false` scores an invalid moment pair `-Inf` rather
#    than raising a `DomainError` mid-gradient, but it is not enough on its
#    own: the log-scale conversion `s² = log1p((sd / mean)²)` overflows to
#    `Inf` once `sd / mean` passes `sqrt(floatmax)`, and `Inf` then reaches the
#    constructor looking like a valid argument. `rand` validates and throws,
#    which is what a `missing` observation hits. So an invalid pair is routed
#    to a REJECTION distribution: `logpdf` is `-Inf` and `rand` returns `Inf`
#    without throwing.
# 2. CLOSED FORMS. `Normal` and `LogNormal` have a closed-form moment solve and
#    a closed-form quantile, and the generic path allocates a distribution per
#    call inside a differentiated loop. Writing the two out is measurably
#    cheaper and keeps the return type concrete.
# 3. TAILS. The generic draw round-trips through `cdf` then `quantile`, and
#    that round trip has no resolution left in the tails: `cdf` saturates at
#    exactly 1 by eight standard deviations out, so the draw comes back as the
#    support's endpoint, and three digits are already gone before it gets
#    there. Samplers do reach eight standard deviations.

# The field type a family is stored under. `typeof(LogNormal)` is `UnionAll`,
# which says nothing about *which* family it is, so every dispatch on the
# family would resolve at run time and whatever the moment core returns would
# infer as `Any` — which then propagates through everything downstream, a
# renewal recursion included. `Type{LogNormal}` has one instance, so the family
# stays in the type domain and the caller stays concretely typed.
#
# A component storing a family calls this from its inner constructor.
_family_type(::Type{F}) where {F} = Type{F}
_family_type(dist) = typeof(dist)

# The rejection distribution: `logpdf == -Inf` everywhere AND a `rand` that
# does not throw, which an invalid `Reparameterised` cannot give (its `rand`
# converts through `native`, which raises). A bare `Float64` sentinel is
# deliberate — giving it the input's type makes AD return `NaN` rather than
# `0.0`. `logcdf` is `-Inf` too, so it stays safe to censor or truncate.
_moment_reject(::Type{Normal}) = Normal(Inf, 1.0)
_moment_reject(::Type{Gamma}) = Gamma(1.0, Inf)
# Every other family falls back to the log-normal sentinel rather than one
# built from its own parameters: what is asked of a sentinel is a `-Inf` score
# and a finite-free draw, not membership of the family it stands in for.
_moment_reject(family) = LogNormal(Inf, 1.0)

# Whether a `(mean, sd)` pair describes a member of the family.
#
# The two closed-form families answer for themselves, in the coordinates their
# own conversion uses: `LogNormal` squares `sd / mean`, so it is that square
# that has to stay finite, and `Normal` takes the mean as its location and so
# accepts any finite one. Every other family defers to the registered
# `valid_moments` predicate, with the finiteness check `reparameterise` cannot
# make from inside a moment guard added on top.
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
# to `s² = log1p((sd / mean)²)` and `μ = log(mean) - s² / 2`, which is exactly
# what `reparameterise` registers for the family — the fast path is the same
# algebra without the wrapper, so the two agree to the last bit and the return
# type stays a plain `LogNormal`. `Normal` is native in these coordinates and
# has no registered reparameterisation at all.
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

The guarded entry point every moment-parameterised component builds its
per-time-point distribution through. An invalid pair scores `-Inf` and samples
`Inf`, so a diverging proposal is rejected rather than raising mid-gradient and
a `missing` observation can still be drawn.

# Arguments

  - `family`: the distribution family (a type, e.g. `LogNormal`).
  - `mean`: the mean the distribution is matched to.
  - `sd`: the standard deviation it is matched to.
"
function _moment_dist(family, mean, sd)
    _moment_valid(family, mean, sd) || return _moment_reject(family)
    return _native_dist(family, mean, sd)
end

# The moment-matched quantile at a standard normal draw `z`, on a pair already
# known valid.
#
# The two closed forms are exact: a `Normal` is its own location-scale family,
# and a `LogNormal` is one on the log scale. Everything else round-trips
# through the normal CDF, which is where the tails go.
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

`NaN` for a `(mean, sd)` pair describing no member of the family, which a
guarded observation model downstream turns into a `-Inf` score. The guard is
what keeps `quantile` — which raises on an invalid pair, as `rand` does — off
an unreachable branch.

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
