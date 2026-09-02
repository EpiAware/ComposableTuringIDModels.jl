# Stochastic infections: the noise specification both parameterisations share.
#
# A deterministic renewal fixes infections at the renewal expectation. This
# gives them their own noise, matched to a negative binomial's first two
# moments. The specification below says which family and how wide; the
# CENTRED form ([`StochasticRenewal`](@ref)) draws `I_t` conditional on the
# expectation inside the loop, and the NON-CENTRED form (the renewal modifier
# in this file) draws standard normals before the scan and transforms them
# against the expectation inside it.

@doc raw"
Stochastic infections for a renewal process: the noise family and its width.

A deterministic renewal fixes infections at the renewal expectation
``\iota_t``. `InfectionNoise` gives them a distribution of their own, matched
to the first two moments of a negative binomial:

```math
\mathbb{E}[I_t] = \iota_t, \qquad
\mathrm{Var}[I_t] = \iota_t \left(1 + \iota_t \xi^2\right),
```

with ``\iota_t = R_t \sum_i I_{t-i} g_i`` the renewal expectation and ``\xi``
the overdispersion (``\xi = 0`` leaves Poisson variance). The coefficient of
variation that follows is ``c_t = \sqrt{1/\iota_t + \xi^2}``, which diverges as
``\iota_t`` approaches zero — an arbitrarily small expectation would otherwise
get an arbitrarily wide draw, which is a funnel. A smooth upper limit
``u - \mathrm{softplus}(u - c_t, k)`` holds it, with `cv_cap` = ``u`` and
`cv_sharpness` = ``k``. An infinite `cv_cap` removes the limit and restores
the exact negative-binomial variance.

Being smooth costs the limit an absolute offset,
``\log(1 + e^{-k(u - c_t)}) / k``, which is 0.2% of the coefficient of
variation at the defaults and a fixed ``6.7 \times 10^{-4}`` as ``c_t``
approaches zero. With `overdispersion = 0` and a large enough expectation that
offset takes the coefficient of variation to zero, leaving no valid moment
pair — which scores ``-\infty`` rather than raising. Pass `cv_cap = Inf`
there.

`dist` is the family those moments are matched onto, which the machinery
leaves free. The default `LogNormal` keeps infections positive by
construction, so no truncation is needed and no proposal can drive a
downstream log-normal observation model to a negative mean. `Normal` gives the
linear form, at the cost of an unbounded support the renewal recursion can
carry below zero. Any family `ReparameterisedDistributions` registers for
`(:mean, :sd)` works.

The same specification drives both parameterisations:

  - [`StochasticRenewal`](@ref) is **centred** — infections are the sampled
    parameter and the likelihood informs them directly. Prefer it.
  - Used as a positional modifier on a [`Renewal`](@ref) it is **non-centred**
    — standard normals are sampled and the location and scale applied inside
    the scan. That is what the modifier seam can express, and it is the worse
    parameterisation when the data are informative.

# Keyword Arguments

  - `dist`: the noise family (default `LogNormal`).
  - `overdispersion`: the prior for ``\xi``, or a fixed scalar. A fixed scalar
    costs no parameter.
  - `cv_cap`: the soft upper limit on the coefficient of variation.
  - `cv_sharpness`: how sharply that limit is approached.

# Examples

The specification on its own carries no priors beyond the overdispersion; it
is a [`StochasticRenewal`](@ref) (or a [`Renewal`](@ref) modifier) that uses
it. A wider family and no cap:

```@example InfectionNoise
using ComposableTuringIDModels, Distributions
noise = InfectionNoise(; dist = Normal, overdispersion = 0.3, cv_cap = Inf)
sr = StochasticRenewal(
    [0.2, 0.3, 0.5]; rt = FixedIntercept(0.0),
    initialisation = Normal(log(500.0), 0.0), noise = noise
)
using Random
as_turing_model(sr, 8)(Xoshiro(1)).I_t
```

## Fields

  - `dist`: the noise family, matched to the negative-binomial moments.
  - `overdispersion`: the prior for ``\xi``, or a fixed scalar.
  - `cv_cap`: the soft upper limit on the coefficient of variation.
  - `cv_sharpness`: the sharpness of that limit.
"
struct InfectionNoise{D, X, C} <: AbstractRenewalModifier
    "The noise family, matched to the negative-binomial moments."
    dist::D
    "Prior for the overdispersion ``\\xi``, or a fixed scalar."
    overdispersion::X
    "Soft upper limit on the coefficient of variation."
    cv_cap::C
    "Sharpness of the soft limit."
    cv_sharpness::C

    # Inner constructor so no default one is generated: the family has to
    # reach the field as `Type{F}` rather than as a `UnionAll` (see
    # `_family_type`), or every draw resolves its family at run time and the
    # incidence it returns infers as `Any` through the whole recursion.
    function InfectionNoise(dist, overdispersion, cv_cap::C, cv_sharpness::C) where {C}
        return new{_family_type(dist), typeof(overdispersion), C}(
            dist, overdispersion, cv_cap, cv_sharpness
        )
    end
end

function InfectionNoise(;
        dist = LogNormal, overdispersion = 0.1, cv_cap = 0.5,
        cv_sharpness = 10.0
    )
    cap, sharp = promote(cv_cap, cv_sharpness)
    return InfectionNoise(dist, overdispersion, cap, sharp)
end

# The standard deviation the negative-binomial moments imply at expectation
# `ι`, under the soft cap: `sqrt(iota (1 + iota ξ²))`, written as `c_t ι` so
# the cap applies to the coefficient of variation.
#
# A non-positive expectation has no such moment pair — `sqrt` of a negative
# number raises, and the renewal has already gone somewhere it cannot come
# back from — so it returns `NaN`, which the moment core's guards turn into a
# `-Inf` score (centred) or a `NaN` draw (non-centred) rather than an error.
function _noise_sd(noise::InfectionNoise, ι, ξ)
    ι > 0 || return convert(typeof(float(ι * ξ * noise.cv_cap)), NaN)
    cv = _soft_upper(sqrt(inv(ι) + ξ^2), noise.cv_cap, noise.cv_sharpness)
    return cv * ι
end

# The per-step distribution the CENTRED form scores `I_t` against.
function _noise_dist(noise::InfectionNoise, ι, ξ)
    return _moment_dist(noise.dist, ι, _noise_sd(noise, ι, ξ))
end

# The per-step draw the NON-CENTRED form transforms a standard normal into.
function _noise_draw(noise::InfectionNoise, ι, ξ, z)
    return _moment_draw(noise.dist, ι, _noise_sd(noise, ι, ξ), z)
end

# The overdispersion slot, drawn once per model. A fixed scalar is a constant,
# not a point-mass parameter: it stays out of the sampled space entirely,
# which is what a fixed overdispersion is for. Kept as a `@model` so both
# parameterisations resolve the slot through the same call with no branch on
# what the slot holds. `n` is the series length, so a process-valued prior
# gives a time-varying overdispersion read per step with `at`.
@model function _noise_overdispersion(noise::InfectionNoise{<:Any, <:Real}, n)
    return noise.overdispersion
end

@model function _noise_overdispersion(noise::InfectionNoise, n)
    ξ ~ as_turing_submodel(noise.overdispersion, n; prefix = true)
    return ξ
end
