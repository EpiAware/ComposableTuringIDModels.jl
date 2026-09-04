# Stochastic infections: the noise specification both parameterisations share,
# and the non-centred renewal modifier that applies it.

@doc raw"
Stochastic infections for a renewal process, giving the noise family and its width.

A deterministic renewal fixes infections at the renewal expectation ``\iota_t``.
`InfectionNoise` gives them a distribution of their own, matched to the first two moments of a negative binomial.

```math
\mathbb{E}[I_t] = \iota_t, \qquad
\mathrm{Var}[I_t] = \iota_t \left(1 + \iota_t \xi^2\right),
```

``\iota_t = R_t \sum_i I_{t-i} g_i`` is the renewal expectation and ``\xi`` the overdispersion, so ``\xi = 0`` leaves Poisson variance.
The coefficient of variation that follows is ``c_t = \sqrt{1/\iota_t + \xi^2}``, which diverges as ``\iota_t`` approaches zero.
An arbitrarily small expectation would then get an arbitrarily wide draw.
A smooth upper limit ``u - \mathrm{softplus}(u - c_t, k)`` holds it, with `cv_cap` = ``u`` and `cv_sharpness` = ``k``.
An infinite `cv_cap` removes the limit and restores the exact negative-binomial variance.

Being smooth costs the limit an absolute offset, ``\log(1 + e^{-k(u - c_t)}) / k``, which is 0.2% of the coefficient of variation at the defaults and ``6.7 \times 10^{-4}`` as ``c_t`` approaches zero.
Subtracting it from a smaller coefficient of variation would leave zero or less, which is no distribution at all, so the limit applies only while it returns a positive value.
Below that the limit is nowhere near binding, 0.5 against ``6.7 \times 10^{-4}``, and the exact coefficient of variation is used.
This is reached with `overdispersion = 0` past a few million infections.

`dist` is the family those moments are matched onto.
The default `LogNormal` keeps infections positive by construction, so no truncation is needed and no proposal can drive a downstream log-normal observation model to a negative mean.
`Normal` gives the linear form, at the cost of an unbounded support the renewal recursion can carry below zero.
Any family `ReparameterisedDistributions` registers for `(:mean, :sd)` works.

The same specification drives both parameterisations.

  - [`StochasticRenewal`](@ref) is **centred**.
    Infections are the sampled parameter and the likelihood informs them directly.
    Prefer it.
  - Used as a positional modifier on a [`Renewal`](@ref) it is **non-centred**.
    Standard normals are sampled and the location and scale applied inside the scan.
    That is what the modifier seam can express, and it is the worse parameterisation when the data are informative.

# Keyword Arguments

  - `dist`: the noise family (default `LogNormal`).
  - `overdispersion`: the prior for ``\xi``, or a fixed scalar that costs no parameter.
  - `cv_cap`: the soft upper limit on the coefficient of variation.
  - `cv_sharpness`: how sharply that limit is approached.

# Examples

A wider family and no cap, used by a [`StochasticRenewal`](@ref).

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
    function InfectionNoise(
            dist, overdispersion, cv_cap::C, cv_sharpness::C
        ) where {C}
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

# The negative-binomial standard deviation `sqrt(iota (1 + iota ξ²))`, written
# as `c_t ι` so the cap applies to the coefficient of variation.
#
# A non-positive expectation has no such moment pair, so it returns `NaN`, which
# the moment core's guards turn into a `-Inf` score or a `NaN` draw.
#
# The smooth cap costs an offset of `log1p(exp(-k(u - c)))/k`, 6.7e-4 at the
# defaults, and a small `ξ` past a few million infections makes the raw
# coefficient of variation smaller than that. The capped value is therefore used
# only while it is positive, and below that a cap of 0.5 is nowhere near
# binding.
function _noise_sd(noise::InfectionNoise, ι, ξ)
    ι > 0 || return convert(typeof(float(ι * ξ * noise.cv_cap)), NaN)
    raw = sqrt(inv(ι) + ξ^2)
    capped = _soft_upper(raw, noise.cv_cap, noise.cv_sharpness)
    return (capped > 0 ? capped : raw) * ι
end

# The per-step distribution the CENTRED form scores `I_t` against.
function _noise_dist(noise::InfectionNoise, ι, ξ)
    return _moment_dist(noise.dist, ι, _noise_sd(noise, ι, ξ))
end

# The per-step draw the NON-CENTRED form transforms a standard normal into.
function _noise_draw(noise::InfectionNoise, ι, ξ, z)
    return _moment_draw(noise.dist, ι, _noise_sd(noise, ι, ξ), z)
end

# The overdispersion slot, drawn once per model. A fixed scalar is a constant
# rather than a point-mass parameter, so it stays out of the sampled space.
# Kept as a `@model` so both parameterisations resolve the slot through the
# same call with no branch on what the slot holds. `n` is the series length, so
# a process-valued prior gives a time-varying overdispersion read per step with
# `at`.
@model function _noise_overdispersion(noise::InfectionNoise{<:Any, <:Real}, n)
    return noise.overdispersion
end

@model function _noise_overdispersion(noise::InfectionNoise, n)
    ξ ~ as_turing_submodel(noise.overdispersion, n; prefix = true)
    return ξ
end

# --- the non-centred modifier -----------------------------------------------

@doc raw"
A resolved [`InfectionNoise`](@ref) modifier: the drawn standard normals,
ready to scan.

This is what an [`InfectionNoise`](@ref) used as a renewal modifier returns
from its pre-scan `as_turing_model` seam, following the [`ImportedCases`](@ref)
→ [`ImportedRate`](@ref) pattern. Its substate is the step counter, so step
``t`` reads its own ``\tilde I_t``.

## Fields

  - `raw`: the standard normal draws ``\tilde I_t``, one per time, or a
    `strata x time` block for a stratified renewal.
  - `noise`: the [`InfectionNoise`](@ref) specification.
  - `overdispersion`: the resolved ``\xi``.
"
struct InfectionNoiseDraws{V, N, X} <: AbstractRenewalModifier
    "The standard normal draws, one per time, or per stratum per time."
    raw::V
    "The noise specification: the family and the soft cap."
    noise::N
    "The resolved overdispersion."
    overdispersion::X
end

# The substate is the step counter. A scan step has no clock of its own, so a
# per-time modifier carries the index it has reached. The window is unused,
# because a step counter has no shape of its own to match.
modifier_init_state(::InfectionNoiseDraws, window) = 0

function apply_modifier(mod::InfectionNoiseDraws, incidence, t)
    drawn = _noise_step(
        mod.noise, incidence, at(mod.overdispersion, t + 1), at(mod.raw, t + 1)
    )
    return drawn, t + 1
end

# `at` reads a per-time draw as a scalar and a `strata x time` draw as that
# step's column, so the two shapes differ only in whether the moment solve is
# broadcast. The scalar method keeps a single series free of broadcast
# machinery.
_noise_step(noise, ι::Real, ξ, z) = _noise_draw(noise, ι, ξ, z)

function _noise_step(noise, ι::AbstractVector, ξ, z)
    return _noise_draw.(Ref(noise), ι, ξ, z)
end

# The incidence entering this modifier is the last point at which it is still
# the renewal expectation.
is_noise(::InfectionNoise) = true
is_noise(::InfectionNoiseDraws) = true

# A stratified shape wraps the `IID` in a `Replicate`, which is the answer the
# shape guard gives any path model asked for a strata axis, so every stratum
# draws its own block instead of sharing one.
_noise_raw(::Int) = IID(Normal())
_noise_raw(::Dims{2}) = Replicate(IID(Normal()))

@doc raw"
Sample the **non-centred** infection noise ahead of the scan.

Draws `n` standard normals through the [`IID`](@ref) seam and resolves the overdispersion slot, returning the [`InfectionNoiseDraws`](@ref) the scan uses.
A fixed scalar `overdispersion` costs no parameter.
The draws are named `I_raw`, prefixed by the modifier's position in the renewal step's modifier tuple.

A `Dims{2}` shape wraps the [`IID`](@ref) seam in a [`Replicate`](@ref), giving a `strata x time` block of standard normals, one per stratum per step, which the scan reads a column at a time with [`at`](@ref).
Each stratum therefore gets its own noise about its own expectation, in the same way a [`Stratify`](@ref) or [`Replicate`](@ref) rate gives [`ImportedCases`](@ref) one exogenous stream per stratum.
The draws are then named `I_raw.stratum<g>.ϵ_t` under the modifier's position, following `Replicate`'s own naming.
A process-valued `overdispersion` still draws one path, so wrap it to vary that across strata too.

The parameterisation is non-centred because the modifier seam gives it no choice.
`apply_modifier` is deterministic and a modifier's priors resolve before the scan, so the sampled quantity has to be the standard normal ``\tilde I_t``, with the location and scale applied against ``\iota_t`` inside the scan.

!!! warning \"Non-centred is the known-worse parameterisation here\"
    When the observations inform the infections, which is the usual case for a renewal model, a non-centred latent path costs sampling efficiency.
    The cost shows up as maximum tree depth rather than as a wrong answer.
    [`StochasticRenewal`](@ref) draws the same noise centred and is what to reach for.
    Use this modifier when a non-centred draw is wanted, or for a stratified renewal, which the centred loop does not cover.

# Arguments

  - `mod`: the [`InfectionNoise`](@ref) specification.
  - `n`: the length of the renewal series.

# Examples
```@example InfectionNoiseModifier
using ComposableTuringIDModels, Distributions, Random
r = Renewal(
    [0.2, 0.3, 0.5], InfectionNoise();
    rt = FixedIntercept(0.0), initialisation = Normal(log(500.0), 0.0)
)
as_turing_model(r, 8)(Xoshiro(1)).I_t
```
"
@model function as_turing_model(mod::InfectionNoise, n)
    I_raw ~ as_turing_submodel(_noise_raw(n), n; prefix = true)
    ξ ~ to_submodel(_noise_overdispersion(mod, n), false)
    return InfectionNoiseDraws(I_raw, mod, ξ)
end
