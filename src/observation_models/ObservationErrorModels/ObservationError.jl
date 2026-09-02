# The generic moment-parameterised observation-error model, and the two
# location-scale families that are common enough to keep a name of their own.
#
# One type per family means adding a family means adding a type, and the types
# differ only in which family they name and in whether the spread is absolute
# or relative to the expected value. `ObservationError` takes both as
# arguments. `NormalError` and `LogNormalError` stay as named components,
# because a name is worth having for the two families most models reach for,
# and both are the same two lines over the shared moment core.

@doc raw"
An observation-error model parameterised by the moments of any family.

Each observed value has the expected value as its mean and a sampled spread as
its standard deviation, with the family's native parameters solved for from
that pair:

```math
\mathbb{E}[y_t] = Y_t, \qquad \mathrm{sd}[y_t] = \sigma \ \text{or} \ \sigma Y_t,
\qquad y_t \sim F(Y_t, \sigma)
```

`relative` chooses between the two: `false` (the default) makes ``\sigma`` an
absolute standard deviation, `true` makes it a coefficient of variation, so
the noise scales with the expected value. A relative log-normal is the usual
choice for a concentration or another strictly positive continuous
measurement; an absolute normal is [`NormalError`](@ref).

The family is any `ReparameterisedDistributions` registers for `(:mean, :sd)`,
which is how a new family becomes a call rather than a new type.
`Normal` and `LogNormal` take a closed-form path instead of the moment solve.
An invalid moment pair — which a diverging sampler reaches — scores `-\infty`
rather than raising mid-gradient, and can still be drawn when the observation
is `missing`.

`sd` sets the prior for ``\sigma`` — a `Distribution` (a constant, one scalar
RV) or a process (a length-`n`, e.g. time-varying, spread). It is drawn
through the single [`as_turing_submodel`](@ref) seam and read per time point
via `at`, so a process makes the observation noise time-varying with no other
change. It is sampled under the name `σ`.

# Arguments

  - `dist`: the observation family, e.g. `LogNormal`.

# Keyword Arguments

  - `sd`: the prior for the spread ``\sigma``.
  - `relative`: whether the spread scales with the expected value.

# Examples
```@example ObservationError
using ComposableTuringIDModels, Distributions
# A relative-noise gamma: 20% coefficient of variation about the expectation.
oe = ObservationError(Gamma; sd = HalfNormal(0.2), relative = true)
rand(as_turing_model(oe, missing, fill(100.0, 5)))
```

## Fields

  - `dist`: the observation family.
  - `sd`: the prior for the spread.
  - `relative`: whether the spread scales with the expected value.
"
struct ObservationError{D, S <: PriorLike} <: AbstractObservationErrorModel
    "The observation family, matched to the expected value and the spread."
    dist::D
    "Prior for the spread."
    sd::S
    "Whether the spread scales with the expected value."
    relative::Bool

    # Inner constructor so no default one is generated: the family has to
    # reach the field as `Type{F}` rather than as a `UnionAll` (see
    # `_family_type`).
    function ObservationError(dist, sd::PriorLike, relative::Bool)
        return new{_family_type(dist), typeof(sd)}(dist, sd, relative)
    end
end

function ObservationError(dist; sd = HalfNormal(0.1), relative::Bool = false)
    return ObservationError(dist, sd, relative)
end

@model function generate_observation_error_priors(
        obs_model::ObservationError, y_t, Y_t
    )
    σ ~ as_turing_submodel(obs_model.sd, length(Y_t); prefix = true)
    return (; σ = σ)
end

function observation_error(obs_model::ObservationError, Y_t, σ)
    return _moment_dist(
        obs_model.dist, Y_t, obs_model.relative ? σ * Y_t : σ
    )
end

@doc raw"
A log-normal observation-error model with an inferred coefficient of variation.

The noise is **relative**: the distribution's real-space mean is the expected
value and its real-space standard deviation is proportional to it, so
``\sigma`` is a coefficient of variation rather than an absolute scale.

```math
\mathbb{E}[y_t] = Y_t, \qquad \mathrm{sd}[y_t] = \sigma Y_t, \qquad
y_t \sim \mathrm{LogNormal}(\mu_t, s)
```

A `LogNormal` is parameterised on the log scale, so those two moments are
converted to its native parameters:

```math
s^2 = \log(1 + \sigma^2), \qquad \mu_t = \log Y_t - s^2 / 2
```

``s`` does not depend on ``Y_t``: a constant coefficient of variation is a
constant variance on the log scale, which is what makes this the
relative-noise family. It is the natural error model for a strictly positive
continuous measurement — a concentration, a prevalence — where the spread
grows with the level, which [`NormalError`](@ref) cannot express.

`cv` sets the prior for ``\sigma``, drawn through the single
[`as_turing_submodel`](@ref) seam and read per time point via `at`, so a
process makes the noise time-varying. It is sampled under the name `cv`.

Equivalent to `ObservationError(LogNormal; sd = cv, relative = true)`, which
is what to reach for with any other family.

# Examples
```@example LogNormalError
using ComposableTuringIDModels, Distributions
lne = LogNormalError()
rand(as_turing_model(lne, missing, fill(100.0, 5)))
```

## Fields

  - `cv`: prior for the coefficient of variation.
"
struct LogNormalError{S <: PriorLike} <: AbstractObservationErrorModel
    "Prior for the coefficient of variation."
    cv::S
end

# A 10% coefficient of variation on average, matching `NormalError`'s own
# default scale.
LogNormalError(; cv = HalfNormal(0.1)) = LogNormalError(cv)

@model function generate_observation_error_priors(
        obs_model::LogNormalError, y_t, Y_t
    )
    cv ~ as_turing_submodel(obs_model.cv, length(Y_t); prefix = true)
    return (; cv = cv)
end

observation_error(::LogNormalError, Y_t, σ) = _moment_dist(LogNormal, Y_t, σ * Y_t)
