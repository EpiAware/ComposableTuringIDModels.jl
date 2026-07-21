# Right-truncation correction observation modifier (nowcasting): scale the
# expected eventual totals by a composable reporting-completeness (CDF) submodel.

@doc raw"
A composable reporting-completeness component: a `Distribution` (or a precomputed
vector) turned into the cumulative reporting proportion by age, for use as the
correction submodel of [`RightTruncate`](@ref).

Given a series length `n`, sampling `as_turing_model(c::ReportingCDF, n)` returns a
length-`n` vector `F` where `F[a + 1]` is the fraction of a reference day's
eventual total reported within `a` days (`a = 0, 1, …, n - 1`), in `[0, 1]`.
Reference days older than the reporting delay's support are fully reported
(`F = 1`), so the vector is padded with ones to length `n`. The CDF built from a
delay distribution is non-decreasing, but `ReportingCDF` does **not** require
monotonicity — a precomputed curve may be non-monotonic, so an over-/under-reporting
correction that recovers can be expressed.

It is the **fixed-delay** default used by [`RightTruncate`](@ref): the completeness
is precomputed once and held constant. Because the correction is supplied to
`RightTruncate` as a submodel, a user can instead pass any latent component that
produces a length-`n` completeness curve — a flexible non-parametric CDF, or even
a non-monotonic correction — without changing `RightTruncate`.

## Constructors

  - `ReportingCDF(distribution; D, Δd)` — discretise a continuous reporting-delay
    distribution via double-interval censoring (CensoredDistributions.jl) and take
    the cumulative sum of the resulting PMF, exactly the released-CD path
    [`LatentDelay`](@ref) uses.
  - `ReportingCDF(cdf)` — from a precomputed completeness vector by age (in
    `[0, 1]`; need not be monotonic).

# Examples
```@example ReportingCDF
using ComposableTuringIDModels, Distributions
c = ReportingCDF(truncated(Normal(5.0, 2.0), 0.0, Inf))
as_turing_model(c, 10)()
```

## Fields

  - `cdf`: the reporting-completeness curve by age (`cdf[a + 1]`), in `[0, 1]`
    (need not be monotonic). Padded with ones up to the requested length when
    shorter.
"
struct ReportingCDF{T <: AbstractVector{<:Real}} <: AbstractLatentModel
    "The reporting-completeness CDF by age."
    cdf::T

    function ReportingCDF(cdf::T) where {T <: AbstractVector{<:Real}}
        @assert all(>=(0), cdf) "The reporting completeness must be non-negative"
        @assert all(<=(1 + 1e-8), cdf) "The reporting completeness must not exceed 1"
        # No monotonicity check: a reporting-delay CDF is non-decreasing, but the
        # correction is deliberately a free completeness curve so a user can supply
        # a non-monotonic correction (e.g. over-/under-reporting that recovers).
        return new{T}(cdf)
    end
end

function ReportingCDF(distribution::C; D = nothing, Δd = 1.0) where {
        C <: ContinuousDistribution}
    # Build the reporting-delay CDF from the released-CD double-interval-censored
    # PMF (the same path `LatentDelay` uses), then accumulate it.
    pmf = _discretised_pmf(distribution; Δd = Δd, D = D)
    return ReportingCDF(cumsum(pmf))
end

@model function as_turing_model(c::ReportingCDF, n)
    F = c.cdf
    nF = length(F)
    # Reference days older than the delay's support (age ≥ `nF`) are fully
    # reported, so pad the completeness with ones up to age `n - 1`.
    return nF >= n ? F[1:n] : vcat(F, ones(eltype(F), n - nF))
end

@doc raw"
Correct an underlying observation model for **right-truncation** (not-yet-reported
counts), the EpiNow2-style CDF-scaling nowcast.

The infection → expected-observation pipeline produces ``Y_t = \mu_t``, the
expected *eventual* total for reference day `t`. At the present time `now` (taken
to be the last reference day, `now = n`) a reference day of age ``a = now - t`` has
only had a fraction ``F[a + 1]`` of its eventual total reported, where ``F`` is the
reporting-completeness CDF. The **expected observed-so-far** is therefore
``\mu_t \cdot F[a + 1]``, and conditioning the inner observation error on that —
rather than on the full ``\mu_t`` — corrects the right-truncation. Recent,
still-maturing reference days are automatically down-weighted, while the model's
`Y_t` remains the eventual total (so the nowcast of the eventual total is just
`Y_t` read out as a generated quantity).

The modifier mirrors [`Ascertainment`](@ref): it wraps an inner observation model,
draws a length-`n` correction series from a **submodel**, transforms the
expected-observation vector, and delegates via `to_submodel(..., false)`. The
correction `cdf_model` is any component producing a length-`n` completeness curve
`F` (by age). The default is a fixed [`ReportingCDF`](@ref) built from a reporting
delay distribution (the released-CD case), but because the correction is a submodel
a user can supply a flexible non-parametric CDF — or even a non-monotonic
correction — without changing `RightTruncate`.

The completeness is indexed by **age**: the most recent reference day (`t = n`,
age `0`) is scaled by `F[1]` and the oldest (`t = 1`, age `n - 1`) by `F[n]`, so
the age-indexed series is **reversed** onto the reference-day axis. A
fully-reported correction (all ones) leaves the inner model unchanged.

This is the **fixed-delay** variant. An estimated / time-varying delay (a latent
correction with sampled parameters) is a planned follow-up, and is exactly the
submodel slot generalising to it.

## Constructors

  - `RightTruncate(model, cdf_model)` — from an inner observation model and a
    correction component (a latent-role submodel producing the length-`n`
    completeness curve, e.g. a [`ReportingCDF`](@ref)).
  - `RightTruncate(model, distribution; D, Δd)` — wrap a continuous reporting-delay
    distribution in the default [`ReportingCDF`](@ref) (the released-CD path).
  - `RightTruncate(model, cdf::AbstractVector)` — wrap a precomputed completeness
    vector in a fixed [`ReportingCDF`](@ref).

# Arguments

  - `obs_model`: the [`RightTruncate`](@ref) model.
  - `y_t`: the observed-so-far series (or `missing` when simulating predictively).
  - `Y_t`: the expected *eventual*-total series.

# Examples
```@example RightTruncate
using ComposableTuringIDModels, Distributions
obs = RightTruncate(NegativeBinomialError(), truncated(Normal(5.0, 2.0), 0.0, Inf))
mdl = as_turing_model(obs, missing, fill(100.0, 30))
rand(mdl)
```

## Fields

  - `model`: the inner observation-error model the corrected expected observations
    are passed to.
  - `cdf_model`: the correction submodel producing the length-`n` reporting
    completeness curve (by age).
"
struct RightTruncate{M <: AbstractObservationModel, C <: AbstractLatentModel} <:
       AbstractObservationModel
    "The inner observation-error model."
    model::M
    "The reporting-completeness correction submodel."
    cdf_model::C
end

function RightTruncate(model::M, distribution::C; D = nothing,
        Δd = 1.0) where {
        M <: AbstractObservationModel, C <: ContinuousDistribution}
    return RightTruncate(model, ReportingCDF(distribution; D = D, Δd = Δd))
end

function RightTruncate(model::M, cdf::V) where {
        M <: AbstractObservationModel, V <: AbstractVector{<:Real}}
    return RightTruncate(model, ReportingCDF(cdf))
end

@model function as_turing_model(obs_model::RightTruncate, y_t, Y_t)
    n = length(Y_t)

    # Draw the reporting completeness `F` (by age) from the correction submodel.
    completeness ~ as_turing_submodel(obs_model.cdf_model, n)
    @assert length(completeness)==n "The reporting-completeness curve must have length $n (the expected-observation series length); got $(length(completeness))"

    # `completeness[a + 1]` is the completeness of a reference day of age `a`. The
    # most recent reference day (`t = n`, age `0`) is least complete and the
    # oldest (`t = 1`, age `n - 1`) most complete, so reverse the age-indexed
    # completeness onto the reference-day axis: `scale[t] = completeness(age = n - t)`.
    scaled_Y_t = Y_t .* reverse(completeness)

    inner ~ as_turing_submodel(obs_model.model, y_t, scaled_Y_t)
    return (; y_t = inner.y_t, expected = inner.expected)
end
