# CDF-scaling nowcasting observation modifier (EpiNow2-style right-truncation
# correction for a fixed reporting delay).

@doc raw"
Scale the expected eventual totals of an underlying observation model by the
reporting-delay CDF, correcting for right-truncation (not-yet-reported counts).

This is the EpiNow2-style **CDF-scaling nowcast**. The infection → expected-observation
pipeline produces ``Y_t = \mu_t``, the expected *eventual* total for reference day
`t`. At the present time `now` (taken to be the last reference day, `now = n`) a
reference day of age ``a = now - t`` has only been partially reported: a fraction
``F[a + 1]`` of its eventual total has been observed, where ``F`` is the
reporting-delay CDF (``F[a+1] = \sum_{d=0}^{a} p[d+1]`` for the delay PMF `p`). The
**expected observed-so-far** is therefore ``\mu_t \cdot F[a+1]``, and conditioning
the inner observation error on that — rather than on the full ``\mu_t`` — is the
whole of CDF-scaling. Recent, still-maturing reference days are automatically
down-weighted, while the model's `Y_t` remains the eventual total (so the nowcast
of the eventual total is just `Y_t` read out as a generated quantity).

The modifier mirrors [`Ascertainment`](@ref) / [`LatentDelay`](@ref): it wraps an
inner observation model, transforms the expected-observation vector, and delegates.
The most recent reference day (`t = n`, age `0`) is scaled by `F[1]`, the oldest
(`t = 1`, age `n - 1`) by `F[n]`, so the scaling is applied by **reversing** the
age-indexed CDF onto the reference-day axis. Reference days older than the delay's
support (age `≥ length(delay_cdf)`) are treated as **fully reported** (scale `1`),
so `delay_cdf` need only cover the maturing tail — a delay PMF/CDF of length
`Dmax` is the natural input, exactly as for [`LatentDelay`](@ref).

This is the **fixed-delay** variant: `delay_cdf` is precomputed once and held
constant (the EpiNow2 default). An estimated / time-varying delay is a planned
follow-up (a later phase of the nowcasting design).

## Constructors

  - `CDFScaledObs(model, delay_cdf)` — from a precomputed reporting-delay CDF
    (non-decreasing, in `[0, 1]`). Reference days older than `length(delay_cdf)`
    are taken as fully reported, so the CDF need only span the maturing tail.
  - `CDFScaledObs(model, distribution; D, Δd)` — discretise a continuous reporting
    delay distribution via double-interval censoring (CensoredDistributions.jl) and
    take the cumulative sum of the resulting PMF, exactly the released-CD path
    [`LatentDelay`](@ref) / [`EpiData`](@ref) use.

# Arguments

  - `obs_model`: the [`CDFScaledObs`](@ref) model.
  - `y_t`: the observed-so-far series (or `missing` when simulating predictively).
  - `Y_t`: the expected *eventual*-total series.

# Examples
```@example CDFScaledObs
using EpiAwarePrototype, Distributions
obs = CDFScaledObs(NegativeBinomialError(), truncated(Normal(5.0, 2.0), 0.0, Inf))
mdl = as_turing_model(obs, missing, fill(100.0, 30))
rand(mdl)
```

## Fields

  - `model`: the inner observation-error model the CDF-scaled expected
    observations are passed to.
  - `delay_cdf`: the reporting-delay CDF `F`, indexed by age (`F[a + 1]` is the
    fraction of an eventual total reported within `a` days); non-decreasing, in
    `[0, 1]`.
"
struct CDFScaledObs{M <: AbstractObservationModel, T <: AbstractVector{<:Real}} <:
       AbstractObservationModel
    "The inner observation-error model."
    model::M
    "The reporting-delay CDF, indexed by age."
    delay_cdf::T

    function CDFScaledObs(model::M,
            delay_cdf::T) where {
            M <: AbstractObservationModel, T <: AbstractVector{<:Real}}
        @assert all(>=(0), delay_cdf) "The delay CDF must be non-negative"
        @assert all(<=(1 + 1e-8), delay_cdf) "The delay CDF must not exceed 1"
        @assert issorted(delay_cdf) "The delay CDF must be non-decreasing"
        return new{M, T}(model, delay_cdf)
    end
end

function CDFScaledObs(model::M, distribution::C; D = nothing,
        Δd = 1.0) where {
        M <: AbstractObservationModel, C <: ContinuousDistribution}
    # Build the reporting-delay CDF from the released-CD double-interval-censored
    # PMF (the same path `LatentDelay` / `EpiData` use), then accumulate it.
    pmf = _discretised_pmf(distribution; Δd = Δd, D = D)
    return CDFScaledObs(model, cumsum(pmf))
end

@model function as_turing_model(obs_model::CDFScaledObs, y_t, Y_t)
    n = length(Y_t)
    F = obs_model.delay_cdf
    nF = length(F)

    # `F[a + 1]` is the completeness of a reference day of age `a`. Reference days
    # older than the delay's support (age ≥ `nF`) are fully reported, so the CDF
    # is conceptually padded with ones up to age `n - 1`. The most recent
    # reference day (`t = n`, age `0`) is least complete and the oldest (`t = 1`,
    # age `n - 1`) most complete, so the age-indexed completeness is reversed onto
    # the reference-day axis: `scale[t] = completeness(age = n - t)`.
    completeness = nF >= n ? F[1:n] : vcat(F, ones(eltype(F), n - nF))
    scale = reverse(completeness)
    scaled_Y_t = Y_t .* scale

    y_t ~ to_submodel(as_turing_model(obs_model.model, y_t, scaled_Y_t), false)
    return y_t
end
