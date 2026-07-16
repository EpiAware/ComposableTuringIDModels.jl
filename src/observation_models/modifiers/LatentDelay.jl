# Reporting-delay observation modifier. Its accumulation step (`LDStep`) lives in
# `src/steps/`.

@doc raw"
Apply a reporting delay to an underlying observation model.

The expected observations are convolved with the (reversed) delay PMF before
being passed to the wrapped `model`. `LatentDelay` shortens the expected
observation vector by the length of the delay PMF to avoid fitting to partially
observed data.

The delay can be **fixed** (a known PMF, or a continuous distribution discretised
once at construction) or **uncertain** (an [`UncertainDelay`](@ref) whose delay
distribution's parameters are prior slots, sampled and rediscretised per draw so
the delay is inferred alongside everything else). The uncertain case threads its
parameters through the same [`as_turing_submodel`](@ref) seam every other
component uses, so a prior on the delay composes exactly like a prior anywhere
else in the model.

## Constructors
  - `LatentDelay(model, pmf)` — from a fixed delay PMF (non-negative, sums to 1).
  - `LatentDelay(model, distribution; D, Δd)` — discretise a fixed continuous
    delay distribution once via double-interval censoring
    (CensoredDistributions.jl).
  - `LatentDelay(model, delay::UncertainDelay)` — an inferred delay whose
    distribution parameters carry priors (see [`UncertainDelay`](@ref)).

## Fields

  - `model`: the wrapped observation model the delayed expected observations are
    passed to.
  - `delay`: the delay specification — either the **reversed** fixed delay PMF (a
    vector), or an [`UncertainDelay`](@ref) component that samples the delay
    parameters and builds the PMF per draw.

# Examples
```@example LatentDelay
using ComposableTuringIDModels, Distributions
obs = LatentDelay(NegativeBinomialError(), truncated(Normal(5.0, 2.0), 0.0, Inf))
mdl = as_turing_model(obs, missing, fill(10, 30))
mdl()
```
"
struct LatentDelay{M <: AbstractObservationModel, D} <: AbstractObservationModel
    model::M
    delay::D

    function LatentDelay{M, D}(model, delay) where {M, D}
        return new{M, D}(model, delay)
    end
end

# Fixed PMF: validate and store it reversed for the `LDStep` convolution (the
# fixed path is behaviourally unchanged from the original modifier).
function LatentDelay(model::AbstractObservationModel, pmf::AbstractVector{<:Real})
    @assert all(>=(0), pmf) "Delay PMF must be non-negative"
    @assert isapprox(sum(pmf), 1) "Delay PMF must sum to 1"
    rev_pmf = reverse(pmf)
    return LatentDelay{typeof(model), typeof(rev_pmf)}(model, rev_pmf)
end

# Fixed continuous distribution: discretise once at construction.
function LatentDelay(model::AbstractObservationModel, distribution::C;
        D = nothing, Δd = 1.0) where {C <: ContinuousDistribution}
    pmf = _discretised_pmf(distribution; Δd = Δd, D = D)
    return LatentDelay(model, pmf)
end

# Uncertain / process-valued delay: hold the component and sample it per draw.
function LatentDelay(model::AbstractObservationModel, delay::AbstractPriorModel)
    return LatentDelay{typeof(model), typeof(delay)}(model, delay)
end

@doc raw"
A delay whose distribution's **parameters are prior slots**, so the reporting
delay is inferred rather than fixed.

`UncertainDelay(family, params; D, Δd)` describes a continuous delay distribution
`family(θ...)` whose positional parameters `θ` are drawn from the priors in
`params` (one prior per parameter). Each draw builds the right-truncated,
double-interval-censored delay PMF (the same [`_discretised_pmf`](@ref
ComposableTuringIDModels._discretised_pmf) path the fixed delay uses), so the
delay carries uncertainty and can be recovered from data. It is a prior-role
component: sampling `as_turing_model(u::UncertainDelay)` returns the (forward)
delay PMF, drawing the parameters through the [`as_turing_submodel`](@ref) seam.

The truncation horizon `D` is **required** and fixed: it holds the PMF length
constant across draws (only the parameters vary), which the convolution relies on.
`Δd` is the discretisation bin width.

## Constructor

  - `UncertainDelay(family, params; D, Δd = 1.0)` — `family` is a distribution
    constructor (e.g. `LogNormal`, `Gamma`) called as `family(θ...)`, and `params`
    is a vector of priors, one per positional parameter of `family`.

## Fields

  - `params`: the priors for the delay distribution's positional parameters.
  - `family`: the distribution constructor built from the sampled parameters.
  - `D`: the fixed right-truncation horizon (keeps the PMF length constant).
  - `Δd`: the discretisation bin width.

# Examples
```@example UncertainDelay
using ComposableTuringIDModels, Distributions
delay = UncertainDelay(
    LogNormal, [Normal(1.5, 0.4), truncated(Normal(0.4, 0.2), 0, Inf)]; D = 20.0)
obs = LatentDelay(NegativeBinomialError(), delay)
mdl = as_turing_model(obs, missing, fill(100.0, 40))
rand(mdl)
```
"
struct UncertainDelay{P <: AbstractVector{<:Distribution}, F, T <: Real} <:
       AbstractPriorModel
    params::P
    family::F
    D::T
    Δd::T
end

function UncertainDelay(family, params::AbstractVector{<:Distribution};
        D, Δd = 1.0)
    @assert !isnothing(D) "UncertainDelay needs a fixed horizon `D` so the delay PMF length is constant across draws"
    @assert Δd>0.0 "Δd must be positive"
    Dp, Δdp = promote(float(D), float(Δd))
    @assert Dp>=Δdp "D can't be shorter than Δd"
    return UncertainDelay{typeof(params), typeof(family), typeof(Dp)}(
        params, family, Dp, Δdp)
end

@model function as_turing_model(u::UncertainDelay)
    # Sample the delay distribution's parameters through the priors seam (the AR
    # damp/init pattern: a vector of per-parameter priors drawn as one slot),
    # then rebuild and discretise the delay per draw.
    θ ~ as_turing_submodel(u.params, length(u.params))
    return _discretised_pmf(u.family(θ...); Δd = u.Δd, D = u.D)
end

@model function as_turing_model(obs_model::LatentDelay, y_t, Y_t)
    if ismissing(y_t)
        y_t = Vector{Missing}(missing, length(Y_t))
    end

    # The delay slot yields the reversed delay PMF used by `LDStep`. A fixed PMF
    # vector is already stored reversed and used directly; an uncertain-delay
    # component samples its distribution parameters through the priors seam and
    # builds the (forward) PMF per draw, which is then reversed.
    if obs_model.delay isa AbstractPriorModel
        delay ~ as_turing_submodel(obs_model.delay; prefix = true)
        rev_pmf = reverse(delay)
    else
        rev_pmf = obs_model.delay
    end

    pmf_length = length(rev_pmf)
    @assert pmf_length<=length(Y_t) "The delay PMF must be no longer than the observation vector"

    expected_obs = accumulate_scan(
        LDStep(rev_pmf),
        (; val = 0, current = Y_t[1:pmf_length]),
        vcat(Y_t[(pmf_length + 1):end], 0.0))

    inner ~ as_turing_submodel(obs_model.model, y_t, expected_obs)
    return (; y_t = inner.y_t, expected = inner.expected)
end
