# Reporting-delay observation modifier. Its accumulation step (`LDStep`) lives in
# `src/steps/`.

@doc raw"
Apply a reporting delay to an underlying observation model.

The expected observations are convolved with the (reversed) delay PMF before
being passed to the wrapped `model`. `LatentDelay` shortens the expected
observation vector by the length of the delay PMF to avoid fitting to partially
observed data.

The delay composes through the same constant-vs-process seam as every other
parameter, so it can be **fixed** (a known PMF, or a continuous distribution
discretised once at construction), **uncertain** (an [`UncertainDelay`](@ref)
whose distribution's parameters are prior slots, sampled and rediscretised per
draw so the delay is inferred), or **time-varying** (a delay that changes with
time). A time-varying delay is either a deterministic per-time sequence of PMFs or
an [`UncertainDelay`](@ref) with a process-valued parameter; it is applied with a
time-indexed convolution (one reversed kernel per step). The uncertain and
process-valued cases thread their parameters through the same
[`as_turing_submodel`](@ref) seam every other component uses, so a prior on the
delay composes exactly like a prior anywhere else in the model.

## Constructors
  - `LatentDelay(model, pmf)` — from a fixed delay PMF (non-negative, sums to 1).
  - `LatentDelay(model, distribution; D, Δd)` — discretise a fixed continuous
    delay distribution once via double-interval censoring
    (CensoredDistributions.jl). `D` is the truncation horizon; if omitted and
    `distribution` has finite support (e.g. it was built with `truncated`),
    the support's upper bound is used as `D` directly, so a caller who has
    already truncated the distribution does not need to repeat the horizon
    as a separate keyword — `truncated(dist, lower, upper)` is enough.
  - `LatentDelay(model, pmfs::AbstractVector{<:AbstractVector})` — a deterministic
    time-varying delay from a per-time sequence of PMFs (one per time point, all
    the same length, each non-negative and summing to 1).
  - `LatentDelay(model, delay::UncertainDelay)` — an inferred delay whose
    distribution parameters carry priors; time-invariant when the parameters are
    `Distribution`s and time-varying when any is a process (see
    [`UncertainDelay`](@ref)).

## Fields

  - `model`: the wrapped observation model the delayed expected observations are
    passed to.
  - `delay`: the delay specification — the fixed delay PMF (a vector), a
    per-time sequence of PMFs, or an [`UncertainDelay`](@ref) component that
    samples the delay parameters and builds the PMF(s) per draw. It is held as
    it was given and reversed for the convolution when the model is built.

# Examples

A fixed delay:

```@example LatentDelay
using ComposableTuringIDModels, Distributions
obs = LatentDelay(NegativeBinomialError(), truncated(Normal(5.0, 2.0), 0.0, Inf))
mdl = as_turing_model(obs, missing, fill(10, 30))
mdl()
```

A fixed delay with an explicit horizon expressed on the distribution itself
(truncating above at 15, rather than passing `D = 15.0` separately):

```@example LatentDelay
bounded = LatentDelay(
    NegativeBinomialError(), truncated(Normal(5.0, 2.0), 0.0, 15.0)
)
length(bounded.delay)
```

A deterministic time-varying delay — a PMF per time point (here sharpening over
time), all the same length:

```@example LatentDelay
n = 30
pmfs = [(w = [0.6 - 0.01t, 0.3, 0.1 + 0.01t]; w ./ sum(w)) for t in 1:n]
tv = LatentDelay(PoissonError(), pmfs)
as_turing_model(tv, missing, fill(100.0, n))().y_t
```
"
struct LatentDelay{M <: AbstractObservationModel, D} <: AbstractObservationModel
    model::M
    delay::D

    function LatentDelay{M, D}(model, delay) where {M, D}
        return new{M, D}(model, delay)
    end
end

# Fixed PMF: validate and store it as given, reversed for the `LDStep`
# convolution when the model is built.
# A reversed PMF is still a valid one, so storing it reversed would leave the
# field's stored and given forms indistinguishable and a PMF set into the slot
# would be convolved backwards with nothing to say so.
function LatentDelay(model::AbstractObservationModel, pmf::AbstractVector{<:Real})
    @assert all(>=(0), pmf) "Delay PMF must be non-negative"
    @assert isapprox(sum(pmf), 1) "Delay PMF must sum to 1"
    return LatentDelay{typeof(model), typeof(pmf)}(model, pmf)
end

# Fixed continuous distribution: discretise once at construction.
function LatentDelay(
        model::AbstractObservationModel, distribution::C;
        D = nothing, Δd = 1.0
    ) where {C <: ContinuousDistribution}
    pmf = _discretised_pmf(distribution; Δd = Δd, D = D)
    return LatentDelay(model, pmf)
end

# A deterministic time-varying delay is one PMF per time point, all the same
# length.
# Each is validated and the sequence is stored forward, then reversed per step in
# the `@model` for the `TimeVaryingLDStep`.
function LatentDelay(
        model::AbstractObservationModel,
        pmfs::AbstractVector{<:AbstractVector{<:Real}}
    )
    @assert !isempty(pmfs) "Delay PMF sequence must be non-empty"
    d = length(first(pmfs))
    for pmf in pmfs
        @assert all(>=(0), pmf) "Each delay PMF must be non-negative"
        @assert isapprox(sum(pmf), 1) "Each delay PMF must sum to 1"
        @assert length(pmf) == d "All delay PMFs must have the same length"
    end
    return LatentDelay{typeof(model), typeof(pmfs)}(model, pmfs)
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
`params` (one prior per parameter). Each prior may be a bare `Distribution` (the
parameter is **constant** — uncertain but time-invariant) or a process (an
[`AbstractPriorModel`](@ref), e.g. a `RandomWalk`, so the parameter is
**time-varying**), through the same constant-vs-process seam every other component
uses. Each draw builds the right-truncated, double-interval-censored delay PMF
(the same `_discretised_pmf` path the fixed delay uses):

  - if every parameter is a `Distribution`, one time-invariant PMF is built and
    `as_turing_model(u::UncertainDelay)` returns it (drawing the parameters as one
    slot through the [`as_turing_submodel`](@ref) seam); and
  - if any parameter is a process, `as_turing_model(u::UncertainDelay, n)` builds
    one PMF **per time point** — each parameter is read at time `t` via
    [`at`](@ref) (a constant stays constant, a process path is indexed), so the
    delay, and its discretised PMF, varies with time. A time-varying delay needs a
    series length, so the no-`n` method raises an error.

The truncation horizon `D` is **required** and fixed: it holds the PMF length
constant across draws and across time points (only the parameters vary), which the
convolution relies on. `Δd` is the discretisation bin width.

## Constructor

  - `UncertainDelay(family, params; D, Δd = 1.0)` — `family` is a distribution
    constructor (e.g. `LogNormal`, `Gamma`) called as `family(θ...)`, and `params`
    is a vector of priors (a `Distribution` or a process), one per positional
    parameter of `family`.

## Fields

  - `params`: the priors for the delay distribution's positional parameters.
  - `family`: the distribution constructor built from the sampled parameters.
  - `D`: the fixed right-truncation horizon (keeps the PMF length constant).
  - `Δd`: the discretisation bin width.

# Examples

An uncertain (time-invariant) delay whose `LogNormal` parameters carry priors:

```@example UncertainDelay
using ComposableTuringIDModels, Distributions
delay = UncertainDelay(
    LogNormal, [Normal(1.5, 0.4), truncated(Normal(0.4, 0.2), 0, Inf)]; D = 20.0)
obs = LatentDelay(NegativeBinomialError(), delay)
mdl = as_turing_model(obs, missing, fill(100.0, 40))
rand(mdl)
```

A time-varying delay: the `LogNormal` meanlog is a `RandomWalk` (it drifts with
time) while the sdlog carries a constant prior:

```@example UncertainDelay
tv = UncertainDelay(
    LogNormal, [RandomWalk(), truncated(Normal(0.4, 0.2), 0, Inf)]; D = 20.0)
tv_obs = LatentDelay(NegativeBinomialError(), tv)
as_turing_model(tv_obs, missing, fill(100.0, 40))().y_t
```
"
struct UncertainDelay{P <: AbstractVector, F, T <: Real, TV} <: AbstractPriorModel
    params::P
    family::F
    D::T
    Δd::T

    # The fields in declaration order, which is the form a field-wise rebuild
    # supplies and which the keyword constructor below has no signature for.
    # `TV` is read off `params` here rather than taken, so it cannot come to
    # describe parameters the delay no longer holds.
    function UncertainDelay(params::AbstractVector, family, D, Δd)
        @assert !isnothing(D) "UncertainDelay needs a fixed horizon `D` so the delay PMF length is constant across draws"
        @assert Δd > 0.0 "Δd must be positive"
        @assert all(p -> p isa Distribution || p isa AbstractPriorModel, params) "Each UncertainDelay parameter must be a `Distribution` (constant → uncertain) or a process (an `AbstractPriorModel` → time-varying)"
        Dp, Δdp = promote(float(D), float(Δd))
        @assert Dp >= Δdp "D can't be shorter than Δd"
        # A process-valued parameter makes the delay time-varying; recorded as a
        # type parameter so the `LatentDelay` convolution branch stays
        # type-stable.
        tv = any(p -> p isa AbstractPriorModel, params)
        return new{typeof(params), typeof(family), typeof(Dp), tv}(
            params, family, Dp, Δdp
        )
    end
end

function UncertainDelay(family, params::AbstractVector; D, Δd = 1.0)
    return UncertainDelay(params, family, D, Δd)
end

# An all-constant delay gives a single time-invariant pmf.
# The parameters are sampled through the priors seam as one `θ` slot, then the
# delay is rebuilt and discretised per draw.
@model function as_turing_model(u::UncertainDelay{P, F, T, false}) where {P, F, T}
    θ ~ as_turing_submodel(u.params, length(u.params))
    return _discretised_pmf(u.family(θ...); Δd = u.Δd, D = u.D)
end

@doc raw"
Sample a **constant-parameter** [`UncertainDelay`](@ref) over an axis,
returning the same pmf at every point.

Every parameter is a `Distribution` (uncertain but constant), so one pmf is
drawn — through the no-`n` method above — and copied across the axis. This is
what lets a stratified renewal's generation-interval slot take a fully
constant `UncertainDelay`: drawn at `n_strata`, every stratum gets the same
(uncertain) interval, with no new component. Give one of the parameters a
[`Hierarchy`](@ref) instead for a **partially pooled** per-stratum interval —
that makes the delay time-varying (see the method below) rather than needing
this one.
"
@model function as_turing_model(
        u::UncertainDelay{P, F, T, false}, n::Int
    ) where {P, F, T}
    pmf ~ to_submodel(as_turing_model(u), false)
    return [pmf for _ in 1:n]
end

# A time-varying delay yields a pmf per time point, so it needs a series length and
# cannot stand in for a single time-invariant pmf.
function as_turing_model(u::UncertainDelay{P, F, T, true}) where {P, F, T}
    throw(
        ArgumentError(
            "A time-varying UncertainDelay (a process-valued parameter) yields a " *
                "pmf per time point and needs a series length `n`; it cannot be used " *
                "where a single time-invariant pmf is required (e.g. a Renewal " *
                "generation interval)."
        )
    )
end

@doc raw"
Sample a **time-varying** [`UncertainDelay`](@ref) as a length-`n` sequence of
delay pmfs.

Each parameter is drawn through the [`as_turing_submodel`](@ref) seam: a
`Distribution` parameter draws a scalar (constant across time), while a process
parameter (an [`AbstractPriorModel`](@ref)) draws a length-`n` path. The pmf at
time `t` is built from each parameter read at `t` via [`at`](@ref), so the delay
distribution — and its discretised pmf — varies with time. The fixed horizon `D`
keeps every pmf the same length.
"
@model function as_turing_model(u::UncertainDelay{P, F, T, true}, n::Int) where {
        P, F, T,
    }
    np = length(u.params)
    # Each parameter's submodel is drawn under its own explicit prefix so a
    # single `~` left-hand side holds the return.
    # `as_turing_model(prior, n)` is used directly rather than through the seam's
    # scalar short-circuit, so a constant parameter is a prefixable submodel too.
    params = Vector{Any}(undef, np)
    for i in 1:np
        drawn ~ to_submodel(
            prefix(as_turing_model(u.params[i], n), Symbol(:param, i)), false
        )
        params[i] = drawn
    end
    # Freeze into a `Tuple`, because a `Vector{Any}` keeps every element boxed
    # and the draws then reach `at` as `Base.RefValue`s it cannot index.
    params_t = Tuple(params)
    return map(1:n) do t
        # `_at_all` reads each element at a static position, because indexing a
        # heterogeneous tuple with a runtime `i` would return a `Union` that
        # Enzyme's type analysis rejects.
        θs = _at_all(params_t, t)
        _discretised_pmf(u.family(θs...); Δd = u.Δd, D = u.D)
    end
end

# Read every drawn parameter at time `t`, one static tuple position at a time.
_at_all(::Tuple{}, t) = ()
function _at_all(params::Tuple, t)
    return (at(first(params), t), _at_all(Base.tail(params), t)...)
end

# Whether a `delay` field yields per-time kernels or a single time-invariant pmf.
# Dispatched on the field type so the `LatentDelay` convolution branch is
# type-stable, with no runtime `isa` on sampled values.
_delay_timevarying(::AbstractVector{<:Real}) = false
_delay_timevarying(::AbstractVector{<:AbstractVector{<:Real}}) = true
_delay_timevarying(::AbstractPriorModel) = false
_delay_timevarying(::UncertainDelay{P, F, T, TV}) where {P, F, T, TV} = TV

# A fixed PMF is stored as it was given and the convolution wants it reversed,
# so it is reversed here, once when the model is built, rather than in the
# `@model` body where it would land on the differentiated path.
# The other specifications build their PMFs per draw and reverse them there.
_reversed_fixed_delay(spec::AbstractVector{<:Real}) = reverse(spec)
_reversed_fixed_delay(spec) = spec

function as_turing_model(obs_model::LatentDelay, y_t, Y_t)
    return _latent_delay(
        obs_model, _reversed_fixed_delay(obs_model.delay), y_t, Y_t
    )
end

@model function _latent_delay(obs_model::LatentDelay, spec, y_t, Y_t)
    # A `missing` series is passed down as it is, so the component that scores
    # it sizes it from the convolved series it actually reads.
    # Sizing it here instead would simulate a series as long as the input, with
    # the lead-in at its head coming back as `missing` rather than not coming
    # back at all.
    #
    # An aggregation is the exception, because its reporting windows are read off
    # the length of `y_t` and it right-aligns a shorter `Y_t` against them.
    # It needs the series on the input calendar, which sizing it here supplies.
    if ismissing(y_t) && _needs_input_calendar(obs_model.model)
        y_t = Vector{Missing}(missing, length(Y_t))
    end
    n = length(Y_t)

    # `_delay_timevarying` is a compile-time-constant trait on the `delay` field
    # type, so this branch is type-stable and the fast time-invariant path is
    # untouched.
    if _delay_timevarying(spec)
        # A reversed kernel per time drives `TimeVaryingLDStep`, with the kernel
        # at time `t` applied to the window ending at `t`.
        if spec isa AbstractPriorModel
            delay ~ as_turing_submodel(spec, n; prefix = true)
            pmfs = delay
        else
            pmfs = spec
        end
        @assert length(pmfs) == n "A per-time delay needs one PMF per time point (length must equal the observation vector)"
        rev_pmfs = reverse.(pmfs)
        d = length(first(rev_pmfs))
        @assert d <= n "The delay PMF must be no longer than the observation vector"
        expected_obs = accumulate_scan(
            TimeVaryingLDStep(),
            (; val = 0, current = Y_t[1:d]),
            collect(zip(vcat(Y_t[(d + 1):end], 0.0), rev_pmfs[d:end]))
        )
    else
        # A single reversed pmf drives `LDStep`, which is the fast path.
        # A fixed PMF arrives already reversed; an all-constant uncertain delay
        # builds its forward PMF per draw, which is reversed here.
        if spec isa AbstractPriorModel
            delay ~ as_turing_submodel(spec; prefix = true)
            rev_pmf = reverse(delay)
        else
            rev_pmf = spec
        end
        d = length(rev_pmf)
        @assert d <= n "The delay PMF must be no longer than the observation vector"
        expected_obs = accumulate_scan(
            LDStep(rev_pmf),
            (; val = 0, current = Y_t[1:d]),
            vcat(Y_t[(d + 1):end], 0.0)
        )
    end

    inner ~ as_turing_submodel(obs_model.model, y_t, expected_obs)
    return (; y_t = inner.y_t, expected = inner.expected)
end
