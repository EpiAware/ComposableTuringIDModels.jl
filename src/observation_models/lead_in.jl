# How much of the head of a series an observation chain drops, and how many
# observations are left to score. Structural accessors over an assembled
# observation model: no sampling, no data needed beyond its length.

@doc raw"
The number of leading time points an observation model drops before it scores
anything — the chain's **lead-in**.

Each [`LatentDelay`](@ref) convolves the expected series with a delay PMF and
returns a series shorter by `length(pmf) - 1`, because the head of a convolution
is only partially observed. The observation-error loop then right-aligns the
data against what is left, so with `as_turing_model(model, y, length(y))` the
first `observation_lead_in(model)` entries of `y` are never scored.

Ask for the lead-in and add it to the series length instead:

```julia
n = length(y) + observation_lead_in(model)
mdl = as_turing_model(model, y, n)
```

`n` stays what it has always been — the length of the infection series — and the
chain now hands the error model exactly `length(y)` expected values.
[`observation_coverage`](@ref) reports the same arithmetic against a specific
`n`.

The walk is structural: it recurses through whichever fields of a component are
themselves observation models, so a modifier that preserves the series length
([`Ascertainment`](@ref), [`RightTruncate`](@ref),
[`TransformObservationModel`](@ref), [`ReportTriangle`](@ref),
[`PrefixObservationModel`](@ref), [`RecordExpectedObs`](@ref)) passes its inner
model's lead-in through unchanged.

## Split streams

[`Split`](@ref) runs its streams in parallel off one expected series, so their
lead-ins do not add up. Streams that drop the same amount report that shared
lead-in as a single number. Streams that drop different amounts (a short
infection→report delay beside a long infection→death one) have no single answer,
so the lead-in comes back as a `NamedTuple` keyed by stream: each stream scores
the last `n - lead_in[stream]` of its own series, and it is the *data* that
differs in length between streams, not `n`.

## Units

A lead-in is counted on the axis the delay sits on. Under an [`Aggregate`](@ref)
that axis is reporting windows rather than time points, so a delay nested inside
an aggregation drops that many **windows**.

## Arguments

  - `model`: an observation model, an [`IDModel`](@ref) or an [`IDProblem`](@ref).

# Examples
```@example observation_lead_in
using ComposableTuringIDModels, Distributions
obs = LatentDelay(LatentDelay(PoissonError(), fill(1 / 15, 15)), fill(1 / 30, 30))
observation_lead_in(obs)
```

Streams with different delays report one lead-in each:

```@example observation_lead_in
streams = Split((
    cases = LatentDelay(PoissonError(), fill(1 / 5, 5)),
    deaths = LatentDelay(PoissonError(), fill(1 / 20, 20))))
observation_lead_in(streams)
```
"
function observation_lead_in end

# The observation model a component wraps, or `nothing` for one that consumes
# the expected series itself. Found by a field walk rather than a method per
# modifier, so a new length-preserving modifier needs no accessor of its own.
# A modifier wraps a single observation model; several in parallel is what a
# `Split` is, and every walk below gives it its own method.
function _wrapped_model(model::AbstractObservationModel)
    for f in fieldnames(typeof(model))
        v = getfield(model, f)
        v isa AbstractObservationModel && return v
    end
    return nothing
end

# Default: a component drops nothing of its own, so its lead-in is that of the
# observation model it wraps, passed straight through.
function observation_lead_in(model::AbstractObservationModel)
    inner = _wrapped_model(model)
    return inner === nothing ? 0 : observation_lead_in(inner)
end

# Add a component's own lead-in to the one it wraps. The wrapped side is a plain
# count until a `Split` with unequal streams reports one per stream; from there
# the surrounding chain's lead-in is added to every stream.
_add_lead_in(a::Int, b::Int) = a + b
_add_lead_in(a::Int, b::NamedTuple) = map(x -> _add_lead_in(a, x), b)

# A delay drops the head of its convolution, on top of whatever it wraps.
function observation_lead_in(model::LatentDelay)
    return _add_lead_in(
        _delay_lead_in(model.delay), observation_lead_in(model.model)
    )
end

# The lead-in of a delay specification: one less than the length of the PMF it
# convolves with. Every supported specification fixes that length before any
# parameter is drawn — a stored (reversed) PMF and a per-time sequence of PMFs
# carry it directly, and an `UncertainDelay` fixes it with its horizon `D` and
# bin width `Δd`, which is exactly why `D` is required.
_delay_lead_in(pmf::AbstractVector{<:Real}) = length(pmf) - 1
_delay_lead_in(pmfs::AbstractVector{<:AbstractVector{<:Real}}) =
    _delay_lead_in(first(pmfs))
_delay_lead_in(u::UncertainDelay) = length(0.0:(u.Δd):(u.D - u.Δd)) - 1

# Anything else is a delay component whose PMF length is not knowable before it
# is sampled. Say so rather than returning a zero that would silently drop
# observations — the failure this accessor exists to prevent.
function _delay_lead_in(spec)
    throw(
        ArgumentError(
            "The lead-in of a $(nameof(typeof(spec))) delay cannot be read " *
                "before it is sampled. Use a fixed PMF, a per-time sequence of " *
                "PMFs, or an `UncertainDelay` (whose horizon `D` fixes the PMF " *
                "length)."
        )
    )
end

# A `Split`'s streams are parallel branches of one expected series, so they share
# a lead-in rather than accumulating one. Streams that drop the same amount
# report it once; streams that drop different amounts report one lead-in each,
# because no single number describes them all. A strata template is one model
# replicated across streams, so it reports its own.
function observation_lead_in(model::Split)
    model.streams isa NamedTuple || return observation_lead_in(model.streams)
    lead_ins = map(observation_lead_in, model.streams)
    return allequal(values(lead_ins)) ? first(values(lead_ins)) : lead_ins
end

# The observation chain a diagnostic runs over. A bare observation model is its
# own chain; the composed wrappers hand theirs over (their methods sit with the
# types, in `compose.jl` and `IDProblem.jl`).
_observation_chain(model::AbstractObservationModel) = model

# The per-stream observation models a chain fans out to, or `nothing` when it
# scores a single series. Follows the same structural path as
# `observation_lead_in`, so a coverage report is keyed by the same streams and
# each stream's data is counted by the model that actually scores it.
_stream_models(m::AbstractObservationModel) = _stream_models(_wrapped_model(m))

# The walk reached the component that consumes the expected series without
# meeting a `Split`, so there is one stream and nothing to key a report by.
_stream_models(::Nothing) = nothing

# A strata template is one model replicated across streams the data names, so
# there are no per-stream models to key a report by either.
_stream_models(m::Split) = m.streams isa NamedTuple ? m.streams : nothing

@doc raw"
How many of the supplied observations an observation chain actually scores, for
a given series length `n`.

Returns `(; n_observations, lead_in, n_scored, n_unscored)`:

  - `n_observations`: the number of observations supplied.
  - `lead_in`: the chain's lead-in (see [`observation_lead_in`](@ref)).
  - `n_scored`: how many of them enter the likelihood.
  - `n_unscored`: how many are dropped from the head of the series.

The observation-error loop right-aligns the data against the expected series,
so the *last* `n - lead_in` observations are scored and any earlier ones are
silently dropped. Call this before sampling to check `n_unscored` is what you
meant it to be; `n = length(y) + observation_lead_in(model)` makes it zero.

A [`Split`](@ref)'s streams are reported one by one — a `NamedTuple` of the
above, keyed by stream. Each stream is counted by the model that scores it, so a
stream taking `NamedTuple` data is not mistaken for a further split.

An [`IDProblem`](@ref) reads its series length from `tspan`, so it takes the
data on its own: `observation_coverage(problem, data)`.

## What counts as an observation

The count is read from the data through the contract of the model that scores
it, on the time axis:

  - a plain vector, or a [`MissingObservations`](@ref) carrier, is the series
    itself; a `strata × time` matrix is counted along its columns;
  - a model taking extra per-time data reads its observations from the `y` field
    of the supplied `NamedTuple`, so [`BinomialError`](@ref)'s
    `(y = successes, N = trials)` has `length(y)` observations, not two;
  - a [`ReportTriangle`](@ref) is counted down the reference days of its
    reference-day × delay matrix, not across its delay columns;
  - a `missing` is simulated at whatever length the chain produces, so it is
    scored in full.

## Reporting triangles do not right-align

Only the per-time-point error families right-align, so only there does a
non-zero `n_unscored` mean observations quietly dropped from the head.
[`ReportTriangle`](@ref) instead asserts that its triangle has exactly
`n - lead_in` reference days, so a non-zero `n_unscored` there is a call that
will *fail*, and the report says so before it does.

## Forecasting

[`forecast`](@ref) rebuilds the model at a longer series, and needs the same `n`
the fit used: pass `forecast(model, y, chain, horizon; n = n)` when `n` is not
`length(y)`. An [`IDProblem`](@ref) records it in `tspan` and needs no keyword.

## Arguments

  - `model`: an observation model, an [`IDModel`](@ref) or an [`IDProblem`](@ref).
  - `y_t`: the observations (or `missing` when simulating predictively).
  - `n`: the series length passed to `as_turing_model` (an `Int`, or the
    `(n_strata, n_time)` shape of a stratified model).

# Examples
```@example observation_coverage
using ComposableTuringIDModels, Distributions
obs = LatentDelay(PoissonError(), fill(1 / 10, 10))
y = fill(10, 30)
observation_coverage(obs, y, length(y))                          # 9 dropped
observation_coverage(obs, y, length(y) + observation_lead_in(obs))
```
"
function observation_coverage(model, y_t, n::ModelShape)
    chain = _observation_chain(model)
    return _coverage(chain, observation_lead_in(chain), y_t, _n_time(n))
end

# A chain that fans out to named streams reports one coverage per stream, each
# read through that stream's own model and data. `lead_in` is shared when the
# streams drop the same amount and keyed by stream when they do not, so it is
# split by the same `_stream` lookup as the data.
function _coverage(model, lead_in, y_t, n_time::Int)
    streams = _stream_models(model)
    streams === nothing &&
        return _series_coverage(model, lead_in, y_t, n_time)
    return NamedTuple{keys(streams)}(
        map(
            k -> _series_coverage(
                streams[k], _stream(lead_in, k), _stream(y_t, k), n_time
            ),
            keys(streams)
        )
    )
end

# One chain, one series: the last `n - lead_in` observations are scored.
function _series_coverage(model, lead_in::Int, y_t, n_time::Int)
    n_expected = max(n_time - lead_in, 0)
    n_observations = _n_observations(model, y_t, n_expected)
    n_scored = min(n_observations, n_expected)
    return (;
        n_observations, lead_in, n_scored,
        n_unscored = n_observations - n_scored,
    )
end

# Pick a stream's share of a per-stream value; a value shared across streams (a
# common lead-in, or one series reaching every stream) passes through unchanged.
_stream(x::NamedTuple, k::Symbol) = x[k]
_stream(x, ::Symbol) = x

# The number of observations a data argument carries on the time axis, read
# through the contract of the model that scores it. The data reaches the end of
# the chain untouched — a modifier reshapes the expected series, not the
# observations — so the walk descends to the component that consumes it.
_n_observations(model::AbstractObservationModel, y_t, n_expected) =
    _n_observations(_wrapped_model(model), y_t, n_expected)

# The end of the walk: the component that consumes the data reads it itself.
_n_observations(::Nothing, y_t, n_expected) =
    _n_series_observations(y_t, n_expected)

_n_series_observations(y_t::AbstractVector, n_expected) = length(y_t)
_n_series_observations(y_t::AbstractMatrix, n_expected) = size(y_t, 2)
_n_series_observations(y_t::MissingObservations, n_expected) = length(y_t.value)
_n_series_observations(::Missing, n_expected) = n_expected
# A model needing more than the counts takes them in the `y` field of a
# `NamedTuple` (see `define_y_t`); the other fields are known per-time-point
# covariates — `BinomialError`'s number of trials — not observations.
_n_series_observations(y_t::NamedTuple, n_expected) =
    _n_series_observations(y_t.y, n_expected)

# A `ReportTriangle` scores the cells of a reference-day × delay matrix, so its
# observations run down the ROWS: the delay columns are one reference day's
# report split by delay, not further time points.
_n_observations(::ReportTriangle, y_t, n_expected) =
    _n_triangle_observations(y_t, n_expected)

_n_triangle_observations(y_t::ReportingTriangle, n_expected) = size(
    y_t.counts, 1
)
_n_triangle_observations(y_t::AbstractMatrix, n_expected) = size(y_t, 1)
# `missing` (simulating) and a long-form `(reference, delay, count)` table are
# both sized from the expected series, so they cover it exactly.
_n_triangle_observations(y_t, n_expected) = n_expected
