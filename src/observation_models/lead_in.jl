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

# Default: a component drops nothing of its own, so its lead-in is that of the
# observation models it wraps. Written as a field walk rather than a method per
# modifier so a new length-preserving modifier is handled without one.
function observation_lead_in(model::AbstractObservationModel)
    total = 0
    for f in fieldnames(typeof(model))
        v = getfield(model, f)
        if v isa AbstractObservationModel
            total = _add_lead_in(total, observation_lead_in(v))
        elseif v isa Union{AbstractVector, Tuple, NamedTuple}
            for el in v
                el isa AbstractObservationModel &&
                    (total = _add_lead_in(total, observation_lead_in(el)))
            end
        end
    end
    return total
end

# Accumulate lead-ins down a chain. Both sides are plain counts until a `Split`
# with unequal streams reports one per stream; from there the surrounding chain's
# lead-in is added to every stream. Two per-stream lead-ins never meet, because a
# modifier wraps one observation model and a `Split` handles its own streams.
_add_lead_in(a::Int, b::Int) = a + b
_add_lead_in(a::Int, b::NamedTuple) = map(x -> _add_lead_in(a, x), b)
_add_lead_in(a::NamedTuple, b::Int) = map(x -> _add_lead_in(x, b), a)

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
function _delay_lead_in(pmfs::AbstractVector{<:AbstractVector{<:Real}})
    return length(first(pmfs)) - 1
end
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
above, keyed by stream — whenever the streams differ in lead-in or the data
arrives as a `NamedTuple` of per-stream series.

An [`IDProblem`](@ref) reads its series length from `tspan`, so it takes the
data on its own: `observation_coverage(problem, data)`.

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
    return _coverage(observation_lead_in(model), y_t, _n_time(n))
end

# One chain, one series: the last `n - lead_in` observations are scored.
function _coverage(lead_in::Int, y_t, n_time::Int)
    n_expected = max(n_time - lead_in, 0)
    n_observations = _n_observations(y_t, n_expected)
    n_scored = min(n_observations, n_expected)
    return (;
        n_observations, lead_in, n_scored,
        n_unscored = n_observations - n_scored,
    )
end

# Per-stream data under a shared lead-in: one report per stream, since the
# streams need not be the same length.
function _coverage(lead_in::Int, y_t::NamedTuple, n_time::Int)
    return map(y -> _coverage(lead_in, y, n_time), y_t)
end

# Per-stream lead-ins: one report per stream, each against that stream's data
# (a shared `missing` or vector reaches every stream unchanged).
_coverage(lead_in::NamedTuple, y_t, n_time::Int) = _stream_coverage(
    lead_in, y_t, n_time
)
_coverage(lead_in::NamedTuple, y_t::NamedTuple, n_time::Int) = _stream_coverage(
    lead_in, y_t, n_time
)

function _stream_coverage(lead_in::NamedTuple, y_t, n_time::Int)
    return NamedTuple{keys(lead_in)}(
        map(k -> _coverage(lead_in[k], _stream(y_t, k), n_time), keys(lead_in))
    )
end

_stream(y_t::NamedTuple, k::Symbol) = y_t[k]
_stream(y_t, ::Symbol) = y_t

# The number of observations a data argument carries, on the time axis. A
# `missing` is simulated at whatever length the chain produces, so it is scored
# in full.
_n_observations(y_t::AbstractVector, n_expected) = length(y_t)
_n_observations(y_t::AbstractMatrix, n_expected) = size(y_t, 2)
_n_observations(y_t::MissingObservations, n_expected) = length(y_t.value)
_n_observations(::Missing, n_expected) = n_expected
