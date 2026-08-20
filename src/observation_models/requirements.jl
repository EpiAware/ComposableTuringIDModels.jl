# What data an observation chain needs, and how much of the infection series its
# delays consume before it scores anything. Structural accessors over an
# assembled observation model: no sampling, and no data beyond its shape.

@doc raw"
The number of leading time points an observation model consumes before it scores
anything — the chain's **lead-in**.

Each [`LatentDelay`](@ref) convolves the expected series with a delay PMF and
returns a series shorter by `length(pmf) - 1`, because the head of a convolution
is only partially observed. The lead-in is what those delays cost, summed down
the chain.

[`as_turing_model`](@ref) adds it to the series it builds, so `n` is the number
of observations and a caller never does this arithmetic. Reading it back is for
saying *why* a model runs its infection process over more time points than the
data covers:

```julia
mdl = as_turing_model(model, y, length(y))   # infections run for
                                             # length(y) + observation_lead_in(model)
```

The walk is structural, over [`wrapped_models`](@ref): a modifier that preserves
the series length ([`Ascertainment`](@ref), [`RightTruncate`](@ref),
[`TransformObservationModel`](@ref), [`ReportTriangle`](@ref),
[`PrefixObservationModel`](@ref), [`RecordExpectedObs`](@ref)) passes its inner
model's lead-in through unchanged.

## Split streams

[`Split`](@ref) runs its streams in parallel off one expected series, so their
lead-ins do not add up. Streams that consume the same amount report that shared
lead-in as a single number. Streams that consume different amounts (a short
infection→report delay beside a long infection→death one) have no single answer,
so the lead-in comes back as a `NamedTuple` keyed by stream. The series is then
built to cover the deepest of them.

## Aggregations

A lead-in is counted on the axis the delay sits on, and an [`Aggregate`](@ref)
moves that axis onto reporting windows. A delay nested *inside* an aggregation
therefore consumes windows rather than time points, which is not a series
length, and this raises rather than converting between the two. Apply the delay
outside the aggregation instead.

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

# The single observation model a component wraps, or `nothing` for one that
# consumes the expected series itself. Read through the traversal seam, so a new
# length-preserving modifier needs no accessor of its own. A `Split` wraps
# several models in parallel rather than one; it is not a link in a chain, and
# every walk below gives it its own method.
_wrapped_model(model::AbstractObservationModel) = _single_wrapped(
    wrapped_models(model)
)

_single_wrapped(::Tuple{}) = nothing
_single_wrapped(wrapped::Tuple{Any}) = only(wrapped)
_single_wrapped(::Tuple) = nothing

# Default: a component consumes nothing of its own, so its lead-in is that of
# the observation model it wraps, passed straight through.
function observation_lead_in(model::AbstractObservationModel)
    inner = _wrapped_model(model)
    return inner === nothing ? 0 : observation_lead_in(inner)
end

# Add a component's own lead-in to the one it wraps. The wrapped side is a plain
# count until a `Split` with unequal streams reports one per stream; from there
# the surrounding chain's lead-in is added to every stream.
_add_lead_in(a::Int, b::Int) = a + b
_add_lead_in(a::Int, b::NamedTuple) = map(x -> _add_lead_in(a, x), b)

# A delay consumes the head of its convolution, on top of whatever it wraps.
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
# a lead-in rather than accumulating one. Streams that consume the same amount
# report it once; streams that consume different amounts report one lead-in each,
# because no single number describes them all. A strata template is one model
# replicated across streams, so it reports its own.
function observation_lead_in(model::Split)
    model.streams isa NamedTuple || return observation_lead_in(model.streams)
    lead_ins = map(observation_lead_in, model.streams)
    return allequal(values(lead_ins)) ? first(values(lead_ins)) : lead_ins
end

# An aggregation moves the axis: what it wraps consumes reporting windows rather
# than time points, and a window is not a time point, so a delay nested inside
# one has no lead-in expressible as a series length. Say so rather than adding
# windows to a count of days. (A delay outside an aggregation is on the time
# axis and is the ordinary case.)
function observation_lead_in(model::Aggregate)
    inner = observation_lead_in(model.model)
    inner == 0 && return 0
    throw(
        ArgumentError(
            "A delay nested inside an `Aggregate` consumes reporting windows, " *
                "not time points, so the series length it implies cannot be " *
                "derived. Apply the delay outside the aggregation: " *
                "`LatentDelay(Aggregate(model, windows), pmf)`."
        )
    )
end

# The lead-in the infection series has to cover. A chain consumes its own lead-in
# from the head of the series; a `Split`'s streams consume theirs in parallel, so
# the series must be long enough for the deepest of them.
_series_lead_in(model) = _max_lead_in(observation_lead_in(model))

_max_lead_in(lead_in::Int) = lead_in
_max_lead_in(lead_in::NamedTuple) = maximum(values(lead_in))

# The observation chain a requirements report runs over. A bare observation
# model is its own chain; the composed wrappers hand theirs over (their methods
# sit with the types, in `compose.jl` and `IDProblem.jl`).
_observation_chain(model::AbstractObservationModel) = model

# The per-stream observation models a chain fans out to, or `nothing` when it
# scores a single series. Follows the same structural path as
# `observation_lead_in`, so a report is keyed by the same streams and each
# stream's data is read through the model that actually scores it.
_stream_models(m::AbstractObservationModel) = _stream_models(_wrapped_model(m))

# The walk reached the component that consumes the expected series without
# meeting a `Split`, so there is one stream and nothing to key a report by.
_stream_models(::Nothing) = nothing

# A strata template is one model replicated across streams the data names, so
# there are no per-stream models to key a report by either.
_stream_models(m::Split) = m.streams isa NamedTuple ? m.streams : nothing

@doc raw"
What a single observation stream requires of the data supplied for it.

## Fields

  - `name`: the stream's name, or `:y_t` for a model scoring a single series.
  - `n_required`: how long the supplied series should be along its own axis —
    time points, or reference days for a [`ReportTriangle`](@ref). This is the
    `n` the model was asked for.
  - `n_max`: the most this stream can score. It exceeds `n_required` only for a
    [`Split`](@ref) stream whose chain consumes a shorter lead-in than its
    neighbours: one series serves every stream, so it is long enough for the
    deepest, and a shallower stream is handed extra leading expected values it
    can score if earlier data happens to exist.
  - `n_scored`: how many of the expected values enter the likelihood. Equal to
    `n_max` except under an [`Aggregate`](@ref), which scores one value per
    reporting window and ignores the days between.
  - `lead_in`: the time points this stream's chain consumes (see
    [`observation_lead_in`](@ref)). Reported for explanation; it is already
    covered by the series [`as_turing_model`](@ref) builds.
  - `shape`: the shape the data takes — `:series` (one entry per time point) or
    `:triangle` (a reference-day × delay matrix).
  - `alignment`: `:right` when the stream right-aligns its data against its
    expected series, so a shorter series is scored at the end and the earlier
    expected values are unobserved; `:exact` when it does not
    ([`Aggregate`](@ref), [`ReportTriangle`](@ref)), so the lengths must match.
  - `fields`: the fields the scoring model reads when it needs more than the
    counts, e.g. `(:y, :N)` for [`BinomialError`](@ref); empty when a plain
    series is enough.
  - `n_supplied`: how much data was actually supplied, or `nothing` when the
    requirements were asked for without data.
"
struct StreamRequirement
    "The stream's name, or `:y_t` for a single series."
    name::Symbol
    "The length the supplied series should have along its own axis."
    n_required::Int
    "The most this stream can score."
    n_max::Int
    "How many of the expected values enter the likelihood."
    n_scored::Int
    "The time points this stream's chain consumes."
    lead_in::Int
    "The shape the data takes: `:series` or `:triangle`."
    shape::Symbol
    "Whether the stream right-aligns its data (`:right`) or not (`:exact`)."
    alignment::Symbol
    "The fields the scoring model reads, when it needs more than the counts."
    fields::Tuple{Vararg{Symbol}}
    "How much data was supplied, or `nothing` if none was."
    n_supplied::Union{Int, Nothing}
end

# A stream is satisfied when every observation supplied for it is scored. A
# right-aligning stream scores the last `n_max` it is given, so anything up to
# that is covered and more than that would drop the head. A stream that does not
# right-align has to match exactly. Asking without data satisfies nothing and
# fails nothing.
function _stream_fits(r::StreamRequirement)
    r.n_supplied === nothing && return true
    r.alignment === :exact && return r.n_supplied == r.n_max
    return r.n_supplied <= r.n_max
end

@doc raw"
What data a model needs, stream by stream.

## Fields

  - `n`: the number of observations the model was asked for.
  - `series_length`: the length of the infection series built internally —
    `n` plus the deepest chain's lead-in.
  - `streams`: one [`StreamRequirement`](@ref) per observation stream.

Indexing by stream name returns that stream's requirement, and iterating runs
over the streams in order.
"
struct DataRequirements
    "The number of observations the model was asked for."
    n::Int
    "The length of the infection series built internally."
    series_length::Int
    "One requirement per observation stream."
    streams::Vector{StreamRequirement}
end

Base.getindex(r::DataRequirements, name::Symbol) = _stream_named(r, name)
Base.getindex(r::DataRequirements, i::Integer) = r.streams[i]
Base.iterate(r::DataRequirements, state...) = iterate(r.streams, state...)
Base.length(r::DataRequirements) = length(r.streams)
Base.eltype(::Type{DataRequirements}) = StreamRequirement

function _stream_named(r::DataRequirements, name::Symbol)
    i = findfirst(s -> s.name === name, r.streams)
    i === nothing && throw(
        KeyError(name)
    )
    return r.streams[i]
end

@doc raw"
Whether a dataset fits a model: a straight yes or no.

`true` when every observation supplied would be scored. A stream that
right-aligns its data (any of the per-time-point error families) is satisfied by
anything up to the number of expected values it is handed, because a shorter
series is scored at the end and the earlier expected values are simply
unobserved. More than that would leave the head of the data unscored, which is
`false`. An [`Aggregate`](@ref) and a [`ReportTriangle`](@ref) do not
right-align at all, so there the length has to match exactly.

Data that is not there to disagree — a `missing` stream, simulated at whatever
length the chain produces — fits by construction.

[`as_turing_model`](@ref) applies one further check this does not: a stream
short by *exactly* its chain's lead-in is the pre-0.2.0 meaning of `n`, where
the caller added the lead-in to the series length by hand, and is rejected
rather than fitted over a longer series than the caller means.

Takes a model and its data, or a report already built by
[`data_requirements`](@ref). When the answer is no, print the report to see
which stream is wrong and by how much.

## Arguments

  - `model`: an observation model, an [`IDModel`](@ref) or an [`IDProblem`](@ref).
  - `y_t`: the data.
  - `n`: the number of observations.

# Examples
```@example data_fits
using ComposableTuringIDModels
obs = LatentDelay(PoissonError(), fill(1 / 10, 10))
data_fits(obs, fill(10, 30), 30), data_fits(obs, fill(10, 40), 30)
```
"
data_fits(r::DataRequirements) = all(_stream_fits, r.streams)
data_fits(r::StreamRequirement) = _stream_fits(r)
data_fits(model, y_t, n::ModelShape) = data_fits(
    data_requirements(model, y_t, n)
)

@doc raw"
What a caller must supply to fit a model to `n` observations.

Returns a [`DataRequirements`](@ref): one [`StreamRequirement`](@ref) per
observation stream, naming the stream, the length its series must have, how many
of those entries are scored, the lead-in its chain consumes, and the shape the
data takes. Pass the data as well and each stream also reports how much was
supplied, which is what [`data_fits`](@ref) reads to answer whether a dataset
fits a model.

Printing it is the intended way to read it; the fields are there for a front end
that builds its own prompts from them.

## What each stream asks for

  - The **length** is `n` for every stream. [`as_turing_model`](@ref) builds the
    infection series long enough to cover the chain's lead-in, so `n` is the
    number of observations and nothing is dropped from the head.
  - **Scored** differs from the length only under an [`Aggregate`](@ref), which
    sums the expected series over reporting windows and scores one value per
    window: a daily series under weekly reporting is supplied in full and scored
    once a week.
  - The **shape** is a plain series unless the stream ends in a
    [`ReportTriangle`](@ref), which takes a reference-day × delay matrix and
    counts its requirement down the reference days.
  - **Fields** are named when the scoring model needs more than the counts:
    [`BinomialError`](@ref) reads `(y = successes, N = trials)`, so the stream is
    supplied as a `NamedTuple` rather than as a bare vector.

## Arguments

  - `model`: an observation model, an [`IDModel`](@ref) or an [`IDProblem`](@ref).
  - `y_t`: (optional) the data, to report what was supplied alongside what is
    required.
  - `n`: the number of observations (an `Int`, or the `(n_strata, n_time)` shape
    of a stratified model). An [`IDProblem`](@ref) reads it from `tspan`.

# Examples
```@example data_requirements
using ComposableTuringIDModels, Distributions
model = IDModel(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
    Split((
        cases = LatentDelay(PoissonError(), fill(1 / 5, 5)),
        deaths = LatentDelay(NegativeBinomialError(), fill(1 / 20, 20)))))
data_requirements(model, 100)
```
"
function data_requirements(model, n::ModelShape)
    return _requirements(model, missing, _n_time(n))
end

function data_requirements(model, y_t, n::ModelShape)
    return _requirements(model, y_t, _n_time(n))
end

function _requirements(model, y_t, n::Int)
    chain = _observation_chain(model)
    lead_in = observation_lead_in(chain)
    deepest = _max_lead_in(lead_in)
    streams = _stream_models(chain)
    entries = if streams === nothing
        [_stream_requirement(:y_t, chain, lead_in, deepest, y_t, n)]
    else
        [
            _stream_requirement(
                    k, streams[k], _stream(lead_in, k), deepest,
                    _stream(y_t, k), n
                ) for k in keys(streams)
        ]
    end
    return DataRequirements(n, n + deepest, entries)
end

function _stream_requirement(
        name::Symbol, model, lead_in::Int, deepest::Int, y_t, n::Int
    )
    contract = _data_contract(_consumer(model))
    # One series serves every stream, so it covers the deepest lead-in and a
    # shallower stream is handed that much more than `n`.
    n_max = n + deepest - lead_in
    return StreamRequirement(
        name, n, n_max, _n_scored(model, n_max), lead_in, contract.shape,
        _alignment(model), contract.fields,
        _n_supplied(contract.shape, y_t)
    )
end

# Whether a stream right-aligns its data against its expected series. The
# per-time-point error families do, so a shorter series is scored at the end and
# the earlier expected values are unobserved run-in. An `Aggregate` indexes its
# data by a presence mask over the expected series, and a `ReportTriangle`
# asserts its reference days, so neither tolerates a length that differs.
_alignment(model::AbstractObservationModel) = _alignment(_wrapped_model(model))
_alignment(::Nothing) = :right
_alignment(::Aggregate) = :exact
_alignment(::ReportTriangle) = :exact

# Pick a stream's share of a per-stream value; a value shared across streams (a
# common lead-in, or one series reaching every stream) passes through unchanged.
_stream(x::NamedTuple, k::Symbol) = x[k]
_stream(x, ::Symbol) = x

# The component that actually reads the data. The data reaches the end of a
# chain untouched — a modifier reshapes the expected series, not the
# observations — so the walk descends to whatever consumes it.
_consumer(model::AbstractObservationModel) = _descend_consumer(
    model, _wrapped_model(model)
)

_descend_consumer(model, ::Nothing) = model
_descend_consumer(_, inner) = _consumer(inner)

# A `ReportTriangle` reads the data itself, as a matrix, rather than passing it
# down to the error model it wraps.
_consumer(model::ReportTriangle) = model

# What the consuming component needs of the data supplied for its stream. The
# default is one entry per time point and nothing more.
_data_contract(::AbstractObservationModel) = (shape = :series, fields = ())

# A model needing more than the counts takes them in the `y` field of a
# `NamedTuple` (see `define_y_t`); the other fields are per-time-point
# covariates, not observations.
_data_contract(::BinomialError) = (shape = :series, fields = (:y, :N))

# A reporting triangle is a reference-day × delay matrix, counted down the
# reference days: the delay columns are one reference day's report split by
# delay, not further time points.
_data_contract(::ReportTriangle) = (shape = :triangle, fields = ())

# How many of a stream's `n` supplied entries enter the likelihood. One per time
# point everywhere except an `Aggregate`, which sums the series into reporting
# windows and hands the error model one value per window; the days its presence
# mask leaves out are never scored.
_n_scored(model::AbstractObservationModel, n::Int) = _n_scored(
    _wrapped_model(model), n
)
_n_scored(::Nothing, n::Int) = n
_n_scored(model::Aggregate, n::Int) =
    _n_scored(model.model, count(_present_mask(model, n)))

# The presence mask an `Aggregate` applies over `n` time points, broadcast by
# the same rule the model itself uses so the two agree on which points report.
_present_mask(ag::Aggregate, n) =
    broadcast_rule(RepeatEach(), ag.present, n, length(ag.present))

# How much data was supplied for a stream, on the axis its requirement is
# counted along. `nothing` when there is nothing to count: a `missing` stream is
# simulated at whatever length the chain produces, and so fits by construction.
_n_supplied(shape::Symbol, y_t) = shape === :triangle ?
    _n_triangle_supplied(y_t) : _n_series_supplied(y_t)

_n_series_supplied(y_t::AbstractVector) = length(y_t)
_n_series_supplied(y_t::AbstractMatrix) = size(y_t, 2)
_n_series_supplied(y_t::MissingObservations) = length(y_t.value)
_n_series_supplied(y_t::NamedTuple) = hasproperty(y_t, :y) ?
    _n_series_supplied(y_t.y) : nothing
_n_series_supplied(::Any) = nothing

_n_triangle_supplied(y_t::ReportingTriangle) = size(y_t.counts, 1)
_n_triangle_supplied(y_t::AbstractMatrix) = size(y_t, 1)
# A long-form `(reference, delay, count)` table, and `missing`, are both sized
# from the expected series, so they cover it exactly.
_n_triangle_supplied(::Any) = nothing

# --- printing ---------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", r::DataRequirements)
    println(
        io, "DataRequirements: ", r.n, " observation", r.n == 1 ? "" : "s",
        " per stream"
    )
    print(io, "  infection series built at length ", r.series_length)
    for stream in r.streams
        println(io)
        print(io, "  ")
        _show_stream(io, stream)
    end
    return
end

function Base.show(io::IO, ::MIME"text/plain", s::StreamRequirement)
    return _show_stream(io, s)
end

# Compact forms, for a report nested inside something else (a tuple in a
# docstring example, a vector of reports). The struct's own field dump says
# nothing a reader wants.
function Base.show(io::IO, r::DataRequirements)
    print(io, "DataRequirements(", r.n, " observations, ")
    print(io, length(r.streams), " stream")
    length(r.streams) == 1 || print(io, "s")
    return print(io, ", series ", r.series_length, ")")
end

Base.show(io::IO, s::StreamRequirement) = _show_stream(io, s)

function _show_stream(io::IO, s::StreamRequirement)
    print(io, s.name, ": ", _supply_description(s))
    s.n_max == s.n_required ||
        print(io, " (up to ", s.n_max, " if earlier data exists)")
    s.n_scored == s.n_max ||
        print(io, ", ", s.n_scored, " of them scored")
    s.lead_in == 0 || print(io, ", after a lead-in of ", s.lead_in)
    s.n_supplied === nothing && return
    _stream_fits(s) && return print(io, " — supplied ", s.n_supplied)
    return print(io, " — SUPPLIED ", s.n_supplied)
end

function _supply_description(s::StreamRequirement)
    s.shape === :triangle &&
        return "a matrix with $(s.n_required) reference-day rows"
    isempty(s.fields) && return "$(s.n_required) values"
    named = join(("`$f`" for f in s.fields), ", ")
    return "$(s.n_required) values in each of $named"
end

# --- the check `as_turing_model` runs --------------------------------------

# What `as_turing_model` refuses. More observations than a chain can score is an
# error, because their head would never enter the likelihood. Fewer is not: the
# data is right-aligned, so a series that starts later is scored at the end.
# The one exception is a shortfall of exactly the chain's lead-in, which is the
# pre-0.2.0 meaning of `n`.
#
# This walks the chain rather than building a `DataRequirements`: an
# `IDProblem`'s model body reassembles its `IDModel` on every evaluation, so
# the check runs per evaluation and must not allocate.
function _check_observation_count(model, y_t, n::ModelShape)
    chain = _observation_chain(model)
    lead_in = observation_lead_in(chain)
    streams = _stream_models(chain)
    n_time = _n_time(n)
    deepest = _max_lead_in(lead_in)
    streams === nothing &&
        return _check_stream(:y_t, chain, lead_in, deepest, y_t, n_time)
    for k in keys(streams)
        _check_stream(
            k, streams[k], _stream(lead_in, k), deepest, _stream(y_t, k),
            n_time
        )
    end
    return nothing
end

function _check_stream(
        name::Symbol, model, lead_in::Int, deepest::Int, y_t, n::Int
    )
    supplied = _n_supplied(_data_contract(_consumer(model)).shape, y_t)
    supplied === nothing && return nothing
    n_max = n + deepest - lead_in
    supplied > n_max && throw(
        ArgumentError(
            "as_turing_model: stream `$name` is scored against $n_max " *
                "expected values but $supplied observations were supplied, " *
                "so the first $(supplied - n_max) would never be scored. `n` " *
                "is the number of observations."
        )
    )
    _check_exact_length(Val(_alignment(model)), name, supplied, n_max)
    return _check_pre_020(name, supplied, lead_in, n)
end

# An exact-length component cannot right-align at all, so a series of the wrong
# length is a call that will fail rather than one that scores what it has.
_check_exact_length(::Val, name, supplied, n_max) = nothing
function _check_exact_length(::Val{:exact}, name, supplied, n_max)
    supplied == n_max && return nothing
    throw(
        ArgumentError(
            "as_turing_model: stream `$name` takes exactly $n_max values but " *
                "$supplied were supplied. An aggregation or a reporting " *
                "triangle is indexed against the expected series rather than " *
                "right-aligned to it."
        )
    )
end

# A shortfall of exactly the chain's lead-in is the pre-0.2.0 idiom, where `n`
# was the infection series length and the caller added the lead-in by hand.
# Every observation would still be scored, but against a model run over a
# longer series than the caller means, so name it.
function _check_pre_020(name, supplied, lead_in, n)
    (lead_in > 0 && n - supplied == lead_in) || return nothing
    throw(
        ArgumentError(
            "as_turing_model: stream `$name` was given $supplied " *
                "observations for a model asked to cover $n, short by " *
                "exactly this chain's lead-in, which is the pre-0.2.0 " *
                "meaning of `n`. `n` is now the number of observations, not " *
                "the length of the infection series: pass `length(y)`. If the " *
                "stream really does have that many fewer observations, pad " *
                "its head with `missing`."
        )
    )
end
