# Sequentially-linked observation models: an ordered cascade where each stream's
# expected (pre-error) output feeds the next stream's expected input.

@doc raw"
Internal error-leaf recorder used by [`SequentialObservationModels`](@ref).

A later stream in a sequential cascade must expose its *expected* (pre-error)
output so the stack can thread it into the next stream as that stream's expected
input. The series that should thread is the **error-leaf input** — the per-time
mean of the stream's observations *after* its transform chain (delay,
ascertainment, ...) has been applied — not the sampled output (which carries
discrete noise and is shortened/`missing`-headed by a delay kernel).

`_SeqExpectedRecorder` wraps the stream's error leaf. Its `as_turing_model`
returns `(; obs, expected)`: `expected` is the `Y_t` handed to the leaf (the
threaded quantity) and `obs` is the leaf's sampled output. The stack inserts it at
the leaf — so it sees the post-transform mean — and reads `.expected` off the
submodel return to feed the next stream. It is **internal**: it is never exported
and never alters the public observation return contract (this is option A of
issue #51; the alternative — making every observation model return
`(; y_t, expected)` — is option B and is deliberately not taken).

## Fields

  - `model`: the wrapped error-leaf observation model.
"
struct _SeqExpectedRecorder{M <: AbstractObservationModel} <: AbstractObservationModel
    "The wrapped error-leaf observation model."
    model::M
end

@model function as_turing_model(rec::_SeqExpectedRecorder, y_t, Y_t)
    # `Y_t` is the error-leaf input: the post-transform expected series. Thread
    # this (not the sampled `obs`) onward.
    obs ~ to_submodel(as_turing_model(rec.model, y_t, Y_t), false)
    return (; obs, expected = Y_t)
end

@doc raw"
A later stream in a [`SequentialObservationModels`](@ref) cascade: a
`(transform_chain, error_leaf)` pair.

`transform_chain` is a function mapping an inner observation model to a wrapped
observation model — the delay/ascertainment/transform stages applied to this
stream's expected input (e.g. `inner -> LatentDelay(Ascertainment(inner, eff),
pmf)`). `error_leaf` is the observation-error model applied at the leaf. The stack
splices a small internal recorder between the chain and the leaf so the stream
exposes the *error-leaf input* (its post-transform expected series) for threading
to the next stream.

A `transform_chain => error_leaf` pair is accepted as shorthand, and a bare
[`AbstractObservationModel`](@ref) is treated as a pure error leaf with an
identity transform chain.

## Fields

  - `transform_chain`: maps an inner observation model to the wrapped stream.
  - `error_leaf`: the observation-error model at the leaf.
"
struct SequentialStream{F, L <: AbstractObservationModel} <: AbstractObservationModel
    "Maps an inner observation model to the wrapped stream."
    transform_chain::F
    "The observation-error model at the leaf."
    error_leaf::L
end

# Normalise the many ways a later stream can be written into a `SequentialStream`.
_to_stream(s::SequentialStream) = s
_to_stream(p::Pair) = SequentialStream(first(p), last(p))
# A bare observation model later in the chain is a pure error leaf (identity
# transform chain) — it scores the threaded expected series directly.
_to_stream(m::AbstractObservationModel) = SequentialStream(identity, m)

# Stream 1 observes `I_t` directly, so it stays a full observation model. A
# `SequentialStream`/`Pair` written for slot 1 is still honoured by realising the
# concrete observation model around its leaf (no recorder — the first stream's
# observations are scored on `I_t`; the cascade threads the post-transform mean of
# the LATER streams).
_as_first_stream(m::AbstractObservationModel) = m
function _as_first_stream(s::Union{SequentialStream, Pair})
    str = _to_stream(s)
    return str.transform_chain(str.error_leaf)
end

# Build a later stream with the internal recorder spliced at its error leaf, so
# `as_turing_model` returns `(; obs, expected)` with `expected` the post-transform
# (error-leaf input) series the stack threads onward.
function _instrument_stream(s)
    str = _to_stream(s)
    return str.transform_chain(_SeqExpectedRecorder(str.error_leaf))
end

@doc raw"
Sequentially link several observation models into an ordered cascade, each
applied to a named component of the data.

`SequentialObservationModels` is the *sequential* counterpart to
[`StackObservationModels`](@ref). Both wrap each stream in a
[`PrefixObservationModel`](@ref) keyed by its name and share the same per-stream
`NamedTuple` data contract for `y_t`. The difference is how the streams are fed:
the parallel stack applies every stream to the **same** expected series; the
sequential stack threads each stream's expected (pre-error) output into the next
stream as *its* expected input. Order is therefore significant — the streams form
a chain `I_t → stream 1 → stream 2 → …`.

The threaded quantity is the **expected** (pre-error) series, captured at each
stream's error leaf — never the sampled output. The sampled output cannot be
threaded: a realistic stream (e.g. `LatentDelay(Ascertainment(error))`) leaves
the head of its sampled series `missing` and shortens it by the delay kernel, so
feeding it forward as an expectation is both statistically wrong (discrete noise
propagated as a mean) and mechanically broken (length/`missing`). The error-leaf
input is the post-transform per-time mean, which differs from the stream's outer
input in both scale and length — exactly the series the next stream should see.

## Stream contract

  - **The first stream** is a full [`AbstractObservationModel`](@ref) observing
    the incoming expected series `I_t` (unchanged observation-model semantics).
  - **Each later stream** is a `(transform_chain, error_leaf)` pair (a
    [`SequentialStream`](@ref), or the `transform_chain => error_leaf` shorthand):
    `transform_chain` maps an inner observation model to a wrapped one (e.g.
    `inner -> LatentDelay(Ascertainment(inner, eff), pmf)`), and `error_leaf` is
    the observation-error model applied at the leaf. The stack inserts a small
    internal recorder at the leaf, so the stream exposes its expected output, and
    threads that `.expected` series into the next stream. A bare
    [`AbstractObservationModel`](@ref) given as a later stream is treated as a
    pure error leaf (an identity transform chain).

# Arguments

  - `obs_model`: the [`SequentialObservationModels`](@ref) model.
  - `y_t`: a `NamedTuple` of observed series, one per stream, in stream order.
  - `Y_t`: the expected-observation series fed to the first stream (a vector), or
    a per-stream `NamedTuple` whose first entry seeds the first stream (the later
    entries are ignored — later streams are fed by the cascade).

# Examples
```@example SequentialObservationModels
using EpiAwarePrototype, Distributions
# infections → cases (delayed) → deaths (ascertained off the delayed-case mean)
obs = SequentialObservationModels((
    cases = LatentDelay(PoissonError(), [0.4, 0.3, 0.2, 0.1]),
    deaths = (inner -> Ascertainment(inner, FixedIntercept(0.1))) => PoissonError()))
mdl = as_turing_model(obs, (cases = missing, deaths = missing), fill(100.0, 12))
rand(mdl)
```

## Fields

  - `models`: the vector of streams (each prefix-wrapped by its name); stream 1 is
    a full observation model, later streams are recorder-instrumented streams.
  - `model_names`: the names identifying each stream, in cascade order.
"
struct SequentialObservationModels{
    M <: AbstractVector, N <: AbstractVector{<:AbstractString}} <:
       AbstractObservationModel
    "The vector of streams (each prefix-wrapped by its name)."
    models::M
    "The names identifying each stream, in cascade order."
    model_names::N

    function SequentialObservationModels(
            models::M, model_names::N) where {
            M <: AbstractVector, N <: AbstractVector{<:AbstractString}}
        @assert length(models)==length(model_names) "The number of streams and stream names must be equal"
        @assert length(models)>=1 "A sequential cascade needs at least one stream"
        wrapped = Vector{AbstractObservationModel}(undef, length(models))
        for i in eachindex(models)
            # Stream 1 is observed directly on `I_t`, so it is a full observation
            # model used as-is. Each later stream is reduced to a
            # (transform-chain, error-leaf) form with the internal recorder placed
            # at the leaf, so it exposes its expected output for threading.
            inner = i == 1 ? _as_first_stream(models[i]) :
                    _instrument_stream(models[i])
            wrapped[i] = PrefixObservationModel(inner, model_names[i])
        end
        return new{typeof(wrapped), N}(wrapped, model_names)
    end
end

function SequentialObservationModels(models::NamedTuple)
    model_names = keys(models) .|> string |> collect
    return SequentialObservationModels(collect(values(models)), model_names)
end

@model function as_turing_model(
        obs_model::SequentialObservationModels, y_t::NamedTuple, Y_t::AbstractVector)
    @assert length(obs_model.models)==length(y_t) "The number of streams must match the number of observed series"
    @assert obs_model.model_names==(keys(y_t) .|> string |> collect) "The stream names must match the keys of the observed series (in order)"

    obs = Vector{Any}(undef, length(obs_model.models))

    # Stream 1: a full observation model on the incoming expected series `I_t`.
    name1 = obs_model.model_names[1]
    obs_1 ~ to_submodel(
        as_turing_model(obs_model.models[1], y_t[Symbol(name1)], Y_t), false)
    obs[1] = obs_1

    # The first stream is a plain observation model returning its sampled `y_t`,
    # so it does not itself expose an expected series. The cascade is seeded with
    # the same expected input the first stream saw; the transform chain inside each
    # later stream (delay, ascertainment, ...) re-shapes it from there.
    expected = Y_t

    # Each later stream is fed the previous stream's expected (pre-error) output
    # and, via its recorder, exposes its own expected output for the next stream.
    for i in 2:length(obs_model.models)
        name = obs_model.model_names[i]
        stream_i ~ to_submodel(
            as_turing_model(obs_model.models[i], y_t[Symbol(name)], expected),
            false)
        obs[i] = stream_i.obs
        expected = stream_i.expected
    end

    return obs
end

@model function as_turing_model(
        obs_model::SequentialObservationModels, y_t::NamedTuple, Y_t::NamedTuple)
    # Only the first stream is seeded externally; later streams are fed by the
    # cascade, so just the leading expected series is used.
    @assert keys(y_t)[1]==keys(Y_t)[1] "The first key of the observed and expected series must match"
    seed = Y_t[keys(Y_t)[1]]
    obs ~ to_submodel(as_turing_model(obs_model, y_t, seed), false)
    return obs
end
