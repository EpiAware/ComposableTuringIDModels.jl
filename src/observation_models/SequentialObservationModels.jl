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

# Normalise the many ways a stream can be written into a `SequentialStream`
# (a `transform_chain`, `error_leaf` split) so the recorder can be placed at the
# error leaf and the leaf's input (the post-transform expected series) threaded on.
_to_stream(s::SequentialStream) = s
_to_stream(p::Pair) = SequentialStream(first(p), last(p))
# A bare error leaf is itself the leaf, with an identity transform chain.
_to_stream(m::AbstractObservationErrorModel) = SequentialStream(identity, m)
# A bare *wrapped* observation model (e.g. `LatentDelay(PoissonError(), pmf)`) is
# peeled into its transform chain + error leaf via `_peel_leaf`, so the recorder
# lands at the leaf rather than wrapping the whole stream — the expected series
# that threads is then the leaf input (post-delay/ascertainment), not the raw
# stream input.
_to_stream(m::AbstractObservationModel) = _peel_leaf(m)

@doc raw"
Peel a wrapped observation model into a [`SequentialStream`](@ref): a
`(transform_chain, error_leaf)` split with the error leaf at the bottom and the
delay/ascertainment/transform wrappers as the chain.

This lets a stream be written as a plain nested observation model (e.g.
`LatentDelay(Ascertainment(PoissonError(), eff), pmf)`) yet still have the
internal recorder placed at its error leaf, so the threaded expected series is the
leaf input (the post-transform per-time mean) rather than the raw stream input.
Dissection is restricted to the package's **return-passing** observation modifiers
(`LatentDelay`, `Ascertainment`, `TransformObservationModel`) — those whose
`as_turing_model` re-returns their inner submodel's value, so the leaf recorder's
expected series propagates back up. A modifier that re-consumes its inner's return
(e.g. `Aggregate`) cannot thread `.expected` and is rejected with a clear error;
any other unrecognised wrapper also raises, so the explicit
`transform_chain => error_leaf` form is used instead of guessing.
"
_peel_leaf(m::AbstractObservationErrorModel) = SequentialStream(identity, m)
function _peel_leaf(m::AbstractObservationModel)
    rebuild = _rewrap(m)
    inner = _peel_leaf(m.model)
    return SequentialStream(leaf -> rebuild(inner.transform_chain(leaf)),
        inner.error_leaf)
end

# Functional re-wrappers for the standard observation modifiers: each returns a
# closure that rebuilds the wrapper around a fresh inner model, preserving the
# wrapper's own configuration (delay PMF, ascertainment latent + transform, ...).
#
# Only **return-passing** modifiers may be re-wrapped here: the recorder spliced at
# the leaf returns a `(; obs, expected)` NamedTuple, and that return value must
# propagate back up the wrapper chain unchanged so the stack can read `.expected`.
# `LatentDelay`, `Ascertainment`, and `TransformObservationModel` all end their
# `as_turing_model` with `y_t ~ to_submodel(inner...); return y_t`, i.e. they
# re-RETURN their inner submodel's value, so the recorder's NamedTuple flows
# through. A modifier that re-CONSUMES its inner's return as a plain value (e.g.
# `Aggregate`, which scatters `pred_obs` into a length-`n` vector and returns THAT)
# cannot be leaf-spliced — see the `_returns_inner_value` guard below.
_rewrap(m::LatentDelay) = inner -> LatentDelay(inner, reverse(m.rev_pmf))
function _rewrap(m::Ascertainment)
    # `m.latent_model` is already the (possibly prefix-wrapped) latent model the
    # original was built with; pass `latent_prefix = ""` so it is reused verbatim
    # rather than wrapped a second time.
    return inner -> Ascertainment(inner, m.latent_model, m.transform, "")
end
function _rewrap(m::TransformObservationModel)
    inner -> TransformObservationModel(inner, m.transform)
end

# Whether a wrapper's `as_turing_model` re-RETURNS its inner submodel's value (so a
# leaf recorder's `(; obs, expected)` propagates up) rather than re-consuming it.
# Only the return-passing modifiers are leaf-spliceable; the default is `false` so
# an unaudited modifier is rejected rather than silently mis-threaded.
_returns_inner_value(::AbstractObservationModel) = false
_returns_inner_value(::LatentDelay) = true
_returns_inner_value(::Ascertainment) = true
_returns_inner_value(::TransformObservationModel) = true

# `Aggregate` re-consumes its inner's return (it scatters the predictions into a
# new length-`n` vector via `_return_aggregate` and returns that), so `.expected`
# cannot propagate through it — splicing the recorder underneath would also crash
# (`_return_aggregate` calls `zero` on the recorder's NamedTuple). Reject it (and
# any other return-consuming modifier) with a clear, actionable error.
#
# Future enhancement: a side-channel (e.g. threading the expected series through a
# tracked `:=` quantity or a dedicated return field on every modifier) could let
# `.expected` flow through return-consuming wrappers. That is option B territory
# (a broader return-contract change) and is deliberately out of scope here; keep
# aggregation outside the cascade for now.
function _rewrap(m::AbstractObservationModel)
    if !_returns_inner_value(m)
        return error("`$(nameof(typeof(m)))` cannot be used as a " *
                     "`SequentialObservationModels` stream transform: it re-consumes " *
                     "its inner observation model's return value rather than " *
                     "re-returning it, so the leaf recorder's expected series cannot " *
                     "propagate up the chain to thread into the next stream. Keep " *
                     "return-consuming modifiers (e.g. `Aggregate`) OUTSIDE the " *
                     "sequential cascade — for example apply aggregation to the final " *
                     "stream's output, or to a single non-sequential observation model.")
    end
    return error("Cannot peel `$(typeof(m))` into a (transform-chain, error-leaf) " *
                 "stream for `SequentialObservationModels`. Write this stream " *
                 "explicitly as `transform_chain => error_leaf` (e.g. " *
                 "`(inner -> $(nameof(typeof(m)))(inner, ...)) => PoissonError()`).")
end

# Build a stream with the internal recorder spliced at its error leaf, so
# `as_turing_model` returns `(; obs, expected)` with `expected` the post-transform
# (error-leaf input) series the stack threads onward to the next stream.
#
# EVERY stream — including the first — is instrumented this way so the cascade
# threads genuinely from stream 1: stream 1's observations are still scored on the
# incoming `I_t`, but its `.expected` (the per-time mean AFTER its own transform
# chain) is what seeds stream 2. This is what distinguishes a sequential
# `cases → deaths` cascade (deaths arise from the expected cases) from the parallel
# `StackObservationModels`, where every stream forks off the same `I_t`.
function _instrument_stream(s)
    str = _to_stream(s)
    instrumented = str.transform_chain(_SeqExpectedRecorder(str.error_leaf))
    # Validate the assembled transform chain end-to-end. This catches the explicit
    # `transform_chain => error_leaf` form too (the bare-nested form is already
    # checked in `_rewrap`): every modifier layer between the top of the stream and
    # the leaf recorder must re-RETURN its inner value, or the recorder's
    # `(; obs, expected)` cannot reach the stack.
    _validate_threading(instrumented)
    return instrumented
end

# Walk the assembled stream from the outside in, confirming every wrapper layer
# down to the `_SeqExpectedRecorder` is return-passing. A return-consuming layer
# (e.g. `Aggregate`) both breaks `.expected` propagation and crashes on the
# recorder's NamedTuple, so reject it with the same clear error.
_validate_threading(::_SeqExpectedRecorder) = nothing
function _validate_threading(m::AbstractObservationModel)
    if !_returns_inner_value(m)
        error("`$(nameof(typeof(m)))` cannot be used as a " *
              "`SequentialObservationModels` stream transform layer: it is not a " *
              "known return-passing modifier, so the leaf recorder's expected series " *
              "cannot be guaranteed to propagate up the chain to thread into the " *
              "next stream. Modifiers that re-consume their inner model's return " *
              "value rather than re-returning it (e.g. `Aggregate`) break this and " *
              "must stay OUTSIDE the sequential cascade — apply aggregation to the " *
              "final stream's output, or to a single non-sequential observation " *
              "model. Stream transforms must be built from return-passing modifiers " *
              "(`LatentDelay`, `Ascertainment`, `TransformObservationModel`).")
    end
    # Every return-passing modifier stores its inner observation model in `.model`;
    # recurse toward the recorder at the leaf.
    return _validate_threading(m.model)
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

Each stream — including the first — is a `(transform_chain, error_leaf)` split:
`transform_chain` maps an inner observation model to a wrapped one (e.g.
`inner -> LatentDelay(Ascertainment(inner, eff), pmf)`) and `error_leaf` is the
observation-error model at the leaf. The stack splices a small internal recorder
at each leaf, so every stream exposes its expected (post-transform) output, and
threads that `.expected` series forward as the next stream's expected input.

A stream may be written as:

  - a **bare error model** (e.g. `PoissonError()`) — a pure leaf, identity chain;
  - a **bare wrapped observation model** (e.g.
    `LatentDelay(Ascertainment(PoissonError(), eff), pmf)`) — automatically peeled
    into its transform chain + error leaf so the recorder lands at the leaf; or
  - an explicit `transform_chain => error_leaf` pair (or a
    [`SequentialStream`](@ref)).

The first stream's observations are scored on the incoming `I_t`, but its
`.expected` (the series after its own transform chain) is what seeds the second
stream — so the cascade threads genuinely from stream 1, not only between later
streams. Use the parallel [`StackObservationModels`](@ref) instead if every stream
should fork off the same `I_t`.

Stream transforms must be **return-passing** modifiers — `LatentDelay`,
`Ascertainment`, and `TransformObservationModel`, which re-return their inner
submodel's value so the leaf recorder's expected series propagates up to thread
into the next stream. A modifier that re-consumes its inner's return rather than
re-returning it (e.g. [`Aggregate`](@ref)) breaks that propagation and is rejected
with a clear error; keep such modifiers outside the cascade (e.g. aggregate the
final stream's output, or a single non-sequential observation model).

# Arguments

  - `obs_model`: the [`SequentialObservationModels`](@ref) model.
  - `y_t`: a `NamedTuple` of observed series, one per stream, in stream order.
  - `Y_t`: the expected-observation series fed to the first stream (a vector), or
    a per-stream `NamedTuple` (same names/order as `y_t`) whose first entry seeds
    the first stream — the later entries are ignored, since later streams are fed
    by the cascade.

# Examples
```@example SequentialObservationModels
using EpiAwarePrototype, Distributions
# infections → cases (delayed) → deaths (ascertained off the delayed-case mean).
# `FixedIntercept(log(0.1))` gives a 0.1× effect: Ascertainment's default
# transform is multiplicative on the exponential scale, so it applies `exp(log 0.1)`.
obs = SequentialObservationModels((
    cases = LatentDelay(PoissonError(), [0.4, 0.3, 0.2, 0.1]),
    deaths = (inner -> Ascertainment(inner, FixedIntercept(log(0.1)))) => PoissonError()))
mdl = as_turing_model(obs, (cases = missing, deaths = missing), fill(100.0, 12))
rand(mdl)
```

## Fields

  - `models`: the vector of streams (each prefix-wrapped by its name and
    recorder-instrumented at its error leaf so it exposes its expected output).
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
            # EVERY stream — including the first — is reduced to a
            # (transform-chain, error-leaf) form with the internal recorder placed
            # at the leaf, so it exposes its expected (pre-error) output. That is
            # what makes the cascade thread genuinely from stream 1: stream 1's
            # observations are still scored on the incoming `I_t`, but its
            # post-transform expected series seeds stream 2.
            wrapped[i] = PrefixObservationModel(
                _instrument_stream(models[i]), model_names[i])
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

    # The cascade is seeded with the incoming expected series `I_t`. Each stream
    # (starting with the first) is fed the previous stream's expected (pre-error)
    # output and, via its recorder, exposes its own expected output for the next
    # stream. Stream 1's observations are scored on `I_t`; its `.expected` (the
    # series after stream 1's own transform chain, e.g. a delay) is what seeds
    # stream 2 — a genuine `I_t → stream 1 → stream 2 → …` chain.
    expected = Y_t
    for i in eachindex(obs_model.models)
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
    # A per-stream `NamedTuple` `Y_t` is accepted for symmetry with the parallel
    # stack, but only the FIRST stream is seeded externally — later streams are fed
    # by the cascade. Validate the whole key set (order included) so a misnamed or
    # misordered `Y_t` errors loudly rather than silently using the wrong column.
    @assert keys(Y_t)==keys(y_t) "The keys of the expected series `Y_t` must match the observed series `y_t` (same names, same order); later streams are still fed by the cascade, only the first seeds it"
    seed = Y_t[keys(Y_t)[1]]
    obs ~ to_submodel(as_turing_model(obs_model, y_t, seed), false)
    return obs
end
