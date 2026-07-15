# Observation composition: chain several named streams in an *ordered* pipeline,
# each observed downstream of the previous one. Where `Split` fans one expected
# series out to every stream in parallel, `Chain` threads each stream's expected
# output into the next stream's expected input, so a later stream arises from the
# earlier stream's full expectation rather than a shared branch point.

@doc raw"
Chain several named observation streams into an *ordered* pipeline, each observed
downstream of the previous one.

Each stream is a full [`AbstractObservationModel`](@ref) (a bare error family or a
delay / ascertainment / truncation pipeline). The first stream observes the
incoming expected series directly (infections, when placed on an [`IDModel`](@ref)).
Every later stream takes the **expected** (pre-error) series produced by the
stream before it as its own expected input, so the whole of an earlier stream's
delay / ascertainment / truncation propagates into the next. Like [`Split`](@ref),
`Chain` prefixes each stream's sampled variables with its name (via
`DynamicPPL.prefix`) and returns the uniform `(; y_t, expected)` contract as
per-stream NamedTuples.

## Chain versus a Split cascade

[`Split`](@ref) can already put a stream *downstream* of a shared upstream layer by
placement: `LatentDelay(Split((cases = leaf, deaths = pipeline)), delay)` branches
every stream off the **same** delayed expectation. That covers the common case, but
each branch sees the *shared* series at the split point — one stream's own delay or
ascertainment never reaches another. `Chain` differs by threading each stream's
**full expected output** into the next, so with

```julia
Chain((cases = Ascertainment(PoissonError(), …), deaths = LatentDelay(…)))
```

`deaths` is observed off the *ascertained* expected cases, not the pre-ascertainment
infections. Reach for `Chain` when a downstream stream should inherit an upstream
stream's own reporting effects; reach for `Split` when several streams are parallel
fractions of one shared expectation.

## Threaded quantity

The threaded series is always a stream's **expected** (pre-error) series, never its
realised noisy draw — matching [`Split`](@ref). Observing a stream off another's
*sampled* counts is out of scope.

## Delay-shortening down the chain

A [`LatentDelay`](@ref) stream shortens its expected series by the delay length. In a
chain that composes naturally: stream `i` sees whatever length stream `i-1`
produced, so successive delays shorten the series step by step. Order the streams so
each stream's delay is no longer than the series reaching it.

## Data contract

`y_t` is a NamedTuple of observed series keyed by stream name (or `missing` to
simulate, shared to every stream). The return value is `(; y_t, expected)`, each a
NamedTuple of per-stream series in chain order. Because that matches the
[`Split`](@ref) contract a `Chain` nests inside a `Split` (and vice versa) unchanged.

# Examples
```@example Chain
using ComposableTuringIDModels, Distributions
# Deaths observed downstream of the ascertained expected cases: the case
# ascertainment propagates into the death stream, which a Split branch could not do.
chain = Chain((
    cases = Ascertainment(PoissonError(), FixedIntercept(log(0.6))),
    deaths = LatentDelay(
        Ascertainment(PoissonError(), FixedIntercept(log(0.1))),
        [0.2, 0.3, 0.5])))
rand(as_turing_model(chain, (cases = missing, deaths = missing), fill(100.0, 12)))
```

## Fields

  - `streams`: an *ordered* NamedTuple of per-stream observation models.
  - `names`: the stream names in chain order.
"
struct Chain{S <: NamedTuple, N} <: AbstractObservationModel
    "Ordered NamedTuple of per-stream observation models."
    streams::S
    "The stream names in chain order."
    names::N
end

function Chain(streams::NamedTuple)
    @assert !isempty(streams) "A Chain needs at least one stream"
    return Chain(streams, collect(string.(keys(streams))))
end

# Pretty-printing: a Chain's children are its ordered streams, one per line, so the
# chain reads as an ordered list of streams rather than an opaque leaf. Mirrors the
# per-stream walk `Split` uses; the shared tree printer handles connectors/nesting.
function _component_children(m::Chain)
    return Tuple{String, AbstractComposableModel}[(string(nm), m.streams[nm])
                                                  for nm in keys(m.streams)]
end

@model function as_turing_model(m::Chain, y_t, Y_t)
    names = m.names
    models = [m.streams[Symbol(nm)] for nm in names]
    # Per-stream data: a NamedTuple splits by name; a shared `missing` fans out.
    yt = _split_y_t(names, y_t)

    ys = Vector{Any}(undef, length(names))
    exps = Vector{Any}(undef, length(names))
    # `input` is the expected series reaching the current stream: the incoming
    # series for stream 1, then each stream's expected output for the next.
    input = Y_t
    for i in eachindex(names)
        nm = names[i]
        res ~ to_submodel(
            prefix(as_turing_model(models[i], yt[i], input), Symbol(nm)), false)
        ys[i] = res.y_t
        exps[i] = res.expected
        input = res.expected
    end

    keysyms = Tuple(Symbol.(names))
    return (; y_t = NamedTuple{keysyms}(Tuple(ys)),
        expected = NamedTuple{keysyms}(Tuple(exps)))
end
