# Observation composition: split one expected series into several named
# observation streams. Each stream sees the expected series arriving at the point
# where `Split` sits in the pipeline, so placement alone gives parallel streams
# (split high, on infections) or a cascade (split low, after a shared delay).

@doc raw"
A strata mapping supplied as the expected series to a [`Split`](@ref): project a
**multi-stratum** expected series onto observation streams through a (possibly
weighted) linear map.

`strata` is an `inf_strata × time` matrix of per-infection-stratum expected
series; `map` is an `obs_strata × inf_strata` weight matrix. Stream `k` receives
``\sum_j \mathrm{map}[k, j]\,\mathrm{strata}[j, :]``, so one `map` covers the
one-to-one (`map = I`), many-to-one (an aggregation row), and many-to-many
(a general weight matrix) infection→observation mappings. The stream count is
`size(map, 1)`, read from the data at model-build time.

## Fields

  - `strata`: the `inf_strata × time` expected series (per infection stratum).
  - `map`: the `obs_strata × inf_strata` weight matrix.

# Examples
```@example StrataMap
using EpiAwarePrototype
M = [10.0 10.0 10.0 10.0 10.0; 4.0 4.0 4.0 4.0 4.0]
W = [1.0 0.0; 0.0 1.0; 1.0 1.0]      # stratum 1, stratum 2, and their sum
sm = StrataMap(M, W)
size(sm.map, 1)                      # 3 observation streams
```
"
struct StrataMap{M <: AbstractMatrix, W <: AbstractMatrix}
    "The `inf_strata × time` expected series (per infection stratum)."
    strata::M
    "The `obs_strata × inf_strata` weight matrix."
    map::W

    function StrataMap(strata::M, map::W) where {M <: AbstractMatrix, W <: AbstractMatrix}
        @assert size(map, 2)==size(strata, 1) "The strata map has $(size(map, 2)) columns but the expected series has $(size(strata, 1)) infection strata (rows)"
        return new{M, W}(strata, map)
    end
end

# --- per-stream expected input from the incoming expected series --------------

# A single shared series is broadcast to every stream (the parallel case).
_split_expected(Y_t::AbstractVector, names) = [Y_t for _ in names]

# A per-stream NamedTuple supplies each stream's series by name.
function _split_expected(Y_t::NamedTuple, names)
    @assert Set(Symbol.(names)) == Set(keys(Y_t)) "The expected series keys $(keys(Y_t)) must match the stream names $(Tuple(Symbol.(names)))"
    return [Y_t[Symbol(nm)] for nm in names]
end

# A multi-stratum matrix: stream `k` is row `k` (one-to-one).
function _split_expected(Y_t::AbstractMatrix, names)
    @assert size(Y_t, 1)==length(names) "The expected series has $(size(Y_t, 1)) strata (rows) but there are $(length(names)) streams; supply a `StrataMap` for a non-identity mapping"
    return [Y_t[k, :] for k in 1:length(names)]
end

# A StrataMap: project infection strata onto streams via the weight matrix.
function _split_expected(Y_t::StrataMap, names)
    W, M = Y_t.map, Y_t.strata
    @assert size(W, 1)==length(names) "The strata map has $(size(W, 1)) observation strata (rows) but there are $(length(names)) streams"
    return [vec(sum(W[k, j] .* M[j, :] for j in 1:size(M, 1))) for k in 1:length(names)]
end

@doc raw"
Split one expected series into several named observation streams — the single
observation-composition construct for **parallel**, **cascade**, and
**data-driven strata** composition.

Each stream is a full [`AbstractObservationModel`](@ref) (a bare error family or
a delay / ascertainment / truncation pipeline). `Split` feeds every stream the
expected series arriving at the point where it sits in the pipeline and
automatically prefixes each stream's sampled variables with its name (via
`DynamicPPL.prefix`), so no manual prefix layer is needed. The uniform
`(; y_t, expected)` return contract exposes each stream's pre-error expected
series so streams can thread on one another.

## Composition by placement

`Split` is itself an observation model consuming an expected series, so *where*
it sits chooses the composition:

  - **Parallel** — placed high, on infections: every stream observes the same
    ``I_t`` (cases and deaths each a delayed, ascertained fraction of the *same*
    infections). `EpiAwareModel(inf, Split((cases = …, deaths = …)))`.
  - **Cascade** — placed low, inside a stream's pipeline: the shared upstream
    layers run first and `Split` branches on their expected output, so a later
    stream is observed *downstream* of an earlier one. For deaths as a delayed
    fraction of the *expected reported cases*, share the case delay then split:
    `LatentDelay(Split((cases = leaf, deaths = pipeline)), case_delay)`.
  - **Data-driven strata** — built from a single **template** model: the number
    and names of streams come from the data at model-build time (one per entry
    of the `y_t` NamedTuple), and a [`StrataMap`](@ref) maps infection strata
    onto observation streams.

The threaded quantity is always a stream's **expected** (pre-error) series, never
its realised noisy draw; observing a downstream stream off another's *sampled*
counts is out of scope.

## Data contract

`y_t` is a NamedTuple of observed series keyed by stream name (or `missing` to
simulate). The return value is `(; y_t, expected)`, each a NamedTuple of
per-stream series. When `Split` is nested inside another modifier the incoming
`missing` reaches it as a shared placeholder; the explicit stream names let it
still fan out.

The incoming expected series may be a single vector (broadcast to every stream),
a per-stream NamedTuple, an `inf_strata × time` matrix (one stream per row), or a
[`StrataMap`](@ref).

## Constructors

  - `Split(streams::NamedTuple)` — explicit named streams.
  - `Split(template::AbstractObservationModel)` — a data-driven strata split
    replicating the template once per `y_t` entry.

# Examples
```@example Split
using EpiAwarePrototype, Distributions
# Parallel: cases and deaths, each a delayed fraction of the SAME infections.
parallel = Split((
    cases = LatentDelay(NegativeBinomialError(), [0.4, 0.3, 0.2, 0.1]),
    deaths = LatentDelay(NegativeBinomialError(), [0.1, 0.2, 0.3, 0.4])))
rand(as_turing_model(parallel, (cases = missing, deaths = missing), fill(100.0, 12)))

# Cascade: share the case delay, then split so deaths sit downstream of cases.
cascade = LatentDelay(
    Split((
        cases = PoissonError(),
        deaths = LatentDelay(
            Ascertainment(PoissonError(), FixedIntercept(log(0.1))),
            [0.2, 0.3, 0.5]))),
    [0.5, 0.3, 0.2])
rand(as_turing_model(cascade, (cases = missing, deaths = missing), fill(100.0, 12)))
```

## Fields

  - `streams`: a NamedTuple of per-stream models, or a single strata template.
  - `names`: the stream names, or `nothing` in strata mode (names come from data).
"
struct Split{S, N} <: AbstractObservationModel
    "Per-stream observation models (NamedTuple) or a single strata template."
    streams::S
    "The stream names, or `nothing` in the data-driven strata mode."
    names::N
end

function Split(streams::NamedTuple)
    @assert !isempty(streams) "A Split needs at least one stream"
    return Split(streams, collect(string.(keys(streams))))
end

# Data-driven strata: a single template replicated per data stream.
Split(template::AbstractObservationModel) = Split(template, nothing)

# Ordered stream names: fixed for explicit streams, else the `y_t` keys.
function _split_names(m::Split, y_t)
    m.names === nothing || return m.names
    y_t isa NamedTuple ||
        error("A strata Split needs a NamedTuple `y_t` to name its streams")
    return collect(string.(keys(y_t)))
end

# Per-stream models: the named entries, or the template replicated per stream.
function _split_models(m::Split, names)
    m.streams isa NamedTuple && return [m.streams[Symbol(nm)] for nm in names]
    return [m.streams for _ in names]
end

# Per-stream data: a NamedTuple splits by name; anything else (a `missing` or a
# nested placeholder vector) is shared to every stream.
function _split_y_t(names, y_t::NamedTuple)
    @assert Set(Symbol.(names)) == Set(keys(y_t)) "The stream names $(Tuple(names)) must match the observed-series keys $(keys(y_t))"
    return [y_t[Symbol(nm)] for nm in names]
end
_split_y_t(names, y_t) = [y_t for _ in names]

@model function as_turing_model(m::Split, y_t, Y_t)
    names = _split_names(m, y_t)
    models = _split_models(m, names)
    expected_in = _split_expected(Y_t, names)
    yt = _split_y_t(names, y_t)

    ys = Vector{Any}(undef, length(names))
    exps = Vector{Any}(undef, length(names))
    for i in eachindex(names)
        nm = names[i]
        res ~ to_submodel(
            prefix(as_turing_model(models[i], yt[i], expected_in[i]),
                Symbol(nm)), false)
        ys[i] = res.y_t
        exps[i] = res.expected
    end

    keysyms = Tuple(Symbol.(names))
    return (; y_t = NamedTuple{keysyms}(Tuple(ys)),
        expected = NamedTuple{keysyms}(Tuple(exps)))
end
