# Unified observation composition: split one expected series into several named
# observation streams, in parallel (shared input), in sequence (one stream's
# expected output feeds the next), or across data-driven strata. Supersedes the
# parallel-only `StackObservationModels` and the sequential cascade of #51/#58.

@doc raw"
A strata mapping supplied as `Y_t` data to a [`Split`](@ref): project a
**multi-stratum** expected series onto a set of observation streams through a
(possibly weighted) linear map.

`strata` is an `inf_strata × time` matrix of per-infection-stratum expected
series; `map` is an `obs_strata × inf_strata` weight matrix. Observation stream
`k` receives the expected series ``\sum_j \mathrm{map}[k, j]\,\mathrm{strata}[j,
:]``. The three infection→observation mappings of issue #45 are all special cases
of one `map`:

  - **one-to-one** — `map = I` (each infection stratum is its own stream);
  - **many-to-one** — an aggregation row sums several infection strata into one
    stream (e.g. age-stratified infections → total hospitalisations);
  - **many-to-many** — a general weight matrix (e.g. age × region → a different
    reporting partition).

The number of observation streams is `size(map, 1)`, read from the **data** at
model-build time (like the series length `n`), so the strata count is not fixed on
the model struct.

## Fields

  - `strata`: the `inf_strata × time` expected series (per infection stratum).
  - `map`: the `obs_strata × inf_strata` weight matrix.

# Examples
```@example StrataMap
using EpiAwarePrototype
# Two infection strata over 5 days, aggregated many-to-one into one stream and
# also kept as a per-stratum stream (a 3 × 2 map: rows = obs streams).
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
        @assert size(map, 2)==size(strata, 1) "The strata map has $(size(map, 2)) columns but the expected series has $(size(strata, 1)) infection strata (rows); they must match"
        return new{M, W}(strata, map)
    end
end

# --- per-stream expected input from the incoming `Y_t` ------------------------

# A single shared expected series is broadcast to every stream (the parallel /
# `StackObservationModels` case: cases and deaths off the same infections).
_split_expected(Y_t::AbstractVector, names) = [Y_t for _ in names]

# A per-stream `NamedTuple` supplies each stream's expected series by name.
function _split_expected(Y_t::NamedTuple, names)
    @assert Set(Symbol.(names)) == Set(keys(Y_t)) "The expected series `Y_t` NamedTuple keys $(keys(Y_t)) must match the stream names $(Tuple(Symbol.(names)))"
    return [Y_t[Symbol(nm)] for nm in names]
end

# A multi-stratum matrix (`inf_strata × time`): stream `k` is row `k` (one-to-one).
function _split_expected(Y_t::AbstractMatrix, names)
    @assert size(Y_t, 1)==length(names) "The expected series matrix has $(size(Y_t, 1)) strata (rows) but there are $(length(names)) streams; supply a `StrataMap` for a non-identity mapping"
    return [Y_t[k, :] for k in 1:length(names)]
end

# A `StrataMap`: project infection strata onto observation streams via the map.
function _split_expected(Y_t::StrataMap, names)
    W = Y_t.map
    M = Y_t.strata
    @assert size(W, 1)==length(names) "The strata map has $(size(W, 1)) observation strata (rows) but there are $(length(names)) streams"
    return [vec(sum(W[k, j] .* M[j, :] for j in 1:size(M, 1))) for k in 1:length(names)]
end

# --- source topology ----------------------------------------------------------

# Resolve the per-stream source list from the constructor keywords. Each entry is
# `:root` (fed the incoming expected series) or the name of an EARLIER stream
# (fed that stream's expected output — the sequential / cascade case).
function _resolve_sources(streamkeys, sequential::Bool, sources)
    ks = collect(streamkeys)
    if sources !== nothing
        @assert !sequential "Pass either `sequential = true` or an explicit `sources`, not both"
        out = Symbol[:root for _ in ks]
        for (k, src) in pairs(sources)
            i = findfirst(==(Symbol(k)), Symbol.(ks))
            i === nothing && error("`sources` names unknown stream `$k`")
            j = findfirst(==(Symbol(src)), Symbol.(ks))
            j === nothing && error("stream `$k` sources unknown stream `$src`")
            @assert j<i "stream `$k` must source an EARLIER stream than itself (streams are evaluated in order); `$src` is not before `$k`"
            out[i] = Symbol(src)
        end
        return out
    elseif sequential
        # A chain: the first stream is fed infections, each later stream is fed the
        # previous stream's expected output.
        return Symbol[i == 1 ? :root : Symbol(ks[i - 1]) for i in eachindex(ks)]
    else
        return Symbol[:root for _ in ks]
    end
end

@doc raw"
Split one expected series into several named observation streams — the single
observation-composition construct for **parallel**, **sequential**, and
**data-driven strata** composition. It replaces the parallel-only
`StackObservationModels` and the sequential cascade sketched in issues #51/#58.

Each stream is a full [`AbstractObservationModel`](@ref) (a bare error family or a
delay / ascertainment / truncation pipeline), wrapped in a
[`PrefixObservationModel`](@ref) keyed by its name so the streams' sampled
variables stay distinct. What differs between the composition modes is only where
each stream's **expected input** comes from — resolved through the uniform
`(; y_t, expected)` observation return contract, which exposes every stream's
pre-error expected series so it can be threaded on or split at any point.

## Composition modes

  - **Parallel** (default) — every stream observes the same incoming expected
    series (e.g. cases and deaths each a delayed fraction of the *same*
    infections). This is the `StackObservationModels` use case, subsumed.
  - **Sequential** (`sequential = true`, or an explicit `sources`) — a stream is
    observed *downstream* of another: its expected input is an earlier stream's
    expected output, so a cascade `I_t → stream 1 → stream 2 → …` forms (e.g.
    deaths as a delayed fraction of *reported cases*). The threaded quantity is
    the earlier stream's **expected** (pre-error, post-transform) series — never
    its sampled output.
  - **Data-driven strata** — constructed from a single **template** observation
    model instead of named streams; the number of streams and their names come
    from the **data** at model-build time (one stream per entry of the `y_t`
    NamedTuple), and a [`StrataMap`](@ref) `Y_t` maps infection strata onto
    observation streams (1:1, many:1, many:many — issue #45).

## Splitting at any point

Because `Split` is itself an observation model that consumes an expected series,
it can be placed **high** (directly on infections) or **low** (after a shared
delay / ascertainment): `Split(streams)` splits on infections, while
`LatentDelay(Split(streams), pmf)` applies a shared delay first and splits the
delayed expectation. The uniform `(; y_t, expected)` contract is what makes this
work at any layer.

## Data contract

`y_t` is a `NamedTuple` of observed series, one per stream, keyed by stream name
(pass `missing` for a stream, or `missing` for the whole model, to simulate). The
return value is the uniform `(; y_t, expected)` pair, where each is a `NamedTuple`
of per-stream series.

`Y_t` (the incoming expected series) may be:

  - a single vector — broadcast to every stream (parallel, shared infections);
  - a per-stream `NamedTuple` — each stream's own expected series;
  - an `inf_strata × time` matrix — one stream per row (one-to-one strata);
  - a [`StrataMap`](@ref) — infection strata projected onto streams (many:many).

## Constructors

  - `Split(streams::NamedTuple; sequential = false, sources = nothing)` — explicit
    named streams. `sequential = true` chains them; `sources = (b = :a, …)` wires
    an explicit dependency DAG (each source must be an earlier stream).
  - `Split(template::AbstractObservationModel)` — a data-driven strata split: the
    template is replicated once per `y_t` entry, prefixed by that entry's name.

# Examples
```@example Split
using EpiAwarePrototype, Distributions
# Parallel: cases and deaths, each a delayed fraction of the SAME infections.
parallel = Split((
    cases = LatentDelay(NegativeBinomialError(), [0.4, 0.3, 0.2, 0.1]),
    deaths = LatentDelay(NegativeBinomialError(), [0.1, 0.2, 0.3, 0.4])))
rand(as_turing_model(parallel, (cases = missing, deaths = missing), fill(100.0, 12)))

# Sequential: deaths as a delayed fraction of REPORTED cases (downstream).
sequential = Split((
    cases = LatentDelay(PoissonError(), [0.5, 0.3, 0.2]),
    deaths = LatentDelay(Ascertainment(PoissonError(), FixedIntercept(log(0.1))),
        [0.2, 0.3, 0.5])); sequential = true)
rand(as_turing_model(sequential, (cases = missing, deaths = missing), fill(100.0, 12)))
```

## Fields

  - `streams`: a `NamedTuple` of per-stream observation models, or a single
    template observation model for the data-driven strata mode.
  - `names`: the stream names (a vector of strings), or `nothing` in strata mode
    (the names come from the `y_t` data).
  - `sources`: the per-stream source list (`:root` or an earlier stream name), or
    `nothing` in strata mode (every stream sources the mapped infections).
"
struct Split{S, N, U} <: AbstractObservationModel
    "Per-stream observation models (`NamedTuple`) or a single strata template."
    streams::S
    "The stream names, or `nothing` in the data-driven strata mode."
    names::N
    "The per-stream source list, or `nothing` in the strata mode."
    sources::U
end

function Split(streams::NamedTuple; sequential::Bool = false, sources = nothing)
    @assert !isempty(streams) "A Split needs at least one stream"
    names = collect(string.(keys(streams)))
    src = _resolve_sources(keys(streams), sequential, sources)
    return Split(streams, names, src)
end

# Data-driven strata: a single template replicated per data stream. The stream
# count and names are supplied by the `y_t` data at model-build time.
Split(template::AbstractObservationModel) = Split(template, nothing, nothing)

# Resolve the ordered stream names: fixed for explicit streams, or the `y_t`
# NamedTuple keys for the data-driven strata mode.
function _split_names(m::Split, y_t::NamedTuple)
    m.names === nothing && return collect(string.(keys(y_t)))
    @assert m.names==collect(string.(keys(y_t))) "The stream names $(m.names) must match the keys of the observed series $(keys(y_t)) (in order)"
    return m.names
end

# The per-stream models: the named entries, or the template replicated per stream.
function _split_models(m::Split, names)
    m.streams isa NamedTuple && return [m.streams[Symbol(nm)] for nm in names]
    return [m.streams for _ in names]   # strata template, one copy per stream
end

# The per-stream source index (0 = root, else the earlier stream's position).
function _split_source_index(m::Split, names)
    m.sources === nothing && return zeros(Int, length(names))
    syms = Symbol.(names)
    return [m.sources[i] === :root ? 0 : findfirst(==(m.sources[i]), syms)
            for i in eachindex(names)]
end

@model function as_turing_model(m::Split, y_t::NamedTuple, Y_t)
    names = _split_names(m, y_t)
    models = _split_models(m, names)
    src_idx = _split_source_index(m, names)
    root_expected = _split_expected(Y_t, names)

    ys = Vector{Any}(undef, length(names))
    exps = Vector{Any}(undef, length(names))
    for i in eachindex(names)
        nm = names[i]
        # `:root` streams read the split expected input; a sourced stream reads its
        # upstream stream's exposed expected (pre-error) output — the cascade hop.
        Yin = src_idx[i] == 0 ? root_expected[i] : exps[src_idx[i]]
        res ~ to_submodel(
            prefix(as_turing_model(models[i], y_t[Symbol(nm)], Yin), Symbol(nm)),
            false)
        ys[i] = res.y_t
        exps[i] = res.expected
    end

    keysyms = Tuple(Symbol.(names))
    return (; y_t = NamedTuple{keysyms}(Tuple(ys)),
        expected = NamedTuple{keysyms}(Tuple(exps)))
end

# Simulate the whole split predictively: expand a `missing` into a per-stream
# NamedTuple of `missing`s (only possible for explicit streams, whose names are
# known without data).
@model function as_turing_model(m::Split, y_t::Missing, Y_t)
    @assert m.names!==nothing "A data-driven strata Split needs a NamedTuple `y_t` to name its streams (pass e.g. `(band1 = missing, …)`)"
    yt = NamedTuple{Tuple(Symbol.(m.names))}(ntuple(_ -> missing, length(m.names)))
    out ~ to_submodel(as_turing_model(m, yt, Y_t), false)
    return out
end
