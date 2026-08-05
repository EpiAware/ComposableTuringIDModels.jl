# Out-of-sample forecasting: fit at length `T`, predict the observations over a
# future horizon `t = T+1 … T+h`, propagating parameter and latent uncertainty.

@doc raw"
Forecast observations over a future horizon from a fitted model.

Given a `model` fit to an observed series `y` of length ``T`` (yielding `chain`),
`forecast` predicts the observations for the next `horizon` time points
``t = T+1, \dots, T+h`` out of sample. Each posterior draw is carried forward:
the fitted parameters and the in-sample latent path are held fixed while the
latent process is extended over the horizon by drawing its future innovations
from the prior, so the returned forecast propagates both parameter and latent
uncertainty.

`y` extends along the TIME axis: a plain vector of length ``T`` (the
single-series case), a `strata x T` matrix (extended by `horizon` further
columns), or a `NamedTuple` of per-stream series (each stream extended in
turn). The strata count, when there is one, is read from the observation model
and `y` exactly as [`IDProblem`](@ref) reads it at build time.

This works because the package's latent processes are non-centred: a
[`RandomWalk`](@ref), [`AR`](@ref) or [`MA`](@ref) accumulates an i.i.d. sequence
of standard innovations, so the future innovations are independent prior draws
that continue — rather than overwrite — the fitted trajectory. `forecast` extends
each draw's innovation stream to the horizon length with fresh prior draws and
then calls [`predict`](https://turinglang.org/) on the model rebuilt at length
``T + h``.

The result is a chain of the same shape as the input containing the predicted
observations `y_t[T+1] … y_t[T+h]` (the in-sample points stay conditioned on the
data and so are not resampled). Pass it to [`generated_observables`](@ref) or
`returned` to recover the extended latent trajectories per draw.

The extension is exact for the package's non-centred processes because their
future innovations are independent of the fitted history. A latent whose *stored*
stream is itself jointly correlated across time (e.g. an exact-GP `MvNormal`)
would instead need its tail drawn conditional on the history; `forecast` detects
that generically and errors rather than returning a mis-calibrated forecast.

# Arguments

  - `model`: the fitted [`IDModel`](@ref) (or an [`IDProblem`](@ref)).
  - `y`: the observed series the model was fit to (length ``T``).
  - `chain`: the posterior samples from fitting `model` to `y`.
  - `horizon`: the number of future time points ``h`` to forecast.

# Keyword arguments

  - `rng`: random number generator for the future innovations and the predictive
    draws.

# Examples
```@example forecast
using ComposableTuringIDModels, Distributions, Turing, Random
Random.seed!(1)
model = IDModel(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal(1.0, 0.5)),
    PoissonError())
y = fill(5, 15)
chain = sample(as_turing_model(model, y, length(y)), Prior(), 40;
    progress = false)
fc = forecast(model, y, chain, 7)
size(fc)
```
"
function forecast(model::IDModel, y, chain, horizon::Integer;
        rng::AbstractRNG = default_rng())
    horizon ≥ 1 ||
        throw(ArgumentError("horizon must be ≥ 1, got $horizon"))
    n_time = _series_time_length(y)
    y_ext = _extend_series(y, horizon)
    shape = _obs_data_shape(model.observation_model, y, n_time + horizon)
    fc_model = as_turing_model(model, y_ext, shape)
    extended = _extend_latent_draws(rng, fc_model, chain)
    return predict(rng, fc_model, extended)
end

# The observed series' length along the time axis, whatever its shape: a plain
# vector's length, a `strata x time` matrix's column count, or (recursively) a
# `NamedTuple` of streams' shared time length (read off its first stream).
_series_time_length(y::AbstractVector) = length(y)
_series_time_length(y::AbstractMatrix) = size(y, 2)
_series_time_length(y::NamedTuple) = _series_time_length(first(y))

# Extend the observed series by `horizon` further TIME points of `missing`,
# preserving its shape: append to a vector, add columns to a `strata x time`
# matrix, or extend every stream of a `NamedTuple` in turn.
_extend_series(y::AbstractVector, horizon) = vcat(y, fill(missing, horizon))
function _extend_series(y::AbstractMatrix, horizon)
    return hcat(y, fill(missing, size(y, 1), horizon))
end
_extend_series(y::NamedTuple, horizon) = map(v -> _extend_series(v, horizon), y)

@doc raw"
Forecast from a fitted [`IDProblem`](@ref); see [`forecast`](@ref) for the model
method. The observation model and infection process are taken from the problem,
and the horizon extends the problem's `tspan`.
"
function forecast(problem::IDProblem, y, chain, horizon::Integer;
        rng::AbstractRNG = default_rng())
    model = IDModel(problem.infection, problem.observation_model)
    return forecast(model, y, chain, horizon; rng = rng)
end

# Extend the latent innovation streams in `chain` to the length the horizon model
# `fc_model` expects, drawing the extra tail entries from the prior so each draw's
# in-sample path is preserved and its future is a genuine prior continuation.
#
# The fitted chain stores each latent innovation stream as a single vector
# parameter (e.g. a `RandomWalk`'s standard-normal increments) whose length is
# tied to the fitting length `T`. Rebuilding the model at `T + h` asks for the
# same parameter at the longer length, so `predict` alone errors on the
# dimension mismatch. Here every resized vector parameter is padded per draw with
# a fresh prior draw's tail; fixed-length parameters (scalars, AR damping) match
# the model length already and are left untouched. The tail is taken from a fresh
# prior draw of the horizon model, so it follows that stream's *actual* prior
# (e.g. `HierarchicalNormal`'s standard-normal increments — the fitted scale and
# any autoregressive/random-walk correlation are re-applied deterministically by
# `predict` using the fitted parameters, not resampled).
#
# This reaches into the FlexiChains storage that Turing's `predict` consumes; if
# FlexiChains gains a public API for length-extending a parameter this should
# move onto it.
#
# Every vector parameter is a candidate, and the forecast model itself decides
# which are actually extended: a stream is only spliced when the forecast model
# draws MORE of it than the chain holds (the `length(full) > fit_len` test
# below). A parameter on an axis the horizon does not grow — a `Hierarchy`'s
# per-stratum effects, say — is the same length in both models, so it is left
# alone without needing to be recognised as such. A length that grows with the
# horizon is what "time-indexed" means here, which is why no heuristic on the
# time-axis length is used or wanted: an innovation stream is typically shorter
# than the series it drives.
function _extend_latent_draws(rng::AbstractRNG, fc_model, chain)
    extended = deepcopy(chain)
    data = extended._data
    resized = Dict{FlexiChains.Parameter, Int}()
    for key in keys(data)
        key isa FlexiChains.Parameter || continue
        sample = data[key][1, 1]
        sample isa AbstractVector && (resized[key] = length(sample))
    end
    isempty(resized) && return extended
    _assert_factorised(rng, fc_model, resized)
    ni, nc = size(chain)
    for j in 1:nc, i in 1:ni

        prior = Dict(vn => val for (vn, val) in pairs(rand(rng, fc_model)))
        for (key, fit_len) in resized
            full = prior[key.name]
            length(full) > fit_len || continue
            data[key][i, j] = vcat(
                data[key][i, j], full[(fit_len + 1):end])
        end
    end
    return extended
end

# Correctness guard for the independent-tail extension above. Splicing an
# independent prior tail onto the fitted head is exact only when the stream
# *factorises* across the forecast boundary — the future entries are independent
# of the fitted head under the prior. Every latent in this package is non-centred
# (its resized stream is a parameter-free i.i.d. innovation sequence), so this
# holds. A latent whose stored vector is itself jointly correlated (e.g. an
# exact-GP `MvNormal`, or any process with dependence in the *stored* stream)
# would instead need its tail drawn conditional on the head; an independent tail
# would be wrong. Detect that case generically — no per-latent code — by checking
# on a batch of prior draws that each resized stream's forecast tail is
# uncorrelated with its fitted head, and refuse rather than silently mis-forecast.
const _FORECAST_INDEP_TOL = 0.5

function _assert_factorised(rng::AbstractRNG, fc_model, resized)
    n_probe = 256
    draws = [Dict(vn => val for (vn, val) in pairs(rand(rng, fc_model)))
             for _ in 1:n_probe]
    for (key, fit_len) in resized
        full_len = length(draws[1][key.name])
        full_len > fit_len || continue
        mat = reduce(vcat, (permutedims(d[key.name]) for d in draws))
        head = vec(Statistics.mean(view(mat, :, 1:fit_len); dims = 2))
        tail = vec(Statistics.mean(
            view(mat, :, (fit_len + 1):full_len); dims = 2))
        adjacent = Statistics.cor(
            view(mat, :, fit_len), view(mat, :, fit_len + 1))
        block = Statistics.cor(head, tail)
        corr = maximum(
            c -> isfinite(c) ? abs(c) : 0.0, (adjacent, block))
        corr < _FORECAST_INDEP_TOL || error(
            "forecast: latent stream `$(key.name)` is correlated across the " *
            "forecast boundary (|corr| ≈ $(round(corr; digits = 2))), so " *
            "extending it with independent prior draws would be incorrect. " *
            "Forecasting a jointly-correlated latent (e.g. an exact-GP) needs " *
            "conditional extension, which is not yet supported.")
    end
    return nothing
end
