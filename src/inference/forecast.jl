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

The observations are drawn rather than inferred. Nothing about `y_t` is a
parameter of the fitted model, so the horizon cannot be read out of the
posterior; it is *generated* from it, by replaying each draw through the model
at length ``T + h`` and sampling the observation error at the expected values
that replay produces. This is the same step a Stan `generated quantities` block
performs, and the same one [`predict`](https://turinglang.org/) on a
fully-missing series performs in-sample.

The result is a chain of the same shape as the input holding the predicted
observations, indexable per time point as `y_t[T+1] … y_t[T+h]`. The in-sample
points are generated too, so the same chain carries the in-sample posterior
predictive; the fitted latent path is unchanged either way, because `y_t` never
feeds back into it. Pass the chain to [`generated_observables`](@ref) or
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
function forecast(
        model::IDModel, y, chain, horizon::Integer;
        rng::AbstractRNG = default_rng()
    )
    horizon ≥ 1 ||
        throw(ArgumentError("horizon must be ≥ 1, got $horizon"))
    n_time = _series_time_length(y)
    # The forecast model sees NO observations: every time point is drawn from
    # the error model at the expected values each posterior draw replays. The
    # in-sample data would only pin points the fit has already accounted for,
    # and conditioning on it cannot change the horizon, since `y_t` never feeds
    # back into the latent path.
    y_blank = _blank_series(y, n_time + horizon)
    shape = _obs_data_shape(model.observation_model, y, n_time + horizon)
    fc_model = as_turing_model(model, y_blank, shape)
    extended = _extend_latent_draws(rng, fc_model, chain)
    return predict(rng, fc_model, extended)
end

# The observed series' length along the time axis, whatever its shape: a plain
# vector's length, a `strata x time` matrix's column count, or (recursively) a
# `NamedTuple` of streams' shared time length (read off its first stream).
_series_time_length(y::AbstractVector) = length(y)
_series_time_length(y::AbstractMatrix) = size(y, 2)
_series_time_length(y::NamedTuple) = _series_time_length(first(y))

# A wholly unobserved series of `n_time` time points, shaped like `y`: a vector,
# a `strata x time` matrix, or a `NamedTuple` with the same stream names (which
# is what tells a data-driven `Split` which streams to build). Every entry is
# `missing`, so the observation model draws every point predictively — the same
# argument `as_turing_model(model, fill(missing, n), n)` takes to sample the
# posterior predictive in-sample.
_blank_series(::AbstractVector, n_time) = Vector{Missing}(missing, n_time)
function _blank_series(y::AbstractMatrix, n_time)
    return Matrix{Missing}(missing, size(y, 1), n_time)
end
_blank_series(y::NamedTuple, n_time) = map(v -> _blank_series(v, n_time), y)

@doc raw"
Forecast from a fitted [`IDProblem`](@ref); see [`forecast`](@ref) for the model
method. The observation model and infection process are taken from the problem,
and the horizon extends the problem's `tspan`.
"
function forecast(
        problem::IDProblem, y, chain, horizon::Integer;
        rng::AbstractRNG = default_rng()
    )
    model = IDModel(problem.infection, problem.observation_model)
    return forecast(model, y, chain, horizon; rng = rng)
end

# Extend the latent innovation streams in `chain` to the length the horizon
# model `fc_model` expects, drawing the extra tail entries from the prior so
# each draw's in-sample path is preserved and its future is a genuine prior
# continuation.
#
# Each stored innovation stream's length is tied to the fitting length `T`;
# rebuilding the model at `T + h` asks for the same parameter at the longer
# length, so `predict` alone errors on the dimension mismatch. Every vector
# parameter is a candidate and the forecast model decides which are extended:
# a stream is spliced only when the model draws MORE of it than the chain
# holds. Parameters on an axis the horizon does not grow (a `Hierarchy`'s
# per-stratum effects, say) are the same length in both models and are left
# alone. No heuristic on the time-axis length is used or wanted: an innovation
# stream is typically shorter than the series it drives.
#
# The tail is a fresh prior draw of the horizon model, so it follows
# that stream's actual prior; fitted scale and correlation are re-applied
# deterministically by `predict`, not resampled.
#
# This reaches into the FlexiChains storage `predict` consumes; a public
# length-extension API there would replace it.
function _extend_latent_draws(rng::AbstractRNG, fc_model, chain)
    extended = deepcopy(chain)
    data = extended._data
    resized = Dict{FlexiChains.Parameter, Int}()
    # A generated quantity (a `:=` value) is recomputed by `predict`, not
    # extended: it is not a model input, so a prior draw of the forecast model
    # never samples it and there is no tail to draw. Probe one prior draw to
    # tell input streams from generated quantities.
    prior = Dict(vn => val for (vn, val) in pairs(rand(rng, fc_model)))
    for key in keys(data)
        key isa FlexiChains.Parameter || continue
        haskey(prior, key.name) || continue
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
                data[key][i, j], full[(fit_len + 1):end]
            )
        end
    end
    return extended
end

# Correctness guard for the independent-tail extension above. Splicing an
# independent prior tail onto the fitted head is exact only when the stream
# *factorises* across the forecast boundary. Every latent in this package is
# non-centred (its resized stream is a parameter-free i.i.d. innovation
# sequence), so this holds. A stored stream with joint correlation (an
# exact-GP `MvNormal`) would need its tail drawn conditional on the head.
# Detect that case generically by checking on a batch of prior draws that each
# resized stream's tail is uncorrelated with its head, and refuse rather than
# silently mis-forecast.
const _FORECAST_INDEP_TOL = 0.5

function _assert_factorised(rng::AbstractRNG, fc_model, resized)
    n_probe = 256
    draws = [
        Dict(vn => val for (vn, val) in pairs(rand(rng, fc_model)))
            for _ in 1:n_probe
    ]
    for (key, fit_len) in resized
        full_len = length(draws[1][key.name])
        full_len > fit_len || continue
        mat = reduce(vcat, (permutedims(d[key.name]) for d in draws))
        head = vec(Statistics.mean(view(mat, :, 1:fit_len); dims = 2))
        tail = vec(
            Statistics.mean(
                view(mat, :, (fit_len + 1):full_len); dims = 2
            )
        )
        adjacent = Statistics.cor(
            view(mat, :, fit_len), view(mat, :, fit_len + 1)
        )
        block = Statistics.cor(head, tail)
        corr = maximum(
            c -> isfinite(c) ? abs(c) : 0.0, (adjacent, block)
        )
        corr < _FORECAST_INDEP_TOL || error(
            "forecast: latent stream `$(key.name)` is correlated across the " *
                "forecast boundary (|corr| ≈ $(round(corr; digits = 2))), so " *
                "extending it with independent prior draws would be incorrect. " *
                "Forecasting a jointly-correlated latent (e.g. an exact-GP) needs " *
                "conditional extension, which is not yet supported."
        )
    end
    return nothing
end
