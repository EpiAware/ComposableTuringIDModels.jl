# Binomial observation-error model.

@doc raw"
A binomial observation-error model: the observed successes are binomially
distributed about a per-time-point number of trials `N` and a success
probability supplied by the expected series.

```math
y_t \sim \mathrm{Binomial}(N_t, p_t)
```

Unlike the count error families ([`PoissonError`](@ref),
[`NegativeBinomialError`](@ref)) — whose expected series `Y_t` is an expected
**count** — the expected series passed to `BinomialError` is the success
**probability** ``p_t \in [0, 1]`` (e.g. a prevalence, test-positivity, or
ascertainment proportion).

## The number-of-trials `N` comes from the data

A binomial likelihood needs a number of trials **per time point**, `N_t`. `N` is
*known data* (it is not inferred), so — like the observed successes — it is
supplied through the observation data `y_t`, **not** stored on the model. The
`BinomialError` struct carries no data.

The observation data is a `NamedTuple` with a `y` field (the observed successes)
and an `N` field (the number of trials):

```julia
y_t = (y = observed_successes, N = trials)
```

where `N` is a scalar `Integer` (the same trials at every time point) or an
`AbstractVector{<:Integer}` (per-time-point trials). To **simulate**, pass
`y = missing` while still supplying `N`, e.g.
`y_t = (y = missing, N = fill(20, n))`.

A per-time-point `N` is right-aligned against the expected series in the same
way the observed successes are, so it can be given with one entry per
observation or one per expected value. Both put the same trials on the same
time points; only the unused head differs. `N` must be long enough to cover
every scored time point (the overlap of the two series). This matters for a
[`Split`](@ref) branch whose lead-in is shorter than its neighbours': it gets
more expected values than the caller has observations, and its `N` does not
have to be padded to match.

This follows the same `NamedTuple`-data pattern as a [`Split`](@ref) stream: the
shared [`define_y_t`](@ref) hook unpacks the `y` field that every error model
scores, and `BinomialError` additionally reads the `N` field it needs.

# Examples
```@example BinomialError
using ComposableTuringIDModels
be = BinomialError()
# 20 trials per time point; the expected series is a success probability.
mdl = as_turing_model(be, (y = missing, N = fill(20, 10)), fill(0.3, 10))
rand(mdl)
```
"
struct BinomialError <: AbstractObservationErrorModel end

# `BinomialError` scores the `y` field of the data NamedTuple (shared with the
# other error families through the default `define_y_t`).
define_y_t(::BinomialError, y_t, Y_t) = define_y_t(PoissonError(), y_t, Y_t)

# Resolve the number of trials carried in the data to a per-time-point vector.
# A scalar is broadcast across the expected series; a vector is taken as given
# and aligned by `_trial_dist`.
_binomial_trials(N::Integer, n) = fill(N, n)
_binomial_trials(N::AbstractVector{<:Integer}, n) = N

# Build the per-time-point trials distribution, right-aligning the trials
# against the expected series exactly as the observations are aligned. `n_diff`
# shifts the trials so their last entry pairs with the last expected step, which
# makes the two natural shapes — one entry per observation, or one per expected
# value — put the same trials on the same time points, and leaves the head of
# whichever runs longer unused. Trials that do not reach back to the first
# scored step are rejected rather than quietly shifted.
function _trial_dist(obs_model, p_t, N, diff_t)
    n = length(p_t)
    N_t = _binomial_trials(N, n)
    first_step = first(_scored_steps(diff_t, p_t))
    n_diff = length(N_t) - n
    @assert first_step + n_diff >= 1 "The number of trials `N` (length $(length(N_t))) must cover every scored time point: $(n - first_step + 1) are needed to score $(n + diff_t) observations against $n expected values"
    return _TrialDist(obs_model, p_t, N_t, n_diff)
end

@model function as_turing_model(obs_model::BinomialError, y_t, Y_t)
    @assert y_t isa NamedTuple&&haskey(y_t, :N) "BinomialError needs `y_t` to be a NamedTuple carrying the number of trials, e.g. `(y = successes, N = trials)` (use `y = missing` to simulate)"

    # Read the number of trials from the data before `y_t` is rebound to the
    # observed successes below.
    N = y_t.N

    # `Y_t` is the success probability; clamp away from 0/1 to avoid a degenerate
    # likelihood, mirroring the count families' `Y_t .+ 1e-6` nudge.
    p_t = clamp.(Y_t, 1.0e-6, 1 - 1.0e-6)

    y = y_t.y
    if y isa MissingObservations
        diff_t = length(y.value) - length(Y_t)
        dist = _trial_dist(obs_model, p_t, N, diff_t)
        y_t, __varinfo__ = _score_missing_observations!!(
            __model__.context, __varinfo__, y, diff_t, Y_t, dist
        )
    else
        # Rebind `y_t` to the observed successes (the same name, so DynamicPPL
        # conditions on it when concrete).
        y_t = define_y_t(obs_model, y_t, Y_t)

        diff_t = length(y_t) - length(Y_t)
        dist = _trial_dist(obs_model, p_t, N, diff_t)

        for i in _scored_steps(diff_t, Y_t)
            y_t[i + diff_t] ~ dist(i)
        end
    end
    return (; y_t, expected = Y_t)
end

observation_error(::BinomialError, p_t, N_t) = Binomial(N_t, p_t)
