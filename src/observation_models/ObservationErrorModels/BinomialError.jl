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

This follows the same `NamedTuple`-data pattern as a [`Split`](@ref) stream: the
shared [`define_y_t`](@ref) hook unpacks the `y` field that every error model
scores, and `BinomialError` additionally reads the `N` field it needs.

# Examples
```@example BinomialError
using EpiAwarePrototype
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
_binomial_trials(N::Integer, n) = fill(N, n)
function _binomial_trials(N::AbstractVector{<:Integer}, n)
    @assert length(N)==n "The number-of-trials vector `N` (length $(length(N))) must match the expected-observation series length ($n)"
    return N
end

@model function as_turing_model(obs_model::BinomialError, y_t, Y_t)
    @assert y_t isa NamedTuple&&haskey(y_t, :N) "BinomialError needs `y_t` to be a NamedTuple carrying the number of trials, e.g. `(y = successes, N = trials)` (use `y = missing` to simulate)"

    # Read the number of trials from the data, then rebind `y_t` to the observed
    # successes (the same name, so DynamicPPL conditions on it when concrete).
    N_t = _binomial_trials(y_t.N, length(Y_t))
    y_t = define_y_t(obs_model, y_t, Y_t)

    diff_t = length(y_t) - length(Y_t)
    @assert diff_t>=0 "The observation vector must be at least as long as the expected observation vector"

    # `Y_t` is the success probability; clamp away from 0/1 to avoid a degenerate
    # likelihood, mirroring the count families' `Y_t .+ 1e-6` nudge.
    p_t = clamp.(Y_t, 1e-6, 1 - 1e-6)

    for i in eachindex(Y_t)
        y_t[i + diff_t] ~ observation_error(obs_model, p_t[i], N_t[i])
    end
    return (; y_t, expected = Y_t)
end

observation_error(::BinomialError, p_t, N_t) = Binomial(N_t, p_t)
