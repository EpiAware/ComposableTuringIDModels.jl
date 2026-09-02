# Shared observation-error-model machinery: the error supertype, the generic
# `as_turing_model` loop, and the `observation_error` /
# `generate_observation_error_priors` interface.

@doc raw"
Internal supertype shared by simple observation-error models (Poisson, negative
binomial).

It exists only so that the generic observation-error `as_turing_model` loop —
which is identical across error families — can be written once and dispatch the
family-specific pieces ([`observation_error`](@ref) and
[`generate_observation_error_priors`](@ref)) on the concrete type. It is the
error sub-role of [`AbstractObservationModel`](@ref); the package keeps no
deeper hierarchy than this.
"
abstract type AbstractObservationErrorModel <: AbstractObservationModel end

# The steps of `Y_t` that are scored. The observed and expected series end at
# the same time point, so they are right-aligned and the head of whichever is
# longer is left out: surplus observations are the chain's lead-in, and surplus
# expected values are unobserved run-in. A `Split` produces the latter whenever
# one branch's chain consumes a shorter lead-in than another's, since a single
# series length has to serve every branch.
_scored_steps(diff_t, Y_t) = max(1, 1 - diff_t):length(Y_t)

@doc raw"
Generate observations from an observation-error model.

Supports missing observations (`y_t === missing`, simulating predictively) and
observed and expected series of different lengths. The two are right-aligned,
so the last `min(length(y_t), length(Y_t))` entries of each are scored against
one another: a longer `y_t` has its head left unscored (the chain's lead-in),
and a longer `Y_t` has its head left unobserved (run-in the data does not
cover). Expected values are nudged by a tiny constant to avoid degenerate error
distributions.

A series with only *some* entries missing is fitted to the entries it has,
the absent ones being marginalised out as described under
[`MissingObservations`](@ref).

The error family supplies [`generate_observation_error_priors`](@ref) (sampled
as a submodel) and [`observation_error`](@ref) (the per-time-point distribution).

Returns the uniform `(; y_t, expected)` tuple: `y_t` is the observed (or
simulated) counts and `expected` is the pre-error series. Exposing `expected`
lets a [`Split`](@ref) thread one stream's expectation into another.
"
@model function as_turing_model(obs_model::AbstractObservationErrorModel, y_t, Y_t)
    priors ~ to_submodel(
        generate_observation_error_priors(obs_model, y_t, Y_t), false
    )

    pad_Y_t = Y_t .+ 1.0e-6
    # A concrete callable rather than a closure: a closure defined in a model
    # body captures boxed locals, costing a dynamic dispatch per scored entry.
    dist = _ErrorDist(obs_model, pad_Y_t, priors)
    y = _scored_series(obs_model, y_t, Y_t)
    if y isa MissingObservations
        diff_t = length(y.value) - length(Y_t)
        y_t, __varinfo__ = _score_missing_observations!!(
            __model__.context, __varinfo__, y, diff_t, Y_t, dist
        )
    else
        # The count series scored by this model (plain vector, or `missing`,
        # replaced by a length-`Y_t` vector of `missing`). Rebinding `y_t`
        # keeps DynamicPPL treating the entries as conditioned observations.
        y_t = y

        diff_t = length(y_t) - length(Y_t)

        # Every entry is scored, including a `missing` one, so the blanks are
        # drawn. That is how a predictive draw is taken from a fully missing
        # series, and how `forecast` fills a horizon appended to a matrix (a
        # matrix is never split into a carrier). A partially-missing vector
        # conditioned through `IDModel` arrives as a `MissingObservations`
        # carrier above instead, which decides per entry.
        for i in _scored_steps(diff_t, Y_t)
            y_t[i + diff_t] ~ dist(i)
        end
    end
    return (; y_t, expected = Y_t)
end

# Score the observed entries of a `MissingObservations` carrier against a
# per-time-point distribution, without tilde-ing against a `Union{Missing,T}`
# value.
#
# An absent entry is missing at random, so it contributes no likelihood term
# and is skipped rather than sampled (imputing it would add a latent HMC
# cannot link for a count family, and an unneeded gradient dimension for a
# continuous one). Predictive values at gaps come
# from replaying the posterior, not from this loop. Skipping needs no
# `Union{Missing,T}` array, which is exactly the array Enzyme's reverse-mode
# type analysis cannot compile inside a model body; driving
# `DynamicPPL.tilde_observe!!` directly off the carrier's concrete `present`
# mask keeps every array the likelihood touches plainly typed.
#
# Returns the scored series (the observed entries as given, the absent ones
# `missing`) and the updated `VarInfo`.
function (d::_ErrorDist)(i)
    return observation_error(
        d.obs_model, d.pad_Y_t[i], map(p -> at(p, i), values(d.priors))...
    )
end

function (d::_TrialDist)(i)
    return observation_error(d.obs_model, d.p_t[i], d.N_t[i + d.n_diff])
end

function _score_missing_observations!!(
        context, varinfo, y::MissingObservations, diff_t, Y_t, dist
    )
    n = length(y.value)
    # `scored` mixes the observed values with `missing` at the gaps, so it is
    # untyped until every entry is in. Narrowing it with the final `identity.()`
    # (the same idiom `concrete_observations` uses) keeps this step off the
    # model's stored arguments — it runs after every tilde call has already
    # accumulated its log-probability, so it is not part of what Enzyme
    # differentiates.
    scored = Vector{Any}(undef, n)
    for i in _scored_steps(diff_t, Y_t)
        idx = i + diff_t
        # A gap reaches no tilde at all, so it costs neither a `VarName` nor a
        # `VarInfo` entry.
        y.present[idx] || continue
        vn = DynamicPPL.VarName{:y_t}(DynamicPPL.Index((idx,), NamedTuple()))
        val = y.value[idx]
        _, varinfo = DynamicPPL.tilde_observe!!(
            context, dist(i), val, vn, y.value, varinfo
        )
        scored[idx] = val
    end
    # Left unassigned above: every marginalised gap, and any entry outside the
    # scored window (reached when `Y_t` is shorter than `y_t`, e.g. a delay
    # convolution, and so never scored). Keep both as given.
    for idx in 1:n
        isassigned(scored, idx) && continue
        scored[idx] = y.present[idx] ? y.value[idx] : missing
    end
    return identity.(scored), varinfo
end

@doc raw"
Unpack the observed count series an observation-error model scores from the data
`y_t`, dispatching on the model type.

The default method covers every count family (Poisson, negative binomial) and the
Gaussian family: it accepts a plain observation vector, a `missing` (replaced by a
length-`Y_t` vector of `missing` for predictive simulation), a
[`MissingObservations`](@ref) carrier (rebuilt into a ragged
`Vector{Union{Missing,T}}`), or a `NamedTuple` carrying the counts in a `y`
field alongside any extra per-time-point data (a model that needs more than
the counts — e.g. [`BinomialError`](@ref), which also needs the number of
trials — reads those extra fields itself). This keeps the simple case
ergonomic (a plain vector just works) while letting a model opt into a richer
`NamedTuple` data contract.

# Arguments

  - `obs_model`: the observation-error model.
  - `y_t`: the observed data — a vector, `missing`, a `MissingObservations`
    carrier, or a `NamedTuple`.
  - `Y_t`: the expected-observation series (used to size a `missing` series).

# Examples
```@example define_y_t
using ComposableTuringIDModels
# A plain vector passes through; a NamedTuple's `y` field is unpacked.
define_y_t(PoissonError(), [1, 2, 3], fill(10.0, 3)),
define_y_t(PoissonError(), (y = [1, 2, 3],), fill(10.0, 3))
```
"
function define_y_t(::AbstractObservationErrorModel, y_t, Y_t)
    # A NamedTuple carries the counts in its `y` field; a plain value is the
    # counts directly. Either way, a `missing` count series becomes a length-`Y_t`
    # vector of `missing` for predictive simulation.
    y = y_t isa NamedTuple ? y_t.y : y_t
    y = _restore_missing(y)
    return ismissing(y) ? Vector{Missing}(missing, length(Y_t)) : y
end

# Rebuild the ragged `Vector{Union{Missing,T}}` a `MissingObservations` carrier
# stands in for, for standalone `define_y_t` calls that fall through to the
# `~`-sugar loop above (the hot loop itself never takes this path — see
# `_score_missing_observations!!`).
_restore_missing(y) = y
_restore_missing(y::MissingObservations) = map((v, p) -> p ? v : missing, y.value, y.present)

# The series an error model scores, detached from the caller's data.
#
# A series reached through a `NamedTuple` field is not a model argument, so
# DynamicPPL neither copies nor promotes it and the `y_t[i] ~ …` sugar would
# write a blank's draw straight back into the caller's array. Narrowing it with
# [`concrete_observations`](@ref) (what an `IDModel` applies to its own `y_t`)
# turns any blank into a [`MissingObservations`](@ref) carrier scored by
# reading only. A series passed as the model's own `y_t` argument already has
# DynamicPPL's copy, and an `IDModel`-narrowed series is already a carrier or
# a concrete vector, so both pass through untouched.
_scored_series(obs_model, y_t, Y_t) = define_y_t(obs_model, y_t, Y_t)
_scored_series(obs_model, y_t::MissingObservations, Y_t) = y_t
function _scored_series(obs_model, y_t::NamedTuple, Y_t)
    y_t.y isa MissingObservations && return y_t.y
    return concrete_observations(define_y_t(obs_model, y_t, Y_t))
end

@doc raw"
Generate the priors required by an observation-error model. Returns a named
tuple consumed by [`observation_error`](@ref). The default is an empty tuple.

# Arguments

  - `obs_model`: the observation-error model whose priors are generated.
  - `y_t`: the observed series (or `missing` when simulating predictively).
  - `Y_t`: the expected-observation series.

# Examples
```@example generate_observation_error_priors
using ComposableTuringIDModels
m = generate_observation_error_priors(NegativeBinomialError(), missing, fill(10.0, 5))
rand(m)
```
"
@model function generate_observation_error_priors(
        obs_model::AbstractObservationErrorModel, y_t, Y_t
    )
    return NamedTuple()
end

@doc raw"
The per-time-point observation-error distribution given an expected value and
the sampled priors. Each error family implements its own method.

# Arguments

  - `obs_model`: the observation-error model.
  - `Y_t`: the expected observation at a single time point.
  - additional positional arguments: any sampled priors produced by
    [`generate_observation_error_priors`](@ref) for the family (e.g. the squared
    cluster factor for [`NegativeBinomialError`](@ref)).

# Examples
```@example observation_error
using ComposableTuringIDModels
observation_error(PoissonError(), 10.0)
```
"
function observation_error end
