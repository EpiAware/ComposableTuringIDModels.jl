# Core architecture: the single light supertype for every model component and
# the generic `as_turing_model` constructor.

@doc raw"
The single light supertype for every model component in `ComposableTuringIDModels`.

Unlike the deep abstract hierarchy used by the original `EpiAware` package, this
package keeps a **shallow** tree: one root supertype, and directly beneath it a
small set of *role* supertypes — [`AbstractLatentModel`](@ref),
[`AbstractInfectionModel`](@ref), [`AbstractObservationModel`](@ref) (and
[`AbstractObservationErrorModel`](@ref) under the last) — that encode the role a
component plays. There is no deeper `AbstractTuring*` tree and there are no
per-concept `generate_*` functions; dispatch happens on the concrete struct
inside the single generic [`as_turing_model`](@ref).

Encoding the role in the type lets the composer and manipulators constrain their
component slots, so passing a wrong-role component (e.g. an observation model
where a latent model is expected) fails at **construction** rather than at
sampling. See [`AbstractLatentModel`](@ref) and its siblings for the interface
each role's `as_turing_model` must satisfy.
"
abstract type AbstractComposableModel end

@doc raw"
Construct a `DynamicPPL.Model` from an `ComposableTuringIDModels` model component.

`as_turing_model` is the single generic entry point of the package. Every
concrete model struct implements exactly one

```julia
@model function as_turing_model(m::MyModel, args...; kwargs...)
    ...
end
```

method, and components are composed by sampling submodels of one another through
the [`as_turing_submodel`](@ref) seam:

```julia
z ~ as_turing_submodel(inner_model, n)
```

`as_turing_submodel` disables automatic variable prefixing by default so that
parameter names stay flat unless prefixing is explicitly requested
(`prefix = true`).

The fallback method below errors with a clear message when a struct does not yet
implement `as_turing_model`, which keeps the public surface honest.

# Arguments

  - `model`: an `ComposableTuringIDModels` model component (a subtype of
    [`AbstractComposableModel`](@ref)).
  - `args...`: positional arguments forwarded to the component's method, such as
    the series length `n` (latent models) or the expected/observed series
    (infection and observation models).
  - `kwargs...`: keyword arguments forwarded to the component's method.

# Examples

```@example
using ComposableTuringIDModels, Distributions
turing_model = as_turing_model(RandomWalk(), 10)
rand(turing_model)
```
"
function as_turing_model(model, args...; kwargs...)
    hint = if model isa AbstractPriorModel
        " expected the prior/latent interface `as_turing_model(m, n)`"
    elseif model isa AbstractInfectionModel
        " expected the infection interface `as_turing_model(m, n)`"
    elseif model isa AbstractObservationModel
        " expected the observation interface `as_turing_model(m, y_t, Y_t)`"
    else
        ""
    end
    throw(
        ArgumentError(
            "no `as_turing_model` method is defined for $(typeof(model)) with " *
                "$(length(args)) positional argument(s);$hint. Each model struct must " *
                "implement `@model function as_turing_model(m::T, ...)`"
        )
    )
end

@doc raw"
A partially-missing observation vector split into a concrete value vector and
a presence mask, so it carries no `Missing` in its type.

`value[i]` is the observed entry at `i` when `present[i]` is `true`, and an
unused placeholder otherwise. Defined here (rather than in `compose.jl`, where
[`concrete_observations`](@ref) builds one) so that it loads before the
observation-error models that score one directly, further down the include
order.

An absent entry is missing at random, so it is marginalised out: neither scored
nor sampled, and costing no parameter. Predictive values at those points come
from replaying the posterior afterwards, not from carrying them through the fit
(see [`forecast`](@ref)).

# Examples
```@example MissingObservations
using ComposableTuringIDModels: MissingObservations
carrier = MissingObservations([1.0, 0.0, 3.0], [true, false, true])
carrier.value[carrier.present]
```
"
struct MissingObservations{V <: AbstractVector, M <: AbstractVector{Bool}}
    value::V
    present::M
end

# Per-time-point error distributions, as concrete callables. A closure defined
# inside a `@model` body captures boxed locals, which costs a dynamic dispatch
# on every entry of a scoring loop; a struct stays inferable.
struct _ErrorDist{M, P, R}
    obs_model::M
    pad_Y_t::P
    priors::R
end

# `n_diff` right-aligns the trials against the expected series, in the same
# way the observations are aligned, so the trials vector may be given at either
# length. See `_trial_dist`.
struct _TrialDist{M, P, N}
    obs_model::M
    p_t::P
    N_t::N
    n_diff::Int
end
