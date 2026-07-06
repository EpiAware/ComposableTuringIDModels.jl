# Core architecture: the single light supertype for every model component and
# the generic `as_turing_model` constructor.

@doc raw"
The single light supertype for every model component in `ComposableTuringIDModels`.

Unlike the deep abstract hierarchy used by the original `EpiAware` package, the
prototype keeps a **shallow** tree: one root supertype, and directly beneath it a
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

`as_turing_model` is the single generic entry point of the prototype. Every
concrete model struct implements exactly one

```julia
@model function as_turing_model(m::MyModel, args...; kwargs...)
    ...
end
```

method, and components are composed by sampling submodels of one another:

```julia
z ~ to_submodel(as_turing_model(inner_model, n), false)
```

The trailing `false` to `to_submodel` disables automatic variable prefixing so
that parameter names stay flat unless prefixing is explicitly requested.

The fallback method below errors with a clear message when a struct does not yet
implement `as_turing_model`, which keeps the public surface honest while the
prototype grows.

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
    hint = if model isa AbstractLatentModel
        " expected the latent interface `as_turing_model(m, n)`"
    elseif model isa AbstractInfectionModel
        " expected the infection interface `as_turing_model(m, n)`"
    elseif model isa AbstractObservationModel
        " expected the observation interface `as_turing_model(m, y_t, Y_t)`"
    else
        ""
    end
    throw(ArgumentError(
        "no `as_turing_model` method is defined for $(typeof(model)) with " *
        "$(length(args)) positional argument(s);$hint. Each model struct must " *
        "implement `@model function as_turing_model(m::T, ...)`"))
end
