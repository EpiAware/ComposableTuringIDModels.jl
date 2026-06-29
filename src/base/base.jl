# Core architecture: the single light supertype for every model component and
# the generic `as_turing_model` constructor.

@doc raw"
A single light supertype for every model component in `EpiAwarePrototype`.

Unlike the deep abstract hierarchy used by the original `EpiAware` package, the
prototype collapses every latent, infection, and observation model under one
supertype. The role a struct plays (latent process, infection process,
observation model, or composed model) is documented and tested behaviourally
rather than encoded in the type tree; dispatch happens on the concrete struct
inside [`as_turing_model`](@ref).
"
abstract type AbstractEpiAwareModel end

@doc raw"
Construct a `DynamicPPL.Model` from an `EpiAwarePrototype` model component.

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

  - `model`: an `EpiAwarePrototype` model component (a subtype of
    [`AbstractEpiAwareModel`](@ref)).
  - `args...`: positional arguments forwarded to the component's method, such as
    the series length `n` (latent models) or the expected/observed series
    (infection and observation models).
  - `kwargs...`: keyword arguments forwarded to the component's method.

# Examples

```@example
using EpiAwarePrototype, Distributions
turing_model = as_turing_model(RandomWalk(), 10)
rand(turing_model)
```
"
function as_turing_model(model, args...; kwargs...)
    throw(ArgumentError(
        "no `as_turing_model` method is defined for $(typeof(model)); each " *
        "model struct must implement `@model function as_turing_model(m::T, ...)`"))
end
