# Structural traversal of an assembled observation chain.
# `wrapped_models` is the contract a component implements, `observation_components`
# is the walk built on it, and `rewrap` puts a component back together around new
# wrapped models.

@doc raw"
The observation models a component wraps directly, as a `Tuple`.

This is the seam an observation chain is traversed through. A component that
consumes the expected series itself — any [`AbstractObservationErrorModel`](@ref)
— wraps nothing and returns `()`. A modifier wraps one model and returns it in a
one-element tuple. A [`Split`](@ref) branches rather than nests, so it returns
every stream at once, in stream order.

## The contract for a new modifier

The default reads the fields of the component that are themselves
`AbstractObservationModel`s, in declaration order, so a modifier holding its
inner model in a plainly typed field needs no method of its own — whatever the
field is called ([`ReportTriangle`](@ref) calls it `error_model`, everything else
`model`).

A component that holds what it wraps inside a container (a `NamedTuple`, a
vector) **must** define its own `wrapped_models`, and a matching
[`rewrap`](@ref). The field's declared type does not matter: the walk tests the
value it finds, so a loosely typed field holding a model directly is read like
any other.
Without one the walk reports it as wrapping nothing, and a traversal reads a
shorter chain than the model actually has.

## Arguments

  - `model`: an observation model.

# Examples
```@example wrapped_models
using ComposableTuringIDModels
obs = LatentDelay(PoissonError(), [0.5, 0.3, 0.2])
ComposableTuringIDModels.wrapped_models(obs)
```
"
function wrapped_models end

# The default is the fields that are themselves observation models, in
# declaration order.
# Read through tuple recursion rather than a generator so the result's type is
# inferred.
function wrapped_models(model::AbstractObservationModel)
    n = fieldcount(typeof(model))
    return _observation_fields(ntuple(i -> getfield(model, i), Val(n)))
end

# Deliberately not the same walk as `_component_children` in
# `base/prettyprinting.jl`, which also descends into `Vector` and `Tuple` fields.
# This is the typed `Tuple` seam `rewrap` rebuilds against, so a component
# holding its wrapped models in a container states its own `wrapped_models`
# rather than having one inferred from field shapes.
# A container-held model is therefore printed but not walked.
_observation_fields(::Tuple{}) = ()
function _observation_fields(values::Tuple)
    rest = _observation_fields(Base.tail(values))
    head = first(values)
    return head isa AbstractObservationModel ? (head, rest...) : rest
end

# An error model consumes the expected series itself, so it is the end of a
# chain.
# Stated rather than inferred, because its prior slots are latent components and
# a future error family holding one must not be mistaken for a wrapper.
wrapped_models(::AbstractObservationErrorModel) = ()

# A `Split`'s streams run in parallel off one expected series, so all of them are
# wrapped at once.
# In the data-driven strata mode `streams` is one template, so one wrapped model.
function wrapped_models(model::Split)
    model.streams isa NamedTuple || return (model.streams,)
    return values(model.streams)
end

@doc raw"
Every component of an observation chain, outermost first.

Returns a `Vector` holding `model` itself, then the models it wraps, depth first
and in order — so a nesting of modifiers comes back in the order the expected
series passes through them, and a [`Split`](@ref) is followed by each of its
streams in turn.

This is what replaces reaching into a chain by field path. Locating a component
is a `filter` over the walk:

```julia
delays = filter(x -> x isa LatentDelay, observation_components(obs))
```

and a quantity accumulated over a chain is a `sum` over it. A quantity that does
*not* simply accumulate reads the branching structure through
[`wrapped_models`](@ref) instead: a [`Split`](@ref)'s streams run in parallel off
one expected series, so what each consumes is a property of the stream rather
than of the chain as a whole.

## Arguments

  - `model`: an observation model, or an [`IDModel`](@ref) (whose observation
    chain is walked).

# Examples
```@example observation_components
using ComposableTuringIDModels, Distributions
obs = LatentDelay(Ascertainment(PoissonError(), FixedIntercept(log(0.1))),
    [0.5, 0.3, 0.2])
ComposableTuringIDModels.observation_components(obs)
```
"
function observation_components(model::AbstractObservationModel)
    return _push_components!(AbstractObservationModel[], model)
end

# A composed model's chain is its observation model's.
function observation_components(model::IDModel)
    return observation_components(model.observation_model)
end

function _push_components!(out::Vector, model::AbstractObservationModel)
    push!(out, model)
    for wrapped in wrapped_models(model)
        _push_components!(out, wrapped)
    end
    return out
end

@doc raw"
Rebuild a component around new wrapped models, returning a new component.

The inverse of [`wrapped_models`](@ref): `rewrap(model, wrapped_models(model))`
reproduces `model`, and passing different models in puts the same wrapper back
around them. Everything else the component holds — a delay PMF, an ascertainment
prior, a [`Split`](@ref)'s stream names and weight map — is carried across
unchanged.

Updating a whole chain is this plus [`wrapped_models`](@ref), via
[`swap`](@ref).

## The contract for a new modifier

The default replaces the component's observation-model fields, in declaration
order, with `models` and rebuilds through `ConstructionBase.constructorof`. A
component that holds what it wraps somewhere the field walk cannot see it (a
[`Split`](@ref)'s `NamedTuple` of streams) must define its own `rewrap`,
alongside its [`wrapped_models`](@ref).

A component whose *public* constructor transforms or derives a field, so that it
cannot accept its own stored fields back, needs no `rewrap` method: it defines
`ConstructionBase.constructorof` instead, which fixes `Accessors` on the type at
the same time. [`LatentDelay`](@ref), [`Ascertainment`](@ref) and
[`Aggregate`](@ref) all do.

## Arguments

  - `model`: the component to rebuild.
  - `models`: a `Tuple` of replacements, one per entry of
    `wrapped_models(model)`, in the same order.

# Examples
```@example rewrap
using ComposableTuringIDModels
obs = LatentDelay(PoissonError(), [0.5, 0.3, 0.2])
ComposableTuringIDModels.rewrap(obs, (NegativeBinomialError(),))
```
"
function rewrap end

function rewrap(model::AbstractObservationModel, models::Tuple)
    _check_rewrap_arity(model, models)
    isempty(models) && return model
    n = fieldcount(typeof(model))
    fields = _substitute_observations(
        ntuple(i -> getfield(model, i), Val(n)), models
    )
    return ConstructionBase.constructorof(typeof(model))(fields...)
end

function _check_rewrap_arity(model, models::Tuple)
    expected = length(wrapped_models(model))
    expected == length(models) && return nothing
    throw(
        ArgumentError(
            "rewrap: $(nameof(typeof(model))) wraps $expected observation " *
                "model(s), but $(length(models)) were supplied."
        )
    )
end

_substitute_observations(::Tuple{}, ::Tuple) = ()
function _substitute_observations(values::Tuple, models::Tuple)
    head = first(values)
    head isa AbstractObservationModel || return (
        head, _substitute_observations(Base.tail(values), models)...,
    )
    return (
        first(models),
        _substitute_observations(Base.tail(values), Base.tail(models))...,
    )
end

# A `Split` holds its streams in a `NamedTuple` rather than in fields of its own,
# so the replacements go back under the same names.
# The stream names and the weight map are the split's own state and are carried
# across.
function rewrap(model::Split, models::Tuple)
    _check_rewrap_arity(model, models)
    model.streams isa NamedTuple ||
        return Split(only(models), model.names, model.map)
    return Split(
        NamedTuple{keys(model.streams)}(models), model.names, model.map
    )
end

@doc raw"
Rebuild every component of a given type in an observation chain.

`swap(f, obs)` applies `f` to the whole chain and to every observation model it
wraps, outermost first, with each component rebuilt by [`rewrap`](@ref) around
its own replacements. `f` returns the replacement for the component it is
given, so a type-keyed swap is a `filter`-like one-liner:

```julia
swap(x -> x isa PoissonError ? NegativeBinomialError() : x, obs)
```

Everything a wrapper holds beyond its inner model — a delay PMF, an
ascertainment prior, a [`Split`](@ref)'s stream names and weight map — is the
`rewrap` contract's to carry across, so it survives the swap. The original
chain is left untouched.

To target a *single* component rather than every one of a type, address it by
position through [`wrapped_models`](@ref) and [`rewrap`](@ref) directly:

```julia
# Replace only the first component an error model wraps.
wrapped = ComposableTuringIDModels.wrapped_models(obs)
inner = rewrap(obs, (NegativeBinomialError(), Base.tail(wrapped)...))
```

## Arguments

  - `f`: a callable mapping an observation model to its replacement, applied to
    the chain and to each component it wraps.
  - `model`: an observation model (or an [`IDModel`](@ref), whose
    observation chain is traversed).

# Examples
```@example swap
using ComposableTuringIDModels
obs = LatentDelay(
    Split((cases = PoissonError(), deaths = PoissonError())),
    [0.5, 0.3, 0.2])
ComposableTuringIDModels.swap(
    x -> x isa PoissonError ? NegativeBinomialError() : x, obs)
```
"
function swap(f, model::AbstractObservationModel)
    return f(rewrap(model, map(x -> swap(f, x), wrapped_models(model))))
end

# A composed model's observation chain is its observation model's.
function swap(f, model::IDModel)
    return IDModel(model.infection_model, swap(f, model.observation_model))
end
