# The headline composition: a full infections → observations model assembled from
# two components, each itself an `as_turing_model`. The latent (parameter) process
# is owned by the infection model, not a separate top-level slot.

@doc raw"
A composed epidemiological model linking an infection process and an observation
model.

The infection process owns its own latent (parameter) process internally, so
a composed model is just two parts: infections, then observations.

Sampling [`as_turing_model(model, y_t, n)`](@ref as_turing_model) runs the two
stages as submodels:

```math
I_t \;\xrightarrow{\text{infections}}\;
y_t \;\xrightarrow{\text{observations}}\; \text{data}
```

The returned generated quantities are `(; generated_y_t, expected_y_t, I_t, Z_t)`.
`generated_y_t` is the observation model's sampled `y_t` (the observed-or-simulated
series, or a `NamedTuple` of streams for a [`Split`](@ref)); `expected_y_t` is its
pre-error `expected` series (the uniform observation return contract). `Z_t` is the
infection model's internal latent draw (e.g. the (log) ``R_t`` path), kept
accessible as a generated quantity, or `nothing` for infection models with no
exposable latent (e.g. [`ODEProcess`](@ref)). Pass `y_t = missing` to simulate
from the prior, or a data vector to condition.

## Fields

  - `infection_model`: the infection process model generating ``I_t`` (and its
    internal latent ``Z_t``).
  - `observation_model`: the observation model mapping ``I_t`` to ``y_t``.

# Examples
```@example IDModel
using ComposableTuringIDModels, Distributions
model = IDModel(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
    PoissonError())
mdl = as_turing_model(model, missing, 20)
rand(mdl)
```
"
struct IDModel{I <: AbstractInfectionModel, O <: AbstractObservationModel} <:
       AbstractComposableModel
    "Infection process model generating ``I_t`` (and its internal latent ``Z_t``)."
    infection_model::I
    "Observation model mapping ``I_t`` to ``y_t``."
    observation_model::O
end

@doc raw"
Narrow a data argument to a concrete element type before conditioning a model
on it, when nothing in it is actually `missing`.

DynamicPPL's predictive container (e.g.
`as_turing_model(m, missing, n)().generated_y_t`) keeps a `Union{Missing, T}`
element type even once every entry is concrete — the `missing` placeholder
stays part of the runtime type. Conditioning on that container leaves the
receiving model's data argument `hasmissing`
(`DynamicPPL.convert_model_argument`), which fires a defensive `deepcopy` of
the argument on every model evaluation; Enzyme forward has no `deepcopy` rule
for a `Union{Missing, T}`-eltype array, so this breaks gradients under
`AutoEnzyme(; mode = Enzyme.Forward)`. `concrete_observations` narrows such an
argument to its concrete element type, avoiding both the wasted `deepcopy` and
the Enzyme-forward failure. Data that is genuinely missing in places (e.g. a
reporting-gap series, or `Aggregate`'s unreported time points) is returned
unchanged — narrowing away a real `missing` would be wrong, not just
unnecessary. That genuinely-missing path stays differentiable via the
`EnzymeRules.inactive` mark on `deepcopy` in
`ext/ComposableTuringIDModelsEnzymeExt.jl`.

[`as_turing_model(::IDModel, y_t, n)`](@ref as_turing_model) applies this to
its `y_t` argument automatically, so most callers never need it directly. It
is exposed (not exported) for the same situation elsewhere: any data argument
built by simulating from a model's prior and fed straight back in as
conditioning data, including a `NamedTuple` of streams (applied field-wise,
e.g. `BinomialError`'s `(y, N)` or a [`Split`](@ref) observation's per-stream
series).

# Arguments

  - `y`: the data argument to narrow — an `AbstractArray`, a `NamedTuple` of
    them, or a scalar (e.g. `missing`, returned unchanged).
"
concrete_observations(y::AbstractArray) = any(ismissing, y) ? y : identity.(y)
concrete_observations(y::NamedTuple) = map(concrete_observations, y)
concrete_observations(y) = y

# The actual `@model`, kept private and wrapped by `as_turing_model` below.
# `y_t` must be narrowed with `concrete_observations` BEFORE it is stored as a
# `DynamicPPL.Model` argument: narrowing inside this body would be too late,
# since `DynamicPPL.convert_model_argument`'s `hasmissing`/`deepcopy` check runs
# on the raw stored argument before this body is ever called.
@model function _as_turing_model_idmodel(model::IDModel, y_t, n)
    infections ~ as_turing_submodel(model.infection_model, n)
    I_t = infections.I_t
    Z_t = infections.Z_t
    obs ~ as_turing_submodel(model.observation_model, y_t, I_t)
    # Uniform observation contract: sampled series and pre-error expected.
    generated_y_t = obs.y_t
    expected_y_t = obs.expected
    return (; generated_y_t, expected_y_t, I_t, Z_t)
end

function as_turing_model(model::IDModel, y_t, n)
    _as_turing_model_idmodel(model, concrete_observations(y_t), n)
end

@doc raw"
Convenience 2-argument form: read the infection process's shape from the data.

The observation model and the data together fix the shape of the infection
process, so nothing about it needs to be supplied explicitly or stored on the
model. `as_turing_model(model, Y)` is `as_turing_model(model, Y, shape)` with
`shape` resolved via [`infection_strata`](@ref): the number of infection
strata the observation model consumes given the data's row count, paired with
the data's time length.

Three age strata observed as one hospitalisation stream is
`Split(NegativeBinomialError(), [1.0 1.0 1.0])`; a `1 x T` data matrix then
builds a 3-stratum infection process.

# Examples
```@example IDModel_shape
using ComposableTuringIDModels, Distributions
model = IDModel(
    DirectInfections(;
        Z = Stratify(RandomWalk(), Hierarchy(; across = IID(Normal(0, 0.5)))),
        initialisation = Normal(log(50), 0.2)),
    Split(PoissonError(), [1.0 1.0 1.0]))
# One observation stream, 12 time steps: the three infection strata come
# from the weight matrix, not from the data.
Ymiss = Matrix{Union{Missing, Float64}}(missing, 1, 12)
sim = as_turing_model(model, Ymiss)()
size(sim.I_t)
```
"
function as_turing_model(model::IDModel, Y::AbstractMatrix)
    return as_turing_model(model, Y, _shape(model, Y))
end

function _shape(model::IDModel, Y::AbstractMatrix)
    (infection_strata(model.observation_model, size(Y, 1)), size(Y, 2))
end

@doc raw"
The number of infection strata an observation model consumes.

The seam an observation model uses to say how many infection strata it
consumes, given the number of observation streams in the data. The default
passes the observation stream count straight through (a one-to-one mapping).
[`Split`](@ref) with a weight `map` overrides it with the map's column count,
so a many-to-one or many-to-many mapping can build the right-shaped infection
process from the data alone.

# Arguments

  - `obs`: the observation model.
  - `n_obs_strata`: the number of observation streams in the data.

# Examples
```@example infection_strata
using ComposableTuringIDModels
ComposableTuringIDModels.infection_strata(
    Split(PoissonError(), [1.0 1.0 1.0]), 1)
```
"
infection_strata(obs, n_obs_strata::Int) = n_obs_strata
function infection_strata(s::Split, n_obs_strata::Int)
    s.map === nothing ? n_obs_strata : size(s.map, 2)
end
