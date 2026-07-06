# Inference-method supertypes and the `IDMethod` optimisation-then-sampler
# combinator.

@doc raw"
Abstract supertype for inference / generative-modelling methods.
"
abstract type AbstractIDMethod end

@doc raw"
Abstract supertype for optimisation-based methods (e.g. variational
initialisation) used as a pre-sampler step.
"
abstract type AbstractIDOptMethod <: AbstractIDMethod end

@doc raw"
Abstract supertype for sampling-based methods (e.g. NUTS).
"
abstract type AbstractIDSamplingMethod <: AbstractIDMethod end

@doc raw"
Combine a sequence of optimisation pre-steps with a sampler.

`apply_method(model, ::IDMethod)` runs each `pre_sampler_steps` entry in turn,
threading the result into the next step and finally into the `sampler` (e.g.
using a [`ManyPathfinder`](@ref) result to initialise a [`NUTSampler`](@ref)).

## Fields

  - `pre_sampler_steps`: optimisation pre-steps (e.g. Pathfinder).
  - `sampler`: the sampler run last (e.g. NUTS).
"
@kwdef struct IDMethod{O <: AbstractIDOptMethod, S <: AbstractIDSamplingMethod} <:
              AbstractIDMethod
    "Optimisation pre-sampler steps."
    pre_sampler_steps::Vector{O}
    "The sampler run last."
    sampler::S
end
