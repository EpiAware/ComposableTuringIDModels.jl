# Inference-method supertypes.

@doc raw"
Abstract supertype for inference / generative-modelling methods.
"
abstract type AbstractIDMethod end

@doc raw"
Abstract supertype for sampling-based methods (e.g. NUTS).
"
abstract type AbstractIDSamplingMethod <: AbstractIDMethod end
