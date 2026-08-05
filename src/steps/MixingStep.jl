# A coupling operator that is itself drawn, resolved through the step seam the
# renewal already uses for prior-carrying modifiers.

@doc raw"
Abstract supertype for **drawn** coupling operators.

A fixed coupling operator is a plain matrix and needs no type. A mixing model is
for the case where the operator is *inferred*, so it has to be rebuilt from
sampled parameters on every draw. Its interface is one method,

```julia
as_turing_model(m::MyMixing, n)  # ⇒ a DynamicPPL.Model returning the operator
```

returning whatever [`renewal_pressure`](@ref) accepts, which is usually a
`strata × strata` matrix. [`Gravity`](@ref) is the worked example.

Handing one to [`Renewal`](@ref)'s `mixing` keyword builds a
[`MixingStep`](@ref) instead of a plain core, and the renewal resolves it before
the scan through the same seam that resolves a prior-carrying modifier.
"
abstract type AbstractMixingModel <: AbstractComposableModel end

@doc raw"
A renewal core whose coupling operator is drawn before the scan.

This is to a coupling operator what [`ImportedCases`](@ref) is to an importation
rate. A scan step is a deterministic function, so an operator built from sampled
parameters cannot be assembled inside the recursion. `MixingStep` holds the
generation interval and an [`AbstractMixingModel`](@ref), and its
`as_turing_model` seam draws the model's parameters and hands back a
[`ConstantRenewalStep`](@ref) carrying the realised operator.

It is generic over every mixing model, so a new movement model needs only its
own `as_turing_model` returning a matrix.

## Fields

  - `rev_gen_int`: the reversed generation interval.
  - `mixing`: the mixing model drawn before the scan.
"
struct MixingStep{T, M <: AbstractMixingModel} <: AbstractConstantRenewalStep
    "The reversed generation interval."
    rev_gen_int::T
    "The mixing model drawn before the scan."
    mixing::M
end

# A `MixingStep` is not scannable: its operator does not exist until the seam
# has drawn it. Name the step and point at the seam rather than failing with a
# bare `MethodError` from inside the recursion.
function _unresolved_mixing(step)
    return error("$(typeof(step)) has no scan interface. A drawn coupling " *
                 "operator must be resolved before the scan through the " *
                 "step's `as_turing_model` seam, i.e. " *
                 "`as_turing_submodel(step, n)` (what `Renewal` does).")
end

(step::MixingStep)(window, Rt) = _unresolved_mixing(step)

function _renewal_init_state(step::MixingStep, I₀, r_approx, len_gen_int)
    return _unresolved_mixing(step)
end

@doc raw"
Draw the coupling operator ahead of the scan.

Samples the [`AbstractMixingModel`](@ref)'s parameters through
[`as_turing_submodel`](@ref) and returns the [`ConstantRenewalStep`](@ref) the
scan uses, carrying the realised operator. The operator's variables are
namespaced under `core.mixing`, so a fixed and an inferred coupling differ in
the chain and nowhere else.

# Arguments

  - `step`: the [`MixingStep`](@ref) to resolve.
  - `n`: the shape the scan will run over.
"
@model function as_turing_model(step::MixingStep, n)
    mixing ~ as_turing_submodel(step.mixing, n; prefix = true)
    return ConstantRenewalStep(step.rev_gen_int, mixing)
end

# The renewal core implied by a generation interval and a coupling operator: a
# plain core for a fixed operator, a `MixingStep` for a drawn one. Both reverse
# the interval the same way, so the two paths cannot drift.
function _renewal_step(gen_int, mixing)
    return ConstantRenewalStep(_reverse_lags(gen_int), mixing)
end

function _renewal_step(gen_int, mixing::AbstractMixingModel)
    return MixingStep(_reverse_lags(gen_int), mixing)
end
