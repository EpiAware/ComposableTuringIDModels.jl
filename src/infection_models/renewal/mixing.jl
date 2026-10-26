# Coupling between the strata of a renewal: the fixed forms ComposableRecurrences
# takes directly, and a coupling operator that is itself drawn, resolved
# before the scan.

@doc raw"
Abstract supertype for **drawn** coupling operators.

A fixed coupling operator is a plain matrix and needs no type. A mixing model is
for the case where the operator is *inferred*, so it has to be rebuilt from
sampled parameters on every draw. Its interface is one method,

```julia
as_turing_model(m::MyMixing, n)  # ⇒ a DynamicPPL.Model returning the operator
```

returning a coupling a [`Renewal`](@ref)'s `mixing` keyword accepts, which is
usually a `strata × strata` matrix. [`Gravity`](@ref) is the worked example.

The renewal draws it before the scan, through the same seam that resolves a
prior-carrying [`AbstractRenewalModifier`](@ref).
Its variables are namespaced under `mixing`, so a fixed and an inferred
coupling differ in the chain and nowhere else.
"
abstract type AbstractMixingModel <: AbstractComposableModel end

# A fixed coupling samples nothing and is used as it is; `n` is the shape a
# drawn one is built at.
@model function _resolve_coupling(mixing, n)
    return mixing
end

@model function _resolve_coupling(mixing::AbstractMixingModel, n)
    mixing ~ as_turing_submodel(mixing, n; prefix = true)
    return mixing
end

# The renewal recurrence a generation interval and a coupling make.
# The interval runs lag 1 first, as ComposableRecurrences takes it.
# A `strata × lags` interval is one row per stratum.
# A `strata × strata × lags` coupling carries a generation interval per pair
# of strata (see `pairwise_gen_int`), so it is the kernel and the shared
# interval only sets the growth rate the seeding window decays at.
_renewal_kernel(gen_int::AbstractVector) = gen_int
_renewal_kernel(gen_int::AbstractMatrix) = PerStratum(gen_int)

function _renewal_recurrence(gen_int, coupling, modifiers)
    return Recurrence(
        _renewal_kernel(gen_int); coupling = coupling, modifiers = modifiers
    )
end

function _renewal_recurrence(
        gen_int, coupling::AbstractArray{<:Any, 3}, modifiers
    )
    return Recurrence(Pairwise(coupling); modifiers = modifiers)
end
