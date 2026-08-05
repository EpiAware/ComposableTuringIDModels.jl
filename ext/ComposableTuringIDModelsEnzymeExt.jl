# Enzyme AD-safety machinery (opt-in package extension).
#
# DynamicPPL defensively `deepcopy`s a model argument whenever its type
# `hasmissing` (`DynamicPPL.convert_model_argument`, `src/compiler.jl`) — a
# guard against a tilde-statement like `x[1] ~ Normal()` mutating the caller's
# array in place. That guard runs on EVERY model evaluation, at EVERY level a
# `Union{Missing, T}`-eltype argument is threaded through (a top-level
# `IDModel` call and any submodel it composes), independently of whether the
# argument is fully concrete
# ([`ComposableTuringIDModels.concrete_observations`](@ref) handles that case
# by narrowing the eltype before conditioning) or genuinely partially
# `missing` (a ragged/reporting-gap series, which cannot be narrowed away —
# some entries really are unobserved).
#
# Enzyme forward has no rule for `deepcopy` on a `Vector{Union{Missing, T}}`
# and raises `MethodError: no method matching forward(...,
# ::Const{typeof(deepcopy)}, ...)`. The `deepcopy` here is of observed
# conditioning DATA, not of anything derived from the differentiated
# parameter vector — it carries no tangent in either AD mode — so it is
# exactly the case `EnzymeRules.inactive` exists for:
# Enzyme runs the primal `deepcopy` unchanged and treats it as contributing no
# derivative information, covering every activity/batch-width/mode permutation
# uniformly (the same pattern `EpiAwareADTools` uses for its own `primal`
# tape-strip, in `EpiAwareADToolsEnzymeExt`).
#
# This is generic AD-safety machinery — a DynamicPPL/Enzyme interaction, not
# anything specific to this package's models — so it belongs upstream in
# `EpiAwareADTools.jl` alongside its other Enzyme rules. It lives here for now
# so this package's own PR is self-contained; see the PR description for the
# exact patch proposed for `EpiAwareADToolsEnzymeExt.jl`.
module ComposableTuringIDModelsEnzymeExt

using Enzyme: Enzyme
using Enzyme.EnzymeRules: EnzymeRules
using LinearAlgebra: dot

# The element-type bound here is load-bearing, not decorative — DO NOT widen
# it to a bare `AbstractArray{>:Missing}` (no `T` bound). `Union{Missing, Any}`
# collapses to `Any` in Julia's type lattice, so an unbounded signature would
# also match `deepcopy(::Vector{Any})`: a much broader, genuinely differentiable
# shape (e.g. the boxed `Vector{Any}` the time-varying reporting-delay scenario
# builds by splatting a runtime-length generator, see `backend_broken_scenarios`
# in `ADFixtures.jl`). Marking THAT inactive would silently drop a real
# tangent instead of raising. Bounding `T <: Real` rules out `T = Any`
# (`Any` is not `<: Real`), so this only ever matches the concrete
# `Union{Missing, Int64}` / `Union{Missing, Float64}`-style shapes
# `DynamicPPL.convert_model_argument`'s `hasmissing` branch actually produces
# for observation data — never the boxed-`Any` shape.
#
# Within this package `deepcopy`'s argument here is always a model's plain
# (non-tilde) constructor argument — captured once into `Model.args` at
# construction time, before any differentiation — so it is structurally
# constant with respect to the flat parameter vector a gradient call
# differentiates, regardless of its concrete `T`. `inactive` is scoped by
# call signature, not call site, so it will also fire for any OTHER
# `deepcopy(::AbstractArray{Union{Missing, T}}) where {T <: Real}` in a
# session with both `Enzyme` and this package loaded; that is judged safe
# because `T <: Real` already excludes the one collision found above, but a
# pathological caller that hands `deepcopy` an array of this shape which IS
# itself part of the differentiated computation would still be mismarked —
# the same caveat the upstream `primal`/`_gamma_cdf` rules in
# `EpiAwareADToolsEnzymeExt` carry for their own signatures.
function EnzymeRules.inactive(::typeof(deepcopy),
        ::AbstractArray{Union{Missing, T}}) where {T <: Real}
    return nothing
end

# Forward-mode `LinearAlgebra.dot` rule (renewal / delay-convolution path).
#
# `RenewalSteps.jl`, `LDStep.jl`, `TimeVaryingLDStep.jl` and
# `infection_models/utils.jl` convolve a recent-incidence/state window with a
# generation-interval or reporting-delay kernel via `dot`. Unlike the AR/MA
# lag-order vectors fixed by `ARStep`/`MAStep` (fixed at 2-3, and rewritten to
# `mapreduce` directly — see those files), a generation interval or delay PMF
# is user-supplied epi data with no architectural length bound (routinely
# 10-30+ discretised days), so rewriting every call site to `mapreduce` would
# trade BLAS for a generic reduction on a length nothing here controls. A
# custom rule fixes the actual gap instead: under Enzyme forward, `dot`
# raises `EnzymeNoDerivativeError: dot Runtime Activity not yet implemented
# for Forward-Mode BLAS calls` (Enzyme reverse already has a working BLAS
# `dot` rule and is untouched by this — only `EnzymeRules.forward` is
# defined below, deliberately no `augmented_primal`/`reverse`, so the
# working reverse path is never overridden).
#
# The derivative of `dot(x, y)` is exact and has no edge cases (`∂/∂x = y`,
# `∂/∂y = x`), so this rule cannot silently produce a wrong-but-plausible
# gradient the way a mismarked `inactive`/`EnzymeRules` shape rule could —
# it can only be right or obviously broken, which is why a hand-written rule
# is judged safe here where it would not be for a broader shape-based rule.
# Bounded tightly to `AbstractVector{<:Real}` on both arguments so it only
# ever intercepts the real-valued vector `dot` calls this package makes, not
# every possible `dot` overload (complex vectors, matrices, etc.) that a
# session with both `Enzyme` and this package loaded might otherwise call.
#
# This is an upstream Enzyme.jl gap (forward-mode BLAS `dot` does not support
# runtime activity), not a defect in this package. Delete this rule once
# Enzyme ships its own forward-mode BLAS `dot` support for runtime activity.

# Un-batched: both `x` and `y` are `Const` or plain `Duplicated`.
function EnzymeRules.forward(config::EnzymeRules.FwdConfig,
        func::Enzyme.Const{typeof(dot)},
        RT::Type{<:Union{Enzyme.Const, Enzyme.DuplicatedNoNeed,
            Enzyme.Duplicated}},
        x::Union{Enzyme.Const{<:AbstractVector{<:Real}},
            Enzyme.Duplicated{<:AbstractVector{<:Real}}},
        y::Union{Enzyme.Const{<:AbstractVector{<:Real}},
            Enzyme.Duplicated{<:AbstractVector{<:Real}}})
    dres = zero(promote_type(eltype(x.val), eltype(y.val)))
    x isa Enzyme.Const || (dres += dot(x.dval, y.val))
    y isa Enzyme.Const || (dres += dot(x.val, y.dval))
    if EnzymeRules.needs_primal(config) && EnzymeRules.needs_shadow(config)
        return Enzyme.Duplicated(func.val(x.val, y.val), dres)
    elseif EnzymeRules.needs_primal(config)
        return func.val(x.val, y.val)
    elseif EnzymeRules.needs_shadow(config)
        return dres
    else
        return nothing
    end
end

# Batched: at least one of `x`/`y` is `BatchDuplicated`
# (DifferentiationInterface's Enzyme-forward gradient batches all `dim(θ)`
# seed directions in one pass).
function EnzymeRules.forward(config::EnzymeRules.FwdConfig,
        func::Enzyme.Const{typeof(dot)},
        RT::Type{<:Union{Enzyme.BatchDuplicatedNoNeed, Enzyme.BatchDuplicated}},
        x::Union{Enzyme.Const{<:AbstractVector{<:Real}},
            Enzyme.BatchDuplicated{<:AbstractVector{<:Real}}},
        y::Union{Enzyme.Const{<:AbstractVector{<:Real}},
            Enzyme.BatchDuplicated{<:AbstractVector{<:Real}}})
    N = EnzymeRules.width(config)
    T = promote_type(eltype(x.val), eltype(y.val))
    dvals = ntuple(Val(N)) do i
        dres = zero(T)
        x isa Enzyme.Const || (dres += dot(x.dval[i], y.val))
        y isa Enzyme.Const || (dres += dot(x.val, y.dval[i]))
        dres
    end
    if EnzymeRules.needs_primal(config) && EnzymeRules.needs_shadow(config)
        return Enzyme.BatchDuplicated(func.val(x.val, y.val), dvals)
    elseif EnzymeRules.needs_primal(config)
        return func.val(x.val, y.val)
    elseif EnzymeRules.needs_shadow(config)
        return dvals
    else
        return nothing
    end
end

end # module ComposableTuringIDModelsEnzymeExt
