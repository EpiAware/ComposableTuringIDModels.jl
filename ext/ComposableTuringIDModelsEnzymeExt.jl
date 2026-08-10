# Enzyme AD-safety machinery: `deepcopy` of missing-bearing conditioning data
# is marked inactive, and forward-mode `dot` gets a rule BLAS does not provide.
# Both are generic DynamicPPL/Enzyme fixes, so they belong in
# `EpiAwareADTools.jl`.
module ComposableTuringIDModelsEnzymeExt

using Enzyme: Enzyme
using Enzyme.EnzymeRules: EnzymeRules
using LinearAlgebra: dot

# `T <: Real` is load-bearing (do not widen): `Union{Missing, Any}` collapses to
# `Any`, so an unbounded signature would also match `deepcopy(::Vector{Any})` —
# a differentiable boxed shape — and mark it inactive, silently dropping a real
# tangent.
function EnzymeRules.inactive(::typeof(deepcopy),
        ::AbstractArray{Union{Missing, T}}) where {T <: Real}
    return nothing
end

# Forward-mode rule for `dot`, which BLAS leaves without a derivative.
# `∂/∂x = y`, `∂/∂y = x`, exact and side-effect free. Bounded to
# `AbstractVector{<:Real}` so unrelated `dot` overloads are not intercepted.

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
