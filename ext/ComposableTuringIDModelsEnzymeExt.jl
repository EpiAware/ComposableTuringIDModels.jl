# Enzyme AD-safety machinery: `deepcopy` of missing-bearing conditioning data
# is marked inactive, and forward-mode `dot` gets a rule BLAS does not provide.
# Both are generic DynamicPPL/Enzyme fixes, so they belong in
# `EpiAwareADTools.jl`.
module ComposableTuringIDModelsEnzymeExt

using Enzyme: Enzyme
using Enzyme.EnzymeRules: EnzymeRules
using LinearAlgebra: dot

# Keep the `T <: Real` bound. `Union{Missing, Any}` collapses to `Any`, so a
# wider signature would also match `deepcopy(::Vector{Any})` and silently zero a
# real tangent.
function EnzymeRules.inactive(::typeof(deepcopy),
        ::AbstractArray{Union{Missing, T}}) where {T <: Real}
    return nothing
end

# Forward-mode rule for `dot`, which BLAS leaves without a derivative.
# `∂/∂x = y`, `∂/∂y = x`, exact and side-effect free. Bounded to
# `AbstractVector{<:Real}` so unrelated `dot` overloads are not intercepted.

# Runtime activity flags an argument that is inactive at run time by aliasing
# its shadow to its primal, so a rule that reads `dval` unconditionally adds the
# primal to the derivative instead of a tangent. Every aliasing argument must
# therefore be treated as `Const`, the same test Enzyme's own `mul!` rules use.
_rt_const(::EnzymeRules.FwdConfig, ::Enzyme.Const) = true
function _rt_const(config::EnzymeRules.FwdConfig, x)
    return EnzymeRules.runtime_activity(config) && x.dval === x.val
end
_rt_const(::EnzymeRules.FwdConfig, ::Enzyme.Const, ::Int) = true
function _rt_const(config::EnzymeRules.FwdConfig, x, i::Int)
    return EnzymeRules.runtime_activity(config) && x.dval[i] === x.val
end

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
    _rt_const(config, x) || (dres += dot(x.dval, y.val))
    _rt_const(config, y) || (dres += dot(x.val, y.dval))
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

# Batched: DifferentiationInterface seeds every direction in one pass.
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
        _rt_const(config, x, i) || (dres += dot(x.dval[i], y.val))
        _rt_const(config, y, i) || (dres += dot(x.val, y.dval[i]))
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
