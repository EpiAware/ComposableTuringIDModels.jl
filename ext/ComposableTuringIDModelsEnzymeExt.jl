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

end # module ComposableTuringIDModelsEnzymeExt
