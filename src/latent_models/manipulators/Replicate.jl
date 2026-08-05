# Independent per-stratum paths: the same path model drawn once per stratum,
# with no cross-stratum relationship at all.

@doc raw"
Draw `n_strata` independent copies of a path model, one per stratum.

Built with `as_turing_model(m, (n_strata, n_time))`, `Replicate` draws `model`
`n_strata` times over the time axis, each draw prefixed by its stratum index
so the variables of stratum `g` never collide with stratum `h`, and stacks the
results into a `n_strata × n_time` matrix. There is no cross-stratum
relationship: every stratum's path is independent of every other's.

As with [`Stratify`](@ref), the stratum count is **not** a field of the
struct: it arrives through the shape argument, so one `Replicate` serves any
panel width.

Two uses:

  - on its own, in a strata slot, e.g. `rt = Replicate(RandomWalk())` gives
    `n_strata` fully independent `R_t` paths with no shared structure at all;
  - in [`Stratify`](@ref)'s `across` slot, e.g.
    `Stratify(RandomWalk(), Replicate(RandomWalk()))` gives per-stratum
    deviations that vary in time rather than a single constant per-stratum
    offset — see [`across_shape`](@ref), which `Replicate` overrides to draw
    the full `(n_strata, n_time)` shape.

The `model` slot takes a raw prior (a bare `Distribution`, or a latent/prior
model), sampled through [`as_turing_submodel`](@ref); a bare `Distribution` is
auto-wrapped in an [`Intercept`](@ref) (a constant path), the same PATH
convention used elsewhere (e.g. [`RandomWalk`](@ref)'s `ϵ_t`).

## Fields

  - `model`: the path model drawn independently for each stratum.

# Examples
```@example Replicate
using ComposableTuringIDModels, Distributions
# 3 fully independent random-walk paths over a 20-step time axis; the strata
# count arrives through the shape, not a field.
rep = Replicate(RandomWalk())
size(as_turing_model(rep, (3, 20))())
```
"
struct Replicate{M <: PriorLike} <: AbstractPriorModel
    "The path model drawn independently for each stratum."
    model::M

    function Replicate(model)
        # A length-`n_time` PATH slot, so a bare `Distribution` is wrapped in
        # an `Intercept` (a constant path per stratum). The wrapping lives in
        # an inner constructor because the default outer one would be more
        # specific than a plain `Replicate(model) = ...` method and so would
        # win, leaving the slot unwrapped.
        wrapped = _path_prior(model)
        return new{typeof(wrapped)}(wrapped)
    end
end

@model function as_turing_model(m::Replicate, n::Dims{2})
    n_strata, n_time = n
    paths = Vector{Any}(undef, n_strata)
    for g in 1:n_strata
        drawn ~ to_submodel(
            prefix(as_turing_model(m.model, n_time), Symbol(:stratum, g)),
            false)
        paths[g] = drawn
    end
    return permutedims(reduce(hcat, paths))
end

# `Replicate` spans both axes itself, so `Stratify`'s `across` slot must draw
# the full shape rather than the bare stratum count the default gives.
across_shape(::Replicate, n::Dims{2}) = n
