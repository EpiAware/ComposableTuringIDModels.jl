# Random walk latent process model. Its accumulation step (`RWStep`) lives in
# `src/steps/`.

@doc raw"
Model the latent process ``Z_t`` as a random walk.

```math
Z_t = Z_0 + \sum_{i=1}^{t} \epsilon_i
```

where ``Z_0`` is drawn from the prior in `init` and the increments
``\epsilon_i`` come from the error model `ϵ_t` (a `HierarchicalNormal` by
default, giving an inferred step standard deviation).

The `init` slot takes a raw prior: pass a bare `Distribution`, or a richer prior
model. It is sampled through [`as_turing_submodel`](@ref).

`ϵ_t` is a length-`n` PATH slot: a process gives time-varying increments, while
a bare `Distribution` is auto-wrapped in an [`Intercept`](@ref), giving a
**constant** increment path (one shared draw broadcast to every step). Use
[`IID`](@ref) for `n` independent increments.

# Examples
```@example RandomWalk
using ComposableTuringIDModels, Distributions
rw = RandomWalk()
mdl = as_turing_model(rw, 10)
rand(mdl)
```
"
struct RandomWalk{D <: PriorLike, E <: PriorLike} <: AbstractLatentModel
    "Prior for the initial value ``Z_0``."
    init::D
    "Error model for the increments."
    ϵ_t::E
end

function RandomWalk(; init = Normal(), ϵ_t = HierarchicalNormal())
    return RandomWalk(init, _path_prior(ϵ_t))
end

@model function as_turing_model(model::RandomWalk, n)
    @assert n>0 "n must be greater than 0"
    rw_init ~ as_turing_submodel(model.init, 1; prefix = true)
    ϵ_t ~ as_turing_submodel(model.ϵ_t, n - 1)
    rw = accumulate_scan(RWStep(), only(rw_init), ϵ_t)
    return rw
end
