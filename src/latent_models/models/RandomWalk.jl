# Random walk latent process model. Its accumulation step (`RWStep`) lives in
# `src/steps/`.

@doc raw"
Model the latent process ``Z_t`` as a random walk.

```math
Z_t = Z_0 + \sum_{i=1}^{t} \epsilon_i
```

where ``Z_0`` is drawn from `init_prior` and the increments ``\epsilon_i`` come
from the error model `ϵ_t` (a `HierarchicalNormal` by default, giving an
inferred step standard deviation).

# Examples
```@example RandomWalk
using ComposableTuringIDModels, Distributions
rw = RandomWalk()
mdl = as_turing_model(rw, 10)
rand(mdl)
```
"
@kwdef struct RandomWalk{D <: Sampleable, E <: AbstractLatentModel} <:
              AbstractLatentModel
    init_prior::D = Normal()
    ϵ_t::E = HierarchicalNormal()
end

@model function as_turing_model(model::RandomWalk, n)
    @assert n>0 "n must be greater than 0"
    rw_init ~ model.init_prior
    ϵ_t ~ to_submodel(as_turing_model(model.ϵ_t, n - 1), false)
    rw = accumulate_scan(RWStep(), rw_init, ϵ_t)
    return rw
end
