# Differenced latent process modifier.

@doc raw"
Model a latent process as a `d`-fold differenced version of an inner process.

If ``\tilde Z_t`` is the inner (undifferenced) latent path supplied via `model`,
then

```math
\Delta^{(d)} Z_t = \tilde Z_t,
```

and ``Z_t`` is recovered by applying `cumsum` `d` times. The `d` initial terms
are inferred from the prior in `init`; `d` equals the length of that prior.

The `init` slot takes a raw prior: pass a vector of `Distribution`s (its length
sets `d`), or a richer prior model. It is sampled through
[`as_turing_submodel`](@ref).

The `model` slot is a length-`n` PATH slot: a bare `Distribution` there is
auto-wrapped in an [`Intercept`](@ref), giving a constant inner path; use
[`IID`](@ref) for `n` independent draws.

The inner process would otherwise draw its own `init` alongside this model's,
so it is namespaced under `diff`. The bare `init` is this model's integration
constants and `diff.damp`, `diff.init` and so on are parameters of the
*differenced* series. Pin one through the nested form,
`fix(mdl, (diff = (damp = 0.1,),))`. Spelling the dotted path as a single
symbol matches nothing, and `fix` ignores a name the model does not have, so
such a call silently pins nothing.

Composing `DiffLatentModel` over an `AR` gives an ARIMA-style latent process.

# Examples
```@example DiffLatentModel
using ComposableTuringIDModels, Distributions
diff = DiffLatentModel(; model = RandomWalk(), init = [Normal(), Normal()])
mdl = as_turing_model(diff, 10)
rand(mdl)
```

The names the composition produces:

```@example DiffLatentModel
using DynamicPPL: VarInfo
keys(VarInfo(mdl))
```
"
struct DiffLatentModel{M <: PriorLike, P <: PriorLike} <: AbstractLatentModel
    "Underlying (undifferenced) latent model."
    model::M
    "Prior for the initial latent variables."
    init::P
    "Number of times differenced."
    d::Int

    function DiffLatentModel(model, init, d::Int)
        @assert d > 0 "d must be greater than 0"
        _assert_prior_length(init, d, "init")
        # `model` is a length-`n` PATH slot: a bare `Distribution` is wrapped in
        # an `Intercept` (a constant inner path), never left as a scalar.
        wrapped = _path_prior(model)
        return new{typeof(wrapped), typeof(init)}(wrapped, init, d)
    end
end

function DiffLatentModel(model, init::Distribution; d::Int)
    return DiffLatentModel(; model = model, init = fill(init, d))
end

function DiffLatentModel(; model, init = [Normal()])
    d = _prior_order(init)
    return DiffLatentModel(model, init, d)
end

@model function as_turing_model(model::DiffLatentModel, n::Int)
    d = model.d
    @assert n > d "n must be longer than d"
    init ~ as_turing_submodel(model.init, d; prefix = true)
    # The inner process draws its own `init`, so namespace it: the bare
    # `init` above is this model's integration constants and `diff.*` are
    # parameters of the differenced series.
    diff ~ as_turing_submodel(model.model, n - d; prefix = true)
    return _combine_diff(init, diff, d)
end

function _combine_diff(init, diff, d)
    combined = vcat(collect(init), collect(diff))
    for _ in 1:d
        combined = cumsum(combined)
    end
    return combined
end
