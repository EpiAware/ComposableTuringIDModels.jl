# Record-the-expected-infections modifier, and the seam it reads through.

@doc raw"
The infection model `model`, rebuilt to return its noise-free expectation as
`exp_I_t` alongside its usual output.

The seam [`RecordExpectedInfections`](@ref) reads through. An infection model
implements a method here when it has a noise layer, so that the expectation is
the value entering that layer rather than the drawn series. The default method
returns the committed `I_t`, which is the honest answer for a model that draws
no infections. A [`DirectInfections`](@ref) path, or a [`Renewal`](@ref) with
no [`is_noise`](@ref) modifier, is its own expectation.

# Arguments

  - `model`: the infection model.
  - `n`: the shape to generate, as [`as_turing_model`](@ref) takes it.

# Examples
```@example with_expected_infections
using ComposableTuringIDModels, Distributions
inf = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
out = ComposableTuringIDModels.with_expected_infections(inf, 5)()
out.exp_I_t == out.I_t
```
"
@model function with_expected_infections(
        model::AbstractInfectionModel, n::ModelShape
    )
    inner ~ as_turing_submodel(model, n)
    return (; inner..., exp_I_t = inner.I_t)
end

@doc raw"
Record an infection model's noise-free expectation as a tracked generated
quantity (`exp_I_t`).

The wrapped model runs unchanged, keeping its own variable names, and the
expectation is tracked with `:=` so a chain recovers it as `chain[:exp_I_t]`.
Where the wrapped model draws infections, as a [`StochasticRenewal`](@ref) does
and as a [`Renewal`](@ref) carrying an [`InfectionNoise`](@ref) modifier does,
`exp_I_t` is the incidence entering the draw, so it is
post-[`SusceptibleDepletion`](@ref) and pre-noise. Where the model draws none it
is the committed series itself.

A model used unwrapped computes no expectation at all.

# Arguments

  - `model`: the infection model whose expectation is recorded.

# Examples
```@example RecordExpectedInfections
using ComposableTuringIDModels, Distributions, Random
rec = RecordExpectedInfections(
    Renewal(
        [0.2, 0.3, 0.5], InfectionNoise(); rt = FixedIntercept(0.0),
        initialisation = Normal(log(100.0), 0.0)
    )
)
out = as_turing_model(rec, 8)(Xoshiro(1))
out.exp_I_t
```

## Fields

  - `model`: the infection model whose expected infections are recorded.
"
struct RecordExpectedInfections{M <: AbstractInfectionModel} <:
    AbstractInfectionModel
    "The infection model whose expected infections are recorded."
    model::M
end

@model function as_turing_model(model::RecordExpectedInfections, n::ModelShape)
    inner ~ to_submodel(with_expected_infections(model.model, n), false)
    exp_I_t := inner.exp_I_t
    return inner
end
