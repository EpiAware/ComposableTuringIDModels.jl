# Renewal infection process model.

@doc raw"
Model unobserved infections via a time-varying renewal process driven by an
internally generated (log) reproduction number.

```math
Z_t \sim \text{latent}, \qquad
\mathcal R_t = g(Z_t), \qquad
I_t = \mathcal R_t \sum_{i=1}^{n-1} I_{t-i} g_i
```

where the latent model `rt` supplies the (log) reproduction number ``Z_t``, ``g``
is `transformation`, ``g_i`` is the discrete generation interval, and the
pre-window infections decay at the growth rate implied by ``\mathcal R_1``. The
``R_t`` process is generated *inside* the model, so `as_turing_model` takes only
the series length `n` and returns the named tuple `(; I_t, Z_t)` with `Z_t` the
(log) ``R_t`` path.

Renewal is the one infection model that needs a generation interval, so it takes
one directly (`gen_int`), either as a discrete probability vector or discretised
from a continuous distribution (see the keyword constructor).

`Renewal` is a step-composing helper: positional [`AbstractRenewalModifier`](@ref)
arguments are composed onto the renewal [`RenewalStep`](@ref). Passing a
[`SusceptibleDepletion`](@ref)`(N)` gives a renewal process with a fixed
population ``N`` and susceptible depletion
```math
I_t = \frac{S_{t-1}}{N} \mathcal R_t \sum_{i=1}^{n-1} I_{t-i} g_i, \qquad
S_t = S_{t-1} - I_t.
```

## Fields

  - `gen_int`: the discrete generation interval vector (non-negative, sums to 1).
  - `transformation`: the transformation between the unconstrained and
    constrained domains (default `exp`).
  - `rt`: the latent process model (an [`AbstractLatentModel`](@ref)) generating
    the (log) reproduction number.
  - `initialisation_prior`: prior for the unconstrained initial infections.
  - `recurrent_step`: the renewal accumulation step, a [`RenewalStep`](@ref)
    carrying the composed modifiers (none by default).

## Constructors

  - `Renewal(gen_int, modifiers...; rt, initialisation_prior, transformation)` —
    from a discrete generation interval vector (must be non-negative and sum to
    1), with any positional [`AbstractRenewalModifier`](@ref)s composed on top.
  - `Renewal(; gen_distribution, D_gen = nothing, Δd = 1.0, transformation = exp,
    rt, initialisation_prior, modifiers = ())` — discretise a continuous
    generation-interval distribution via double-interval censoring
    (CensoredDistributions.jl), with an optional tuple of `modifiers`.

# Examples
```@example Renewal
using ComposableTuringIDModels, Distributions
renewal = Renewal([0.2, 0.3, 0.5]; rt = RandomWalk(),
    initialisation_prior = Normal())
rand(as_turing_model(renewal, 20))

# With a fixed population and susceptible depletion.
depleting = Renewal([0.2, 0.3, 0.5], SusceptibleDepletion(1000.0);
    rt = RandomWalk(), initialisation_prior = Normal())
rand(as_turing_model(depleting, 20))
```
"
struct Renewal{T <: Real, F <: Function, L <: AbstractLatentModel,
    S <: Sampleable, A <: AbstractConstantRenewalStep} <: AbstractInfectionModel
    "Discrete generation interval."
    gen_int::Vector{T}
    "Transformation between unconstrained and constrained domains."
    transformation::F
    "Latent process model generating the (log) reproduction number."
    rt::L
    "Prior for the unconstrained initial infections."
    initialisation_prior::S
    "The renewal accumulation step."
    recurrent_step::A
end

function Renewal(gen_int::AbstractVector,
        modifiers::AbstractRenewalModifier...;
        rt::AbstractLatentModel = RandomWalk(),
        initialisation_prior = Normal(), transformation::Function = exp)
    @assert all(gen_int .>= 0) "Generation interval must be non-negative"
    @assert sum(gen_int)≈1 "Generation interval must sum to 1"
    core = ConstantRenewalStep(reverse(gen_int))
    recurrent_step = RenewalStep(core, modifiers)
    return Renewal(gen_int, transformation, rt, initialisation_prior,
        recurrent_step)
end

function Renewal(; gen_distribution::ContinuousDistribution, D_gen = nothing,
        Δd = 1.0, transformation::Function = exp,
        rt::AbstractLatentModel = RandomWalk(), initialisation_prior = Normal(),
        modifiers = ())
    # Drop the delay-0 bin (a generation interval has no mass at lag 0) and
    # renormalise, as the original EpiAware did.
    gen_int = _discretised_pmf(gen_distribution; Δd = Δd, D = D_gen) |>
              p -> p[2:end] ./ sum(p[2:end])
    return Renewal(gen_int, modifiers...; rt = rt,
        initialisation_prior = initialisation_prior,
        transformation = transformation)
end

# Initial renewal state from sampled I₀ and R₀, decaying at the implied rate.
function _make_renewal_init(infection::Renewal, I₀, Rt₀)
    r_approx = R_to_r(Rt₀, infection)
    return _renewal_init_state(
        infection.recurrent_step, I₀, r_approx, length(infection.gen_int))
end

@model function as_turing_model(infection::Renewal, n)
    Z_t ~ to_submodel(as_turing_model(infection.rt, n), false)
    init_incidence ~ infection.initialisation_prior
    I₀ = infection.transformation(init_incidence)
    Rt = infection.transformation.(Z_t)
    init = _make_renewal_init(infection, I₀, Rt[1])
    I_t = accumulate_scan(infection.recurrent_step, init, Rt)
    return (; I_t, Z_t)
end
