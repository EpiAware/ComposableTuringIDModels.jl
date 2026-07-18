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
one directly through the `generation_time` keyword, which dispatches on the value:
a discrete probability vector is used as-is, while a continuous `Distribution` is
discretised internally (see the constructor).

## Fields

  - `gen_int`: the discrete generation interval vector (non-negative, sums to 1).
  - `transformation`: the transformation between the unconstrained and
    constrained domains (default `exp`).
  - `rt`: the latent process model (an [`AbstractLatentModel`](@ref)) generating
    the (log) reproduction number.
  - `initialisation`: prior for the unconstrained initial infections (a
    `Distribution` or prior model, sampled through [`as_turing_submodel`](@ref)).
  - `recurrent_step`: the renewal accumulation step (an
    [`AbstractConstantRenewalStep`](@ref)).

## Constructor

  - `Renewal(; generation_time, rt, initialisation, transformation = exp,
    D_gen = nothing, Δd = 1.0)` — one keyword constructor that dispatches on
    `generation_time`:

      + a discrete probability **vector** (non-negative, sums to 1) is used
        directly as the generation interval; and
      + a continuous **`Distribution`** is discretised via double-interval
        censoring (CensoredDistributions.jl), using `D_gen`/`Δd`, with the
        delay-0 bin dropped and the remainder renormalised.

# Examples
```@example Renewal
using ComposableTuringIDModels, Distributions
renewal = Renewal(; generation_time = [0.2, 0.3, 0.5], rt = RandomWalk(),
    initialisation = Normal())
rand(as_turing_model(renewal, 20))
```
"
struct Renewal{T <: Real, F <: Function, L <: PriorLike, S <: PriorLike,
    A <: AbstractConstantRenewalStep} <: AbstractInfectionModel
    "Discrete generation interval."
    gen_int::Vector{T}
    "Transformation between unconstrained and constrained domains."
    transformation::F
    "Latent process model generating the (log) reproduction number."
    rt::L
    "Prior for the unconstrained initial infections."
    initialisation::S
    "The renewal accumulation step."
    recurrent_step::A
end

function Renewal(; generation_time, rt = RandomWalk(),
        initialisation = Normal(), transformation::Function = exp,
        D_gen = nothing, Δd = 1.0)
    gen_int = _renewal_gen_int(generation_time; D_gen = D_gen, Δd = Δd)
    recurrent_step = ConstantRenewalStep(reverse(gen_int))
    return Renewal(gen_int, transformation, rt, initialisation,
        recurrent_step)
end

# `generation_time` as a discrete PMF: use it directly (must be a valid pmf).
function _renewal_gen_int(gen_int::AbstractVector; D_gen = nothing, Δd = 1.0)
    @assert all(gen_int .>= 0) "Generation interval must be non-negative"
    @assert sum(gen_int)≈1 "Generation interval must sum to 1"
    return collect(gen_int)
end

# `generation_time` as a continuous distribution: discretise via double-interval
# censoring, drop the delay-0 bin (a generation interval has no mass at lag 0) and
# renormalise, as the original EpiAware did.
function _renewal_gen_int(gen_distribution::ContinuousDistribution;
        D_gen = nothing, Δd = 1.0)
    return _discretised_pmf(gen_distribution; Δd = Δd, D = D_gen) |>
           p -> p[2:end] ./ sum(p[2:end])
end

# Initial renewal state from sampled I₀ and R₀, decaying at the implied rate.
function _make_renewal_init(infection::Renewal, I₀, Rt₀)
    r_approx = R_to_r(Rt₀, infection)
    return _renewal_init_state(
        infection.recurrent_step, I₀, r_approx, length(infection.gen_int))
end

@model function as_turing_model(infection::Renewal, n)
    Z_t ~ as_turing_submodel(infection.rt, n)
    init_incidence ~ as_turing_submodel(
        infection.initialisation, 1; prefix = true)
    I₀ = infection.transformation(only(init_incidence))
    Rt = infection.transformation.(Z_t)
    init = _make_renewal_init(infection, I₀, Rt[1])
    I_t = accumulate_scan(infection.recurrent_step, init, Rt)
    return (; I_t, Z_t)
end
