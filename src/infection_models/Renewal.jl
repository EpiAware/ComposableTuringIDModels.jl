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
a discrete probability vector is used as-is, a continuous `Distribution` is
discretised internally (see the constructor), and a pmf-producing prior model
(e.g. an [`UncertainDelay`](@ref)) lets the generation interval itself be
**inferred** — its distribution's parameters carry priors and the interval is
rediscretised per draw through the [`as_turing_submodel`](@ref) seam.

`Renewal` is a step-composing helper: positional [`AbstractRenewalModifier`](@ref)
arguments are composed onto the renewal [`RenewalStep`](@ref). Passing a
[`SusceptibleDepletion`](@ref)`(N)` gives a renewal process with a fixed
population ``N`` and susceptible depletion
```math
I_t = \frac{S_{t-1}}{N} \mathcal R_t \sum_{i=1}^{n-1} I_{t-i} g_i, \qquad
S_t = S_{t-1} - I_t.
```

## Fields

  - `gen_int`: the discrete generation interval vector (non-negative, sums to 1),
    or, for an inferred generation interval, the pmf-producing prior model.
  - `transformation`: the transformation between the unconstrained and
    constrained domains (default `exp`).
  - `rt`: the latent process model (an [`AbstractLatentModel`](@ref)) generating
    the (log) reproduction number. A length-`n` PATH slot: a bare `Distribution`
    here is auto-wrapped in an [`Intercept`](@ref), giving a constant path (one
    shared draw broadcast to length `n`); use [`IID`](@ref) for `n` independent
    draws.
  - `initialisation`: prior for the unconstrained initial infections (a
    `Distribution` or prior model, sampled through [`as_turing_submodel`](@ref)).
  - `recurrent_step`: the renewal accumulation step (an
    [`AbstractConstantRenewalStep`](@ref)), or `nothing` when the generation
    interval is inferred and the step is built per draw.

## Constructor

  - `Renewal(; generation_time, rt, initialisation, transformation = exp,
    D_gen = nothing, Δd = 1.0)` — one keyword constructor that dispatches on
    `generation_time`:

      + a discrete probability **vector** (non-negative, sums to 1) is used
        directly as the generation interval;
      + a continuous **`Distribution`** is discretised via double-interval
        censoring (CensoredDistributions.jl), using `D_gen`/`Δd`, with the
        delay-0 bin dropped and the remainder renormalised; and
      + a pmf-producing **prior model** (an [`AbstractPriorModel`](@ref) such as
        an [`UncertainDelay`](@ref)) is held as-is and sampled per draw, giving an
        **inferred** generation interval (an uncertain discretised distribution
        used as the generation interval). Its fixed horizon keeps the interval
        length constant across draws; the lag-0 bin is dropped and the remainder
        renormalised per draw, exactly as for the fixed distribution.

# Examples

A fixed generation interval:

```@example Renewal
using ComposableTuringIDModels, Distributions
renewal = Renewal(; generation_time = [0.2, 0.3, 0.5], rt = RandomWalk(),
    initialisation = Normal())
rand(as_turing_model(renewal, 20))
```

An inferred generation interval — an uncertain discretised distribution used as
the generation interval, whose `LogNormal` parameters carry priors:

```@example Renewal
gen = UncertainDelay(
    LogNormal, [Normal(1.9, 0.2), truncated(Normal(0.5, 0.2), 0, Inf)]; D = 14.0)
renewal = Renewal(; generation_time = gen, rt = RandomWalk(),
    initialisation = Normal())
rand(as_turing_model(renewal, 20))

# With a fixed population and susceptible depletion.
depleting = Renewal([0.2, 0.3, 0.5], SusceptibleDepletion(1000.0);
    rt = RandomWalk(), initialisation_prior = Normal())
rand(as_turing_model(depleting, 20))
```
"
struct Renewal{G, F <: Function, L <: PriorLike, S <: PriorLike, A} <:
       AbstractInfectionModel
    "Discrete generation interval, or a pmf-producing prior model when inferred."
    gen_int::G
    "Transformation between unconstrained and constrained domains."
    transformation::F
    "Latent process model generating the (log) reproduction number."
    rt::L
    "Prior for the unconstrained initial infections."
    initialisation::S
    "The renewal accumulation step (`nothing` when the interval is inferred)."
    recurrent_step::A
end

function Renewal(; generation_time, rt = RandomWalk(),
        initialisation = Normal(), transformation::Function = exp,
        D_gen = nothing, Δd = 1.0)
    gen_int, recurrent_step = _renewal_fields(
        generation_time; D_gen = D_gen, Δd = Δd)
    return Renewal(gen_int, transformation, _path_prior(rt), initialisation,
        recurrent_step)
end

# Fixed generation interval (a pmf vector or a continuous distribution): bake the
# discretised interval and its reversed renewal step at construction, exactly as
# before.
function _renewal_fields(generation_time; D_gen = nothing, Δd = 1.0)
    gen_int = _renewal_gen_int(generation_time; D_gen = D_gen, Δd = Δd)
    return gen_int, ConstantRenewalStep(reverse(gen_int))
end

# Inferred generation interval: hold the pmf-producing prior model and build the
# renewal step per draw inside the `@model`, so no interval or step is baked.
function _renewal_fields(generation_time::AbstractPriorModel; D_gen = nothing,
        Δd = 1.0)
    return generation_time, nothing
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

# Initial renewal state from sampled I₀ and R₀, decaying at the implied rate. The
# generation interval and its step are passed in so the fixed (baked) and
# inferred (per-draw) paths share one initialiser.
function _make_renewal_init(step::AbstractConstantRenewalStep, gen_int, I₀, Rt₀)
    r_approx = R_to_r(Rt₀, gen_int)
    return _renewal_init_state(step, I₀, r_approx, length(gen_int))
end

@model function as_turing_model(infection::Renewal, n)
    Z_t ~ as_turing_submodel(infection.rt, n)
    init_incidence ~ as_turing_submodel(
        infection.initialisation, 1; prefix = true)
    I₀ = infection.transformation(only(init_incidence))
    Rt = infection.transformation.(Z_t)

    # The generation interval is either fixed (its baked renewal step used
    # directly) or inferred: a pmf-producing prior model sampled through the
    # single seam, with the lag-0 bin dropped and the remainder renormalised per
    # draw (the generation-interval convention the fixed distribution path applies
    # at construction). The renewal step is built per draw so the gradient flows
    # through the discretisation.
    if infection.gen_int isa AbstractPriorModel
        gen ~ as_turing_submodel(infection.gen_int; prefix = true)
        gen_int = gen[2:end] ./ sum(gen[2:end])
        step = ConstantRenewalStep(reverse(gen_int))
    else
        gen_int = infection.gen_int
        step = infection.recurrent_step
    end

    init = _make_renewal_init(step, gen_int, I₀, Rt[1])
    I_t = accumulate_scan(step, init, Rt)
    return (; I_t, Z_t)
end
