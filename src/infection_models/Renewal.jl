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
is `data.transformation`, ``g_i`` is the discrete generation interval, and the
pre-window infections decay at the growth rate implied by ``\mathcal R_1``. The
``R_t`` process is generated *inside* the model, so `as_turing_model` takes only
the series length `n` and returns the named tuple `(; I_t, Z_t)` with `Z_t` the
(log) ``R_t`` path.

Renewal is the one infection model that needs a generation interval, so it alone
keeps an [`EpiData`](@ref) object.

## Fields

  - `data`: the [`EpiData`](@ref) object (generation interval + transformation).
  - `rt`: the latent process model (an [`AbstractLatentModel`](@ref)) generating
    the (log) reproduction number.
  - `initialisation_prior`: prior for the unconstrained initial infections.
  - `recurrent_step`: the renewal accumulation step (an
    [`AbstractConstantRenewalStep`](@ref)).

# Examples
```@example Renewal
using EpiAwarePrototype, Distributions
data = EpiData([0.2, 0.3, 0.5], exp)
renewal = Renewal(data; rt = RandomWalk(), initialisation_prior = Normal())
rand(as_turing_model(renewal, 20))
```
"
struct Renewal{E <: EpiData, L <: AbstractLatentModel, S <: Sampleable,
    A <: AbstractConstantRenewalStep} <: AbstractInfectionModel
    "`EpiData` object."
    data::E
    "Latent process model generating the (log) reproduction number."
    rt::L
    "Prior for the unconstrained initial infections."
    initialisation_prior::S
    "The renewal accumulation step."
    recurrent_step::A
end

function Renewal(data::EpiData; rt::AbstractLatentModel = RandomWalk(),
        initialisation_prior = Normal())
    recurrent_step = ConstantRenewalStep(reverse(data.gen_int))
    return Renewal(data, rt, initialisation_prior, recurrent_step)
end

function Renewal(; data::EpiData, rt::AbstractLatentModel = RandomWalk(),
        initialisation_prior = Normal())
    return Renewal(data; rt = rt, initialisation_prior = initialisation_prior)
end

# Initial renewal state from sampled I₀ and R₀, decaying at the implied rate.
function _make_renewal_init(epi_model::Renewal, I₀, Rt₀)
    r_approx = R_to_r(Rt₀, epi_model)
    return _renewal_init_state(
        epi_model.recurrent_step, I₀, r_approx, epi_model.data.len_gen_int)
end

@model function as_turing_model(epi_model::Renewal, n)
    Z_t ~ to_submodel(as_turing_model(epi_model.rt, n), false)
    init_incidence ~ epi_model.initialisation_prior
    I₀ = epi_model.data.transformation(init_incidence)
    Rt = epi_model.data.transformation.(Z_t)
    init = _make_renewal_init(epi_model, I₀, Rt[1])
    I_t = accumulate_scan(epi_model.recurrent_step, init, Rt)
    return (; I_t, Z_t)
end
