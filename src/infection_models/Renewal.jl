# Renewal infection process model.

@doc raw"
Model unobserved infections via a time-varying renewal process.

```math
\mathcal R_t = g(Z_t), \qquad
I_t = \mathcal R_t \sum_{i=1}^{n-1} I_{t-i} g_i
```

where the latent path supplies (log) ``\mathcal R_t``, ``g`` is
`data.transformation`, ``g_i`` is the discrete generation interval, and the
pre-window infections decay at the growth rate implied by ``\mathcal R_1``.

# Arguments

  - `epi_model`: the [`Renewal`](@ref) model.
  - `_Rt`: the latent path of (log) reproduction numbers.

# Examples
```@example Renewal
using EpiAwarePrototype, Distributions
data = EpiData([0.2, 0.3, 0.5], exp)
renewal = Renewal(data; initialisation_prior = Normal())
rand(as_turing_model(renewal, randn(20) * 0.05))
```

## Fields

  - `data`: the [`EpiData`](@ref) object.
  - `initialisation_prior`: prior for the unconstrained initial infections.
  - `recurrent_step`: the renewal accumulation step (an
    [`AbstractConstantRenewalStep`](@ref)).
"
struct Renewal{E <: EpiData, S <: Sampleable, A <: AbstractConstantRenewalStep} <:
       AbstractInfectionModel
    "`EpiData` object."
    data::E
    "Prior for the unconstrained initial infections."
    initialisation_prior::S
    "The renewal accumulation step."
    recurrent_step::A
end

function Renewal(data::EpiData; initialisation_prior = Normal())
    recurrent_step = ConstantRenewalStep(reverse(data.gen_int))
    return Renewal(data, initialisation_prior, recurrent_step)
end

function Renewal(; data::EpiData, initialisation_prior = Normal())
    return Renewal(data; initialisation_prior = initialisation_prior)
end

# Initial renewal state from sampled I₀ and R₀, decaying at the implied rate.
function _make_renewal_init(epi_model::Renewal, I₀, Rt₀)
    r_approx = R_to_r(Rt₀, epi_model)
    return _renewal_init_state(
        epi_model.recurrent_step, I₀, r_approx, epi_model.data.len_gen_int)
end

@model function as_turing_model(epi_model::Renewal, _Rt)
    init_incidence ~ epi_model.initialisation_prior
    I₀ = epi_model.data.transformation(init_incidence)
    Rt = epi_model.data.transformation.(_Rt)
    init = _make_renewal_init(epi_model, I₀, Rt[1])
    return accumulate_scan(epi_model.recurrent_step, init, Rt)
end
