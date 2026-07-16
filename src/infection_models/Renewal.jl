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
keeps an [`IDData`](@ref) object.

## Fields

  - `data`: the [`IDData`](@ref) object (generation interval + transformation).
  - `rt`: the latent process model (an [`AbstractLatentModel`](@ref)) generating
    the (log) reproduction number.
  - `initialisation`: prior for the unconstrained initial infections (a
    `Distribution` or prior model, sampled through [`as_turing_submodel`](@ref)).
  - `recurrent_step`: the renewal accumulation step (an
    [`AbstractConstantRenewalStep`](@ref)).

# Examples
```@example Renewal
using ComposableTuringIDModels, Distributions
data = IDData([0.2, 0.3, 0.5], exp)
renewal = Renewal(data; rt = RandomWalk(), initialisation = Normal())
rand(as_turing_model(renewal, 20))
```
"
struct Renewal{E <: IDData, L <: PriorLike, S <: PriorLike,
    A <: AbstractConstantRenewalStep} <: AbstractInfectionModel
    "`IDData` object."
    data::E
    "Latent process model generating the (log) reproduction number."
    rt::L
    "Prior for the unconstrained initial infections."
    initialisation::S
    "The renewal accumulation step."
    recurrent_step::A
end

function Renewal(data::IDData; rt = RandomWalk(),
        initialisation = Normal())
    recurrent_step = ConstantRenewalStep(reverse(data.gen_int))
    return Renewal(data, rt, initialisation, recurrent_step)
end

function Renewal(; data::IDData, rt = RandomWalk(),
        initialisation = Normal())
    return Renewal(data; rt = rt, initialisation = initialisation)
end

# Initial renewal state from sampled I₀ and R₀, decaying at the implied rate.
function _make_renewal_init(infection::Renewal, I₀, Rt₀)
    r_approx = R_to_r(Rt₀, infection)
    return _renewal_init_state(
        infection.recurrent_step, I₀, r_approx, infection.data.len_gen_int)
end

@model function as_turing_model(infection::Renewal, n)
    Z_t ~ as_turing_submodel(infection.rt, n)
    init_incidence ~ as_turing_submodel(
        infection.initialisation, 1; prefix = true)
    I₀ = infection.data.transformation(only(init_incidence))
    Rt = infection.data.transformation.(Z_t)
    init = _make_renewal_init(infection, I₀, Rt[1])
    I_t = accumulate_scan(infection.recurrent_step, init, Rt)
    return (; I_t, Z_t)
end
