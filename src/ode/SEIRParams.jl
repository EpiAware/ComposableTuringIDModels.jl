# SEIR compartmental model parameters (latent component of an `ODEProcess`).

function _seir_vf(du, u, p, t)
    S, E, Iv, R = u
    β, α, γ = p
    du[1] = -β * S * Iv
    du[2] = β * S * Iv - α * E
    du[3] = α * E - γ * Iv
    du[4] = γ * Iv
    return nothing
end

function _seir_jac(J, u, p, t)
    S, E, Iv, R = u
    β, α, γ = p
    J[1, 1] = -β * Iv
    J[1, 3] = -β * S
    J[2, 1] = β * Iv
    J[2, 2] = -α
    J[2, 3] = β * S
    J[3, 2] = α
    J[3, 3] = -γ
    J[4, 3] = γ
    return nothing
end

const _seir_function = ODEFunction(_seir_vf; jac = _seir_jac)

@doc raw"
SEIR compartmental model parameters and priors, usable as the latent component of
an [`ODEProcess`](@ref).

```math
\frac{dS}{dt} = -\beta S I, \quad
\frac{dE}{dt} = \beta S I - \alpha E, \quad
\frac{dI}{dt} = \alpha E - \gamma I, \quad
\frac{dR}{dt} = \gamma I
```

The sampled initial infected proportion is split between the exposed and
infectious compartments using the constant-incidence equilibrium proportions
``\gamma/(\alpha+\gamma)`` and ``\alpha/(\alpha+\gamma)``.

# Arguments

  - `params`: the [`SEIRParams`](@ref) struct.
  - `n`: unused size argument; accepted for the common `as_turing_model`
    signature.

# Keyword Arguments

  - `tspan`: the ODE solution time span.
  - `infectiousness`: prior for ``\beta``.
  - `incubation_rate`: prior for ``\alpha``.
  - `recovery_rate`: prior for ``\gamma``.
  - `initial_prop_infected`: prior for the initial infected proportion.

# Examples
```@example SEIRParams
using EpiAwarePrototype, OrdinaryDiffEq, Distributions
seirparams = SEIRParams(
    tspan = (0.0, 30.0),
    infectiousness = LogNormal(log(0.3), 0.05),
    incubation_rate = LogNormal(log(0.1), 0.05),
    recovery_rate = LogNormal(log(0.1), 0.05),
    initial_prop_infected = Beta(1, 99))
rand(as_turing_model(seirparams, nothing))
```

## Fields

  - `prob`: the `ODEProblem` instance for the SEIR model.
  - `infectiousness`: prior for ``\beta``.
  - `incubation_rate`: prior for ``\alpha``.
  - `recovery_rate`: prior for ``\gamma``.
  - `initial_prop_infected`: prior for the initial infected proportion.
"
struct SEIRParams{P <: ODEProblem, D <: Sampleable, E <: Sampleable,
    F <: Sampleable, G <: Sampleable} <: AbstractLatentModel
    "The `ODEProblem` instance for the SEIR model."
    prob::P
    "Prior for the infectiousness parameter."
    infectiousness::D
    "Prior for the incubation rate parameter."
    incubation_rate::E
    "Prior for the recovery rate parameter."
    recovery_rate::F
    "Prior for the initial infected proportion."
    initial_prop_infected::G
end

function SEIRParams(; tspan, infectiousness::Distribution, incubation_rate::Distribution,
        recovery_rate::Distribution, initial_prop_infected::Distribution)
    seir_prob = ODEProblem(_seir_function, [0.99, 0.05, 0.05, 0.0], tspan)
    return SEIRParams{typeof(seir_prob), typeof(infectiousness),
        typeof(incubation_rate), typeof(recovery_rate),
        typeof(initial_prop_infected)}(seir_prob, infectiousness,
        incubation_rate, recovery_rate, initial_prop_infected)
end

@model function as_turing_model(params::SEIRParams, n)
    β ~ params.infectiousness
    α ~ params.incubation_rate
    γ ~ params.recovery_rate
    initial_infs ~ params.initial_prop_infected
    u0 = [1.0 - initial_infs, initial_infs * γ / (α + γ),
        initial_infs * α / (α + γ), 0.0]
    p = [β, α, γ]
    return (u0, p)
end
