# SIR compartmental model parameters (latent component of an `ODEProcess`).

# Vector field of the density/per-capita SIR model.
function _sir_vf(du, u, p, t)
    S, Iv, R = u
    β, γ = p
    du[1] = -β * S * Iv
    du[2] = β * S * Iv - γ * Iv
    du[3] = γ * Iv
    return nothing
end

# Jacobian of the SIR vector field (speeds up stiff solves).
function _sir_jac(J, u, p, t)
    S, Iv, R = u
    β, γ = p
    J[1, 1] = -β * Iv
    J[1, 2] = -β * S
    J[2, 1] = β * Iv
    J[2, 2] = β * S - γ
    J[3, 2] = γ
    return nothing
end

const _sir_function = ODEFunction(_sir_vf; jac = _sir_jac)

@doc raw"
SIR compartmental model parameters and priors, usable as the latent component of
an [`ODEProcess`](@ref).

```math
\frac{dS}{dt} = -\beta S I, \quad
\frac{dI}{dt} = \beta S I - \gamma I, \quad
\frac{dR}{dt} = \gamma I
```

# Arguments

  - `params`: the [`SIRParams`](@ref) struct.
  - `n`: unused size argument (the ODE dimension is fixed); accepted for the
    common `as_turing_model` signature.

# Keyword Arguments

  - `tspan`: the ODE solution time span.
  - `infectiousness`: prior for ``\beta``.
  - `recovery_rate`: prior for ``\gamma``.
  - `initial_prop_infected`: prior for the initial infected proportion.

# Examples
```@example SIRParams
using ComposableTuringIDModels, OrdinaryDiffEq, Distributions
sirparams = SIRParams(
    tspan = (0.0, 30.0),
    infectiousness = LogNormal(log(0.3), 0.05),
    recovery_rate = LogNormal(log(0.1), 0.05),
    initial_prop_infected = Beta(1, 99))
rand(as_turing_model(sirparams, nothing))
```

## Fields

  - `prob`: the `ODEProblem` instance for the SIR model.
  - `infectiousness`: prior for ``\beta``.
  - `recovery_rate`: prior for ``\gamma``.
  - `initial_prop_infected`: prior for the initial infected proportion.
"
struct SIRParams{P <: ODEProblem, D <: Sampleable, E <: Sampleable, F <: Sampleable} <:
       AbstractLatentModel
    "The `ODEProblem` instance for the SIR model."
    prob::P
    "Prior for the infectiousness parameter."
    infectiousness::D
    "Prior for the recovery rate parameter."
    recovery_rate::E
    "Prior for the initial infected proportion."
    initial_prop_infected::F
end

function SIRParams(; tspan, infectiousness::Distribution, recovery_rate::Distribution,
        initial_prop_infected::Distribution)
    sir_prob = ODEProblem(_sir_function, [0.99, 0.01, 0.0], tspan)
    return SIRParams{typeof(sir_prob), typeof(infectiousness),
        typeof(recovery_rate), typeof(initial_prop_infected)}(
        sir_prob, infectiousness, recovery_rate, initial_prop_infected)
end

@model function as_turing_model(params::SIRParams, n)
    β ~ params.infectiousness
    γ ~ params.recovery_rate
    I₀ ~ params.initial_prop_infected
    u0 = [1.0 - I₀, I₀, 0.0]
    p = [β, γ]
    return (u0, p)
end
