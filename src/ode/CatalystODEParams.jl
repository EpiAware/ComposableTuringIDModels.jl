# Public, model-agnostic Catalyst ODE parameter component.
#
# The struct and its docstring live here in `src/` so `CatalystODEParams` is a
# first-class, exported, `@ref`-able public component — not something users have
# to dig out of `Base.get_extension`. Its fields are fully type-parameterised and
# name no `Catalyst`/`ModelingToolkit` types, so this file has no heavy-stack
# dependency. Everything that DOES need Catalyst — the `ReactionSystem`
# constructor, the `as_turing_model` sampling method, and the `remake_ode_problem`
# hook — lives in the `EpiAwarePrototypeCatalystExt` extension and is only defined
# once a user runs `using Catalyst, ModelingToolkit`.

@doc raw"
Declarative, model-agnostic ODE parameter component built from **any** `Catalyst`
`ReactionSystem`, usable as the parameter component of an [`ODEProcess`](@ref) in
place of the hand-coded [`SIRParams`](@ref) / [`SEIRParams`](@ref).

You declare a reaction network and give priors for its initial conditions and
rate parameters; `Catalyst` + `ModelingToolkit` generate the ODE system **and a
symbolic Jacobian** (`jac = true`), so there is no hand-written vector field or
Jacobian to keep in sync, and nothing here is specialised to a particular
compartmental model. Construct it for an SIR network, an SEIR network, or any
other network the same way — only the reactions change.

Sampling and problem rebuilding are **symbolic**: `as_turing_model` samples each
supplied prior into a flat Turing variable named after its species / parameter
symbol (e.g. `β`, `S`) and returns symbolic `symbol => value` maps, which
`remake` places into the problem by name. There is no positional-index
bookkeeping, so species / parameter ordering inside the compiled problem is never
assumed. Index the resulting solution symbolically too, with the network's own
handles: `sol2infs = sol -> sol[rn.I, :]`.

!!! note \"Optional extension\"
    The constructor and sampling logic load only when `Catalyst` and
    `ModelingToolkit` are present (`using EpiAwarePrototype, Catalyst,
    ModelingToolkit`). The heavy symbolic stack stays out of the default install;
    the hand-coded models remain the zero-latency default. Constructing a
    `CatalystODEParams` before loading `Catalyst` raises an informative error.

# Arguments

  - `rn`: the `Catalyst` `ReactionSystem` (e.g. from `@reaction_network`).

# Keyword Arguments

  - `tspan`: the ODE solution time span.
  - `u0_priors`: the initial conditions, as symbolic-handle ⇒ spec pairs
    (`[rn.S => Beta(...), rn.R => 0.0, ...]`). Each spec is either a
    `Distribution` (sampled, as a flat variable named after the species) or a
    plain `Real` (a fixed initial value, not sampled). Every species of `rn` must
    appear.
  - `p_priors`: the rate parameters, as symbolic-handle ⇒ spec pairs
    (`[rn.β => LogNormal(...), ...]`), each a `Distribution` (sampled) or a fixed
    `Real`. Every parameter of `rn` must appear.

# Fields

  - `prob`: the `ODEProblem` built from `rn` (auto symbolic Jacobian).
  - `u0_specs`: per-species specs (symbolic handle, flat name, prior-or-fixed).
  - `p_specs`: per-parameter specs (symbolic handle, flat name, prior-or-fixed).

# Examples
```julia
using EpiAwarePrototype, Catalyst, ModelingToolkit, OrdinaryDiffEq, Distributions
sir = @reaction_network begin
    β, S + I --> 2I
    γ, I --> R
end
params = CatalystODEParams(sir;
    tspan = (0.0, 30.0),
    u0_priors = [sir.S => Beta(99, 1), sir.I => Beta(1, 99), sir.R => 0.0],
    p_priors = [sir.β => LogNormal(log(0.3), 0.1), sir.γ => LogNormal(log(0.1), 0.1)])
process = ODEProcess(params = params, sol2infs = sol -> sol[sir.I, :])
```
"
struct CatalystODEParams{P, U, R} <: AbstractLatentModel
    "The `ODEProblem` built from the reaction network (auto symbolic Jacobian)."
    prob::P
    "Per-species specs (symbolic handle, flat name, prior-or-fixed)."
    u0_specs::U
    "Per-parameter specs (symbolic handle, flat name, prior-or-fixed)."
    p_specs::R
end

# Fallback constructor: the real `ReactionSystem` method lives in the Catalyst
# extension. Reaching this means `Catalyst`/`ModelingToolkit` are not loaded.
function CatalystODEParams(rn; kw...)
    throw(ArgumentError(
        "CatalystODEParams requires the Catalyst extension. Run " *
        "`using Catalyst, ModelingToolkit` and pass a ReactionSystem."))
end
