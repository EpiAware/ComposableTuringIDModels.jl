# [A declarative SEIR model with Catalyst](@id case-study-catalyst-seir)

The [SIR compartmental model](@ref case-study-sir) builds its dynamics from a
hand-written vector field. For that case study the package also hand-writes the
ODE's Jacobian, which the stiff/auto solver uses for speed and stability. Writing
a Jacobian by hand is error-prone and has to be re-derived, and kept in sync, for
every new compartmental model.

[Catalyst.jl](https://github.com/SciML/Catalyst.jl)
[loman2023catalyst](@citep) removes that step. You *declare* a reaction network,
and Catalyst together with
[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) generate the
ODE system **and a symbolic Jacobian** for you, kept consistent by construction.
The extension's `CatalystODEParams` is **model-agnostic**: it takes *any* Catalyst
network plus priors for its initial conditions and rate parameters, so the same
type builds an SIR model, the SEIR model below, or a network with a vaccinated
class or a second strain — only the reactions change.

!!! note "Optional extension"
    The Catalyst path lives in an optional package extension. It loads only when
    you add and import `Catalyst` (and `ModelingToolkit`) alongside
    `EpiAwarePrototype`, which keeps the heavy symbolic stack out of the default
    install. The hand-coded [`SIRParams`](@ref) / [`SEIRParams`](@ref) remain the
    zero-latency default; the declarative path is opt-in for users building new
    or more complex compartmental models.

## The reaction network

We declare SEIR as three reactions. Each reads like the transmission diagram: a
rate constant, then the species that react.

```math
\begin{aligned}
S + I &\xrightarrow{\beta} E + I & &\text{infection (}I\text{ catalytic)} \\
E &\xrightarrow{\alpha} I & &\text{incubation} \\
I &\xrightarrow{\gamma} R & &\text{recovery}
\end{aligned}
```

Catalyst expands these into exactly the SEIR drift used in the hand-coded model,

```math
\frac{dS}{dt} = -\beta S I, \quad
\frac{dE}{dt} = \beta S I - \alpha E, \quad
\frac{dI}{dt} = \alpha E - \gamma I, \quad
\frac{dR}{dt} = \gamma I,
```

and generates its Jacobian symbolically — there is no `_seir_jac` to write or
maintain.

```@example catalyst
using EpiAwarePrototype, Catalyst, ModelingToolkit, OrdinaryDiffEq
using Distributions, Random, Turing, LogExpFunctions, ADTypes
Random.seed!(1066)

seir = @reaction_network begin
    β, S + I --> E + I
    α, E --> I
    γ, I --> R
end
nothing # hide
```

## The infection process

Loading `Catalyst` activates the extension that backs `CatalystODEParams` (the
type itself is a normal, exported `EpiAwarePrototype` component). We hand it the
network, a solver time span, and priors for the initial conditions and rates. A
prior can be a `Distribution` (sampled, and named after its symbol in the chain —
`β`, `α`, `γ`) or a plain number (a fixed value, not sampled — here the
fully-susceptible start `S(0)` and the empty recovered class `R(0)`). The
`(u0, p)` sampling contract is the same as the hand-coded parameter models, so it
drops straight into an [`ODEProcess`](@ref).

```@example catalyst
N = 763          # children in the boarding school
n_days = 14

seir_params = CatalystODEParams(seir;
    tspan = (0.0, Float64(n_days)),
    u0_priors = [seir.S => 0.99, seir.E => Beta(2, 200),
        seir.I => Beta(2, 200), seir.R => 0.0],
    p_priors = [seir.β => LogNormal(-0.5, 0.4),
        seir.α => Gamma(8, 0.05), seir.γ => Gamma(8, 0.03125)])
nothing # hide
```

### A note on species ordering

Catalyst **sorts** the species and parameters when it builds the problem, so the
internal layout is generally not the order you wrote the network in (here the
state vector is `[S, I, E, R]` — the infectious compartment `I` is at index 2,
not 3). We never rely on that order: `CatalystODEParams` samples into **symbolic**
`symbol => value` maps that `remake` places by name, and the `sol2infs` link
indexes the solution **symbolically** with the network's own handle, so it pulls
the infectious compartment out by identity rather than by a hard-coded position:

```@example catalyst
seir_process = ODEProcess(
    params = seir_params,
    sol2infs = sol -> sol[seir.I, :],
    solver_options = Dict(:saveat => 1.0))
nothing # hide
```

## Composing and fitting

From here nothing is Catalyst-specific. We scale the infectious proportion to
expected counts with a population [`TransformObservationModel`](@ref) and a
[`PoissonError`](@ref), assemble with [`EpiAwareModel`](@ref), simulate from the
prior, and fit.

```@example catalyst
observation = TransformObservationModel(PoissonError(), x -> softplus.(N .* x))
model = EpiAwareModel(seir_process, observation)
nothing # hide
```

```@example catalyst
sim = as_turing_model(model, fill(missing, n_days + 1), n_days + 1)()
y_obs = sim.generated_y_t
y_obs
```

!!! warning "Use forward-mode autodiff for ODE models"
    The rest of these docs recommend Mooncake as the default AD backend, but ODE
    infection models are the exception: they sample under **ForwardDiff** today.
    Reverse-mode **Mooncake-driven NUTS through the ODE solver is not yet
    supported** — for the hand-coded *or* the Catalyst model — a pre-existing
    Turing/Mooncake/SciMLSensitivity integration gap rather than anything
    introduced by Catalyst. We therefore pass `AutoForwardDiff()` to NUTS
    explicitly.

```@example catalyst
chain = sample(
    as_turing_model(model, y_obs, n_days + 1),
    NUTS(; adtype = AutoForwardDiff()), 100; progress = false)
nothing # hide
```

`sample` returns a [FlexiChains](https://github.com/penelopeysm/FlexiChains.jl)
chain, which `summarystats` summarises directly — no conversion step — giving the
transmission, incubation, and recovery rates with their uncertainty alongside the
effective sample size and ``\hat R`` convergence diagnostic:

```@example catalyst
using MCMCChains
summarystats(chain)
```

The basic reproduction number ``R_0 = \beta / \gamma`` is a deterministic function
of the rates, formed by indexing the chain for the ``\beta`` and ``\gamma`` draws:

```@example catalyst
using Statistics
β = vec(chain[@varname(β)])
γ = vec(chain[@varname(γ)])
(R0 = mean(β ./ γ),)
```

Adding a fourth compartment, a vaccinated class, or a second strain is now a
matter of writing a different reaction network and passing it to the same
`CatalystODEParams` — the vector field and Jacobian follow automatically, rather
than re-deriving a Jacobian by hand. That is the trade the Catalyst extension
offers: a one-off symbolic-compilation cost and a heavier dependency tree, in
exchange for declarative, model-agnostic, self-consistent dynamics.

## References

```@bibliography
Pages = ["catalyst-seir.md"]
Canonical = false
```
