# [An SIR compartmental model](@id case-study-sir)

The renewal equation is one way to generate infections, but it is not the only
one. Mechanistic compartmental models describe transmission with a system of
ordinary differential equations (ODEs). [chatzilena2019contemporary](@citet)
showed how to embed such an ODE in a Bayesian model and infer its parameters,
using a classic influenza outbreak in an English boarding school as their
example.

This case study swaps the renewal infection process for an
[`ODEProcess`](@ref) built from [`SIRParams`](@ref), keeping the same
composable observation machinery. Infections come from solving the SIR equations
with the [SciML](https://sciml.ai) stack
[rackauckas2017differentialequations](@citep) rather than from a bespoke solver.

## The model

```math
\begin{aligned}
\frac{dS}{dt} &= -\beta S I, &
\frac{dI}{dt} &= \beta S I - \gamma I, &
\frac{dR}{dt} &= \gamma I, \\
R_0 &= \beta / \gamma, &
\lambda_t &= \mathrm{softplus}\!\big(N\, I(t)\big), &
y_t &\sim \mathrm{Poisson}(\lambda_t).
\end{aligned}
```

``S``, ``I``, ``R`` are population proportions; ``\beta`` is the transmission
rate, ``\gamma`` the recovery rate, and ``N`` the population size. The softplus
link smoothly scales the infected proportion to expected counts while staying
positive even if the solver returns a small negative value near zero.

## The infection process

[`SIRParams`](@ref) declares priors for the transmission rate, recovery rate,
and initial infected proportion, over a solver time span. We use weakly
informative priors that keep the basic reproduction number ``R_0 = \beta/\gamma``
in a plausible range for influenza and bounded away from the ``\gamma \to 0``
(``R_0 \to \infty``) singularity.

```@example sir
using EpiAwarePrototype, Distributions, Random, Turing, LogExpFunctions
Random.seed!(1978)

N = 763          # children in the school
n_days = 14

sir_params = SIRParams(
    tspan = (0.0, Float64(n_days)),
    infectiousness = LogNormal(-0.5, 0.5),
    recovery_rate = Gamma(8, 0.03125),
    initial_prop_infected = Beta(2, 200))
nothing # hide
```

[`ODEProcess`](@ref) composes those parameters with a solver and a `sol2infs`
link that pulls the infected compartment out of the ODE solution. This is the
standard SciML pattern — a problem definition composed with a solution method —
specialised to probabilistically sampled parameters. The default solver
switches between explicit and implicit methods, which keeps the solve robust
when the sampler proposes stiff parameter values.

```@example sir
sir_process = ODEProcess(
    params = sir_params,
    sol2infs = sol -> sol[2, :],
    solver_options = Dict(:saveat => 1.0))
nothing # hide
```

## The observation model

The ODE returns the infected proportion ``I(t)``; we scale it to counts with the
population size and a softplus transform using
[`TransformObservationModel`](@ref), then link to data with a
[`PoissonError`](@ref).

```@example sir
observation = TransformObservationModel(PoissonError(), x -> softplus.(N .* x))
nothing # hide
```

A compartmental model needs no time-varying latent ``R_t`` process — the
dynamics are fully determined by the ODE parameters — so the latent slot is
[`Null`](@ref). [`EpiAwareModel`](@ref) assembles the three parts exactly as in
the renewal examples.

```@example sir
model = EpiAwareModel(Null(), sir_process, observation)
nothing # hide
```

## Simulate and fit

Simulating from the prior produces an outbreak curve; fitting recovers the SIR
parameters. Differentiating through the ODE solution works with the default
(forward-mode) autodiff for a system this small.

```@example sir
sim = as_turing_model(model, fill(missing, n_days + 1), n_days + 1)()
y_obs = sim.generated_y_t
y_obs
```

```@example sir
chain = sample(as_turing_model(model, y_obs, n_days + 1), NUTS(), 100; progress = false)
nothing # hide
```

The posterior gives the transmission and recovery rates directly, and the basic
reproduction number is a deterministic function of them:

```@example sir
using MCMCChains, Statistics
mc = MCMCChains.Chains(chain)
β = vec(mc[:β])
γ = vec(mc[:γ])
R0 = β ./ γ
(β = mean(β), γ = mean(γ), R0 = mean(R0))
```

## Adding a stochastic ascertainment process

[chatzilena2019contemporary](@citet) also consider a stochastic variant where an
autoregressive process absorbs noise from model mis-specification. That is again
a composition: wrap the Poisson link in an [`Ascertainment`](@ref) modifier
carrying an [`AR`](@ref) latent process on the log scale, then re-apply the
population transform. No part of the infection model changes.

```@example sir
ascertainment = AR(
    damp_priors = [HalfNormal(0.005)],
    init_priors = [Normal(0, 0.001)],
    ϵ_t = HierarchicalNormal(std_prior = HalfNormal(0.02)))

stochastic_obs = TransformObservationModel(
    Ascertainment(model = PoissonError(), latent_model = ascertainment),
    x -> softplus.(N .* x))

stochastic_model = EpiAwareModel(Null(), sir_process, stochastic_obs)
length(rand(as_turing_model(stochastic_model, fill(missing, n_days + 1), n_days + 1)))
```

Swapping the deterministic observation model for the stochastic one is, once
more, a single structural change — the SIR infection process is reused
untouched.

## References

```@bibliography
Pages = ["sir-ode.md"]
Canonical = false
```
