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
dynamics are fully determined by the ODE parameters — so the [`ODEProcess`](@ref)
carries no latent process at all (its `Z_t` generated quantity is `nothing`).
[`EpiAwareModel`](@ref) assembles the infection and observation parts exactly as
in the renewal examples.

```@example sir
model = EpiAwareModel(sir_process, observation)
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

The deterministic model assumes the SIR equations describe the data exactly up to
Poisson counting noise. Real outbreaks rarely oblige: the compartmental model is
an approximation, and reporting intensity drifts over time. [chatzilena2019contemporary](@citet)
therefore also consider a stochastic variant in which a latent autoregressive
process on the log scale modulates the expected counts, absorbing variation the
mechanistic part cannot explain:

```math
\begin{aligned}
\kappa_t &= \rho\, \kappa_{t-1} + \epsilon_t, & \epsilon_t &\sim \mathrm{Normal}(0, \sigma), \\
\lambda_t &= \mathrm{softplus}\!\big(N\, I(t)\big)\,\exp(\kappa_t), &
y_t &\sim \mathrm{Poisson}(\lambda_t).
\end{aligned}
```

Setting ``\kappa_t = 0`` for all ``t`` recovers the deterministic model, so the
two are nested. In this package the ``\kappa_t`` process is exactly the [`AR`](@ref)
latent model already used for ``\log R_t`` in the renewal examples — here it
modulates the observation process rather than infections. An [`Ascertainment`](@ref)
modifier wraps the Poisson link and carries that latent process; the population
[`TransformObservationModel`](@ref) is re-applied on the outside. No part of the
infection model changes. The priors are weakly informative: damping near zero
(highly autocorrelated increments), an initial state near zero (no baseline
adjustment), and a small innovation standard deviation.

```@example sir
ascertainment = AR(
    damp_priors = [HalfNormal(0.005)],
    init_priors = [Normal(0, 0.001)],
    ϵ_t = HierarchicalNormal(std_prior = HalfNormal(0.02)))

stochastic_obs = TransformObservationModel(
    Ascertainment(model = PoissonError(), latent_model = ascertainment),
    x -> softplus.(N .* x))

stochastic_model = EpiAwareModel(sir_process, stochastic_obs)
nothing # hide
```

Swapping the deterministic observation model for the stochastic one is a single
structural change — the SIR infection process is reused untouched — and the
composed model is fit exactly as before. The ascertainment process adds latent
parameters, so we raise the NUTS target acceptance rate a little to keep the
sampler stable through the ODE solve.

```@example sir
stochastic_chain = sample(
    as_turing_model(stochastic_model, y_obs, n_days + 1),
    NUTS(0.9), 100; progress = false)
nothing # hide
```

The SIR parameters keep their flat names (`β`, `γ`, `I₀`); the ascertainment
process contributes its own block, prefixed `Ascertainment.` because modifiers
that introduce a named sub-process prefix their variables to keep them distinct.
The basic reproduction number is recovered as before, and the posterior for the
ascertainment innovation scale ``\sigma`` quantifies how much observation-level
noise the latent process absorbed:

```@example sir
stochastic_mc = MCMCChains.Chains(stochastic_chain)
βs = vec(stochastic_mc[:β])
γs = vec(stochastic_mc[:γ])
(R0 = mean(βs ./ γs),
    ascertainment_sigma = mean(vec(stochastic_mc[Symbol("Ascertainment.std")])))
```

Because the deterministic model is the ``\kappa_t = 0`` special case, the two
fits are directly comparable. Here the data were simulated from the deterministic
model, so the stochastic variant should recover a similar ``R_0`` while shrinking
its extra latent process towards zero — which is what the small fitted
``\sigma`` above and the close agreement below show:

```@example sir
(deterministic_R0 = mean(R0), stochastic_R0 = mean(βs ./ γs))
```

On a real, mis-specified outbreak the stochastic process would instead soak up
systematic departures from the SIR dynamics, guarding the mechanistic parameters
against that bias — the reason [chatzilena2019contemporary](@citet) introduce it.

## References

```@bibliography
Pages = ["sir-ode.md"]
Canonical = false
```
