# [Renewal model with negative-binomial reporting](@id case-study-renewal)

The renewal equation is the workhorse of real-time epidemic estimation: it
expresses new infections as a function of past infections weighted by the
generation interval, scaled by a time-varying reproduction number ``R_t``
[cori2013new](@citep). [mishra2020derivation](@citet) showed that this renewal
construction follows from an age-dependent branching process and pairs naturally
with a negative-binomial observation model to give a Bayesian hierarchical model
for reported case counts.

This case study builds that model from two composed parts — a [`Renewal`](@ref)
infection process that carries an autoregressive latent process for ``\log R_t``,
and a [`NegativeBinomialError`](@ref) observation model — simulates a synthetic
outbreak from it, and fits it back. The latent ``\log R_t`` process is *folded
into* the renewal model rather than supplied as a separate top-level component:
the reproduction number is the renewal model's own parameter process.

## The model

```math
\begin{aligned}
Z_t &= \rho\, Z_{t-1} + \epsilon_t, & \epsilon_t &\sim \mathrm{Normal}(0, \sigma), \\
R_t &= \exp(Z_t), \\
I_t &= R_t \sum_{s \ge 1} g_s\, I_{t-s}, \\
y_t &\sim \mathrm{NegBinomial}(I_t, \phi).
\end{aligned}
```

``g_s`` is the discretised generation interval, ``\rho`` the autoregressive
damping, ``\sigma`` the innovation standard deviation, and ``\phi`` the
observation overdispersion.

## Components

The latent process is a first-order autoregressive model on ``\log R_t`` with a
[`HierarchicalNormal`](@ref) innovation term. Strong autocorrelation in the
reproduction number is encoded by a damping prior concentrated near one. This
process is the renewal model's reproduction-number process — it is folded into
the infection model below rather than composed separately.

```@example renewal
using EpiAwarePrototype, Distributions, Random, Turing, Mooncake
using ADTypes: AutoMooncake
Random.seed!(1234)

latent = AR(
    damp_priors = [truncated(Normal(0.8, 0.05), 0, 1)],
    init_priors = [Normal(0.0, 0.25)],
    ϵ_t = HierarchicalNormal(std_prior = HalfNormal(0.1)))
```

The infection process needs a discrete generation interval. [`EpiData`](@ref)
takes a continuous distribution and discretises it with double interval
censoring [charniga2024best](@citep), using
[CensoredDistributions.jl](https://github.com/EpiAware/CensoredDistributions.jl).
Following [mishra2020derivation](@citet) we use a ``\mathrm{Gamma}(6.5, 0.62)``
serial interval as a proxy for the generation interval.

```@example renewal
data = EpiData(gen_distribution = Gamma(6.5, 0.62))
data.gen_int
```

The stored `gen_int` is a probability vector — the continuous serial interval
binned into daily weights that sum to one. Double interval censoring is not the
same as evaluating the continuous density at integer days: it accounts for both
the primary and secondary events falling anywhere within their days, which
shifts and spreads the mass relative to the underlying ``\mathrm{Gamma}``
[charniga2024best](@citep).

```@example renewal
sum(data.gen_int), length(data.gen_int)
```

The [`Renewal`](@ref) process couples that generation interval to the latent
``\log R_t`` process (its `rt` slot) and a prior for the initial infections.
[`Renewal`](@ref) is the only infection model that carries an [`EpiData`](@ref),
because it is the only one that uses a generation interval.

```@example renewal
renewal = Renewal(data; rt = latent, initialisation_prior = Normal(log(1.0), 0.25))
nothing # hide
```

## The infection process in isolation

Because the renewal model is a model in its own right, it can be exercised on its
own — without an observation model — and we can isolate the contribution of the
renewal equation by *pinning* its reproduction-number process to a known
trajectory. With the latent folded in, the way to do that is to build a renewal
model whose `rt` slot is a deterministic [`FixedIntercept`](@ref) latent, giving
a constant ``\log R_t``, and to fix the initial-infections parameter. The same
[`as_turing_model`](@ref) call that composes into the full model then runs the
infection process standalone, returning its infections `I_t` and the internal
latent draw `Z_t`.

```@example renewal
fixed_logR = log(1.4)
renewal_fixed = Renewal(data;
    rt = FixedIntercept(fixed_logR), initialisation_prior = Normal())
demo = fix(as_turing_model(renewal_fixed, 60), (init_incidence = 0.0,))()
(constant_Rt = round(exp(first(demo.Z_t)), digits = 2),
    grows = demo.I_t[end] > demo.I_t[1])
```

A constant ``R_t > 1`` grows incidence; a ``\log R_t`` path that declined through
zero would instead produce the textbook turn-over (incidence growing,
decelerating as ``R_t \to 1``, and falling once ``R_t < 1``). Driving the renewal
model with a richer fixed path is just a matter of swapping the
[`FixedIntercept`](@ref) latent for a deterministic latent of the desired shape.
Nothing here is conditioned on data; the component is inspected in isolation
before it is assembled into the full model with its sampled ``R_t`` process.

Reported cases are overdispersed counts of the latent infections. The prior is
placed on the cluster factor ``\sqrt{1/\phi}``, which is roughly the coefficient
of variation of the observation noise and so easier to reason about a priori.

```@example renewal
obs = NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1))
nothing # hide
```

[`EpiAwareModel`](@ref) assembles the two parts — the renewal infection process
(which already carries the latent ``R_t`` process) and the observation model —
into one composed model.

```@example renewal
model = EpiAwareModel(renewal, obs)
```

## Simulate

Passing `missing` observations turns the model into a prior simulator. The
composed model returns its generated quantities — the reported cases
`generated_y_t`, the latent infections `I_t`, and the latent process `Z_t`.

```@example renewal
n = 60
sim = as_turing_model(model, fill(missing, n), n)()
(; generated_y_t, I_t, Z_t) = sim
first(generated_y_t, 10)
```

The implied reproduction number is ``\exp(Z_t)``:

```@example renewal
extrema(exp.(Z_t))
```

We treat one simulated path as our observed data.

```@example renewal
y_obs = generated_y_t
sum(y_obs)
```

## Fit

Conditioning on the observed counts and sampling with NUTS recovers the
posterior. A short run keeps the page quick to build; the slightly raised target
acceptance rate keeps the sampler stable on the hierarchical innovation scale.
We differentiate with [Mooncake](https://chalk-lab.github.io/Mooncake.jl/), the
recommended backend for this package (see
[Automatic differentiation backend](@ref ad-backend)).

```@example renewal
posterior = as_turing_model(model, y_obs, n)
chain = sample(
    posterior, NUTS(0.9; adtype = AutoMooncake(; config = nothing)), 100;
    progress = false)
nothing # hide
```

Sampling returns a chain whose parameters keep their flat component names
(prefixing is disabled throughout the package). `sample` returns a
[FlexiChains](https://github.com/penelopeysm/FlexiChains.jl) chain, which
`summarystats` summarises directly — no conversion step — giving point estimates
*and* their uncertainty alongside the effective sample size and ``\hat{R}``
convergence diagnostic. The autoregressive damping ``\rho`` (`damp_AR[1]`), the
innovation scale ``\sigma`` (`std`), and the observation cluster factor
``\sqrt{1/\phi}`` (`cluster_factor`) are all identified from the single simulated
series:

```@example renewal
using MCMCChains
summarystats(chain)
```

The reproduction number ``R_t = \exp(Z_t)`` is a *generated quantity* rather than
a sampled parameter. [`generated_observables`](@ref) re-runs the fitted model
over the chain to recover the latent and infection trajectories per draw, from
which a posterior ``R_t`` band can be summarised:

```@example renewal
post = generated_observables(model, y_obs, chain)
typeof(post)
```

This page stops at the parameter summary rather than plotting the posterior
``R_t`` and posterior-predictive ``y_t`` ribbons: the docs build runs no plotting
stack, and turning the per-draw generated quantities into dated credible-interval
bands needs quantile-reduction helpers that would not earn their length here. The
trajectories are all present in `post` for a reader who wants them.

## Swap a component

Because the parts share one interface, an alternative observation assumption is
a one-line change. Swapping the negative-binomial reporting for a
[`PoissonError`](@ref) leaves the renewal infection process — and its latent
``R_t`` process — untouched:

```@example renewal
poisson_model = EpiAwareModel(renewal, PoissonError())
length(rand(as_turing_model(poisson_model, fill(missing, n), n)))
```

## References

```@bibliography
Pages = ["renewal-negbin.md"]
Canonical = false
```
