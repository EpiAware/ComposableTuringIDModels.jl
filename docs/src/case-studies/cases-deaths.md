# [Joint cases and deaths from shared infections](@id case-study-cases-deaths)

A single surveillance stream rarely tells the whole story. During an epidemic the
same latent infections are seen through several lenses at once: confirmed cases,
hospital admissions, deaths. Each stream reports a different *fraction* of
infections, after a different *delay*, with a different *noise* characteristic —
cases are plentiful but ascertainment-dependent, deaths are rarer and lag much
further behind but are comparatively reliably recorded. Estimating the
reproduction number from several streams together uses more of the available
signal and guards against the biases of any single stream — different surveillance
streams can imply materially different reproduction numbers, so reconciling them
is part of the modelling task [sherratt2021exploring](@citep).

This case study fits one [`Renewal`](@ref) infection process to **two**
observation streams at once. Both streams watch the *same* latent infections
``I_t``; each carries its own reporting delay, ascertainment fraction, and
observation-error model. The two pipelines are composed into one model with
[`StackObservationModels`](@ref) — the package's construct for observing a shared
latent signal through several named streams in parallel. We simulate a synthetic
outbreak from the joint model, fit it back with Turing, and recover the shared
reproduction-number process together with each stream's own parameters.

## The model

```math
\begin{aligned}
Z_t &= \rho\, Z_{t-1} + \epsilon_t, & \epsilon_t &\sim \mathrm{Normal}(0, \sigma), \\
R_t &= \exp(Z_t), \qquad
I_t = R_t \sum_{s \ge 1} g_s\, I_{t-s}, \\[2pt]
C_t &= \alpha_{\mathrm c} \sum_{s \ge 1} d^{\mathrm c}_s\, I_{t-s}, &
y^{\mathrm c}_t &\sim \mathrm{NegBinomial}(C_t, \phi_{\mathrm c}), \\
D_t &= \alpha_{\mathrm d} \sum_{s \ge 1} d^{\mathrm d}_s\, I_{t-s}, &
y^{\mathrm d}_t &\sim \mathrm{NegBinomial}(D_t, \phi_{\mathrm d}).
\end{aligned}
```

A shared infection process (top line) feeds two parallel observation pipelines.
``\alpha_{\mathrm c}`` is the case-ascertainment ratio and ``\alpha_{\mathrm d}``
the infection-fatality ratio; ``d^{\mathrm c}`` and ``d^{\mathrm d}`` are the
infection-to-case and infection-to-death delay kernels; ``\phi_{\mathrm c}`` and
``\phi_{\mathrm d}`` are stream-specific overdispersions. Crucially ``R_t``,
``I_t``, ``\rho`` and ``\sigma`` are **shared** — both streams inform them — while
everything subscripted by a stream is its own.

## The shared infection process

The infection core is the renewal model of the
[first case study](@ref case-study-renewal): an autoregressive ``\log R_t``
process folded into a [`Renewal`](@ref) process driven by a discretised
generation interval. Nothing about it is multi-stream — it generates one latent
infection trajectory, which the observation layer then views twice.

```@example casesdeaths
using EpiAwarePrototype, Distributions, Random, Turing, Mooncake
using ADTypes: AutoMooncake
Random.seed!(20240515)

latent = AR(
    damp_priors = [truncated(Normal(0.8, 0.05), 0, 1)],
    init_priors = [Normal(0.0, 0.25)],
    ϵ_t = HierarchicalNormal(std_prior = HalfNormal(0.1)))

data = EpiData(gen_distribution = Gamma(6.5, 0.62))
renewal = Renewal(data;
    rt = latent, initialisation_prior = Normal(log(1000.0), 0.25))
nothing # hide
```

The initialisation prior is centred at a thousand infections rather than one: we
want an outbreak large enough that even the rare death stream carries information,
since a low infection-fatality ratio turns thousands of infections into only tens
of deaths.

## Two observation pipelines

Each stream is a full observation model in its own right — the same kind of
layered object built in the [delays case study](@ref case-study-delays), one per
stream. We build them from the inside out.

### Cases

Cases are a large fraction of infections reported after a short delay. Reading the
composition from the inside:

  - a [`NegativeBinomialError`](@ref) link, because reported counts are
    overdispersed relative to Poisson;
  - an [`Ascertainment`](@ref) modifier scaling expected infections down to the
    reported fraction. Here the case-ascertainment ratio is *known* (around 60%),
    so the latent is a deterministic [`FixedIntercept`](@ref) on the log scale and
    we pass `latent_prefix = ""` to keep it from introducing a named sub-process —
    it contributes no sampled parameters;
  - a [`LatentDelay`](@ref) convolving expected infections with a short
    infection-to-report delay, discretised from a continuous distribution by
    double interval censoring [charniga2024best](@citep).

```@example casesdeaths
cases_obs = LatentDelay(
    Ascertainment(
        NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1)),
        FixedIntercept(log(0.6)); latent_prefix = ""),
    LogNormal(1.6, 0.5))
```

Every component has a `show` method, so the constructed object prints its
structure: the outer [`LatentDelay`](@ref) tree, the wrapped
[`Ascertainment`](@ref)/[`NegativeBinomialError`](@ref) it carries, and the
discretised delay kernel (`rev_pmf`) that the convolution will use.

### Deaths

Deaths share that structure but with three differences that make them a genuinely
distinct view of the same infections: a much smaller ascertainment fraction (an
infection-fatality ratio of roughly 1.5%), a longer and more dispersed
infection-to-death delay, and — as a consequence of being rarer — a stream that
is far less informative on its own.

```@example casesdeaths
deaths_obs = LatentDelay(
    Ascertainment(
        NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1)),
        FixedIntercept(log(0.015)); latent_prefix = ""),
    LogNormal(2.8, 0.4))
```

The same tree prints for the death stream — note its longer, more spread-out
delay kernel and the smaller fixed ascertainment intercept (``\log 0.015``).

Fixing the two ascertainment fractions rather than inferring them is a deliberate
identifiability choice. A single shared infection trajectory observed through two
*unknown* multiplicative fractions is only identified up to an overall scale; with
the case-ascertainment ratio and infection-fatality ratio supplied as known
external quantities, the absolute infection level — and hence ``R_t`` — is pinned.
Replacing either [`FixedIntercept`](@ref) with a prior is a one-line change if a
stream's fraction is to be learned alongside the rest.

## Stacking the streams

[`StackObservationModels`](@ref) takes a `NamedTuple` of observation models and
applies each to the *same* expected-infection series, keeping every stream's
variables distinct by prefixing them with the stream name. This is the parallel
multi-stream construct: one shared latent signal in, several named observed
series out.

```@example casesdeaths
stacked = StackObservationModels((cases = cases_obs, deaths = deaths_obs))
```

The printed `model_names` (`["cases", "deaths"]`) are the keys that prefix each
stream's parameters and key its data.

[`EpiAwareModel`](@ref) then assembles the whole model from exactly two parts, as
in every other case study — the infection process and *an* observation model. The
observation model just happens to be a stack of two pipelines rather than a single
one; the composition does not change.

```@example casesdeaths
model = EpiAwareModel(renewal, stacked)
```

## Simulate

The stacked observation model expects its data as a `NamedTuple` with one entry
per stream, matching the names given to [`StackObservationModels`](@ref). Passing
`missing` for both turns the model into a joint prior simulator.

```@example casesdeaths
n = 70
y_missing = (cases = fill(missing, n), deaths = fill(missing, n))
sim = as_turing_model(model, y_missing, n)()
(; I_t, Z_t) = sim
extrema(exp.(Z_t))   # the simulated reproduction-number range
```

The simulator returns one series per stream in `generated_y_t`. Each
[`LatentDelay`](@ref) shortens its own expected series by the length of its delay
kernel, so the leading entries — for which no fully observed expectation exists —
come back as `missing`; the deaths stream, with the longer delay, has more leading
`missing`s than cases. The counts themselves reflect the two ascertainment
fractions: many cases, far fewer deaths.

```@example casesdeaths
cases_obs_series = sim.generated_y_t[1]
deaths_obs_series = sim.generated_y_t[2]
(cases_total = sum(skipmissing(cases_obs_series)),
    deaths_total = sum(skipmissing(deaths_obs_series)))
```

We take the two simulated paths as our observed data, packaged in the same
named-tuple shape the model expects.

```@example casesdeaths
y_obs = (cases = cases_obs_series, deaths = deaths_obs_series)
nothing # hide
```

## Fit

Conditioning on both streams and sampling with NUTS recovers the joint posterior.
We differentiate with [Mooncake](https://chalk-lab.github.io/Mooncake.jl/), the
recommended backend for this package, passed to NUTS through
[`AutoMooncake`](https://github.com/SciML/ADTypes.jl); this reverse-mode backend
scales better than the default forward mode as the latent ``\log R_t`` vector
grows with the series length. A short run keeps the page quick to build, and the
slightly raised target acceptance rate keeps the sampler stable on the
hierarchical innovation scale.

```@example casesdeaths
posterior = as_turing_model(model, y_obs, n)
chain = sample(
    posterior, NUTS(0.9; adtype = AutoMooncake(; config = nothing)), 100;
    progress = false)
nothing # hide
```

`sample` returns a [FlexiChains](https://github.com/penelopeysm/FlexiChains.jl)
chain, which we index by variable name directly — no conversion step. The stack
prefixes each stream's own parameters with the stream name, so the two
overdispersions are reached as `@varname(cases.cluster_factor)` and
`@varname(deaths.cluster_factor)` and stay distinct. The shared infection
parameters — the autoregressive damping ``\rho`` (`damp_AR[1]`), the innovation
scale ``\sigma`` (`std`), and the initial infections (`init_incidence`) — keep
their flat names, because they belong to the single shared infection process
rather than to either stream:

```@example casesdeaths
using Turing: @varname
using Statistics

posterior_draws(vn) = vec(chain[vn])
posterior_summary(vn) = (mean = mean(posterior_draws(vn)),
    std = std(posterior_draws(vn)))

(init_incidence = posterior_summary(@varname(init_incidence)),
    damp = posterior_summary(@varname(damp_AR[1])),
    sigma = posterior_summary(@varname(std)),
    cases_cluster = posterior_summary(@varname(cases.cluster_factor)),
    deaths_cluster = posterior_summary(@varname(deaths.cluster_factor)))
```

The shared parameters are recovered from the two streams jointly: `init_incidence`
sits close to its simulated ``\log 1000 \approx 6.9``, and the damping and
innovation scale match the priors that generated the data. The two cluster factors
are estimated *separately* and tell the two-stream story directly — the
case-stream overdispersion is pinned down tightly by thousands of counts (a small
posterior `std`), while the death-stream overdispersion has a visibly wider
posterior because the sparse death counts carry less information. This is the
payoff of the stack: each stream contributes what it can to the shared infection
process while keeping its own noise model.

The reproduction number ``R_t = \exp(Z_t)`` remains a *generated quantity* of the
shared infection process, exactly as in the single-stream case.
[`generated_observables`](@ref) bundles the fitted model, data, and chain so the
per-draw latent and infection trajectories can be recovered for a posterior
``R_t`` band:

```@example casesdeaths
post = generated_observables(model, y_obs, chain)
typeof(post)
```

As in the other sampling case studies we stop at the parameter summary rather than
plotting ribbons: the docs build runs no plotting stack, and the trajectories are
all available in `post` for a reader who wants them.

## Adding a third stream

Because the streams are composed rather than hard-wired, a third surveillance
signal — hospital admissions, say — is one more entry in the named tuple. Its
pipeline is built exactly like the other two (its own delay, ascertainment, and
error model), and nothing about the shared infection process changes:

```@example casesdeaths
admissions_obs = LatentDelay(
    Ascertainment(
        NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1)),
        FixedIntercept(log(0.08)); latent_prefix = ""),
    LogNormal(2.2, 0.45))

three_stream = StackObservationModels((
    cases = cases_obs, admissions = admissions_obs, deaths = deaths_obs))
three_model = EpiAwareModel(renewal, three_stream)
length(rand(as_turing_model(three_model,
    (cases = fill(missing, n), admissions = fill(missing, n),
        deaths = fill(missing, n)), n)))
```

## Parallel versus sequential composition

[`StackObservationModels`](@ref) composes streams **in parallel**: every stream
observes the same expected-infection series independently, and the streams never
interact. That is the right structure for cases and deaths when both are read as
fractions of *infections*, which is the model fitted above.

A different, equally common assumption is **sequential** structure, where one
stream arises *downstream of* another — for example deaths modelled as a delayed
fraction of *reported cases* rather than of infections, so that whatever is
reflected in the case series (a reporting artefact, an ascertainment dip)
propagates into the death series. Expressing that needs a stream's *expected
output* to feed the next stream's *expected input*, which the parallel stack does
not do — it broadcasts one shared expected series to all streams and never threads
one stream's expectation into another.

That ordered cascade is a genuinely different construct rather than a parameter of
the stack, so the two are provided as a matched pair: the parallel
[`StackObservationModels`](@ref) demonstrated here, and a sequential
`SequentialObservationModels` (the cascade `infections → stream 1 → stream 2 → …`)
added separately in
[#58](https://github.com/EpiAware/EpiAwarePrototype.jl/pull/58), which resolves the
follow-up issue [#51](https://github.com/EpiAware/EpiAwarePrototype.jl/issues/51)
this page originally filed. Both share the per-stream pipelines built above — each
stream is the same delay/ascertainment/error object either way; only how the
streams are wired to the shared infection signal differs. The parallel model here
covers the common cases-and-deaths-from-infections setting; reach for the
sequential construct when a downstream stream should inherit an upstream stream's
reporting structure.

## References

```@bibliography
Pages = ["cases-deaths.md"]
Canonical = false
```
