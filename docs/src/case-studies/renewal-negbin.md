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
and a [`NegativeBinomialError`](@ref) observation model — and fits it to the
test-confirmed COVID-19 cases from South Korea that [mishra2020derivation](@citet)
analysed. The latent ``\log R_t`` process is *folded into* the renewal model
rather than supplied as a separate top-level component: the reproduction number
is the renewal model's own parameter process.

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

The latent process is a second-order autoregressive model on ``\log R_t`` with a
[`HierarchicalNormal`](@ref) innovation term, matching [mishra2020derivation](@citet).
Strong autocorrelation in the reproduction number is encoded by a first damping
prior concentrated near one (``\rho_1 \sim \mathrm{Normal}(0.8, 0.05)`` on
``[0,1]``) with a weaker second lag. This process is the renewal model's
reproduction-number process — it is folded into the infection model below rather
than composed separately.

```@example renewal
using ComposableTuringIDModels, Distributions, Random, Turing, Mooncake
using ADTypes: AutoMooncake
Random.seed!(1234)

latent = AR(
    damp = [truncated(Normal(0.8, 0.05), 0, 1),
        truncated(Normal(0.1, 0.05), 0, 1)],
    init = [Normal(0.0, 0.2), Normal(0.0, 0.2)],
    ϵ_t = HierarchicalNormal(std = HalfNormal(0.1)))
```

The infection process needs a discrete generation interval. [`IDData`](@ref)
takes a continuous distribution and discretises it with double interval
censoring [charniga2024best](@citep), using
[CensoredDistributions.jl](https://github.com/EpiAware/CensoredDistributions.jl).
Following [mishra2020derivation](@citet) we use a ``\mathrm{Gamma}(6.5, 0.62)``
serial interval as a proxy for the generation interval.

```@example renewal
data = IDData(gen_distribution = Gamma(6.5, 0.62))
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
[`Renewal`](@ref) is the only infection model that carries an [`IDData`](@ref),
because it is the only one that uses a generation interval.

```@example renewal
renewal = Renewal(data; rt = latent, initialisation = Normal(log(1.0), 0.1))
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
    rt = FixedIntercept(fixed_logR), initialisation = Normal())
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
obs = NegativeBinomialError(cluster_factor = HalfNormal(0.1))
nothing # hide
```

[`IDModel`](@ref) assembles the two parts — the renewal infection process
(which already carries the latent ``R_t`` process) and the observation model —
into one composed model.

```@example renewal
model = IDModel(renewal, obs)
```

Before fitting, the composed model is also a prior simulator: passing `missing`
observations makes [`as_turing_model`](@ref) return generated quantities — the
reported cases `generated_y_t`, the latent infections `I_t`, and the latent
process `Z_t` — instead of conditioning on data. That is the mechanism used for
the prior checks above; here we go straight to real data.

## The data

[mishra2020derivation](@citet) fit this model to daily test-confirmed COVID-19
cases in South Korea over the first wave of 2020. The series is stored with the
docs and read with [CSV](https://csv.juliadata.org)/[DataFrames](https://dataframes.juliadata.org).

```@example renewal
using CSV, DataFrames
datapath = joinpath(pkgdir(ComposableTuringIDModels),
    "docs", "src", "case-studies", "data", "south_korea_data.csv")
south_korea = CSV.read(datapath, DataFrame)
first(south_korea, 5)
```

We fit the growth-and-decline window of the first wave, matching the span used by
[mishra2020derivation](@citet), and take the reported cases over it as the
observed series.

```@example renewal
tspan = (45, 80)
y_obs = south_korea.cases_new[first(tspan):last(tspan)]
n = length(y_obs)
(n = n, total_cases = sum(y_obs),
    from = south_korea.date[first(tspan)], to = south_korea.date[last(tspan)])
```

## Fit

Conditioning on the observed counts and sampling with NUTS recovers the
posterior. We draw two chains in parallel with `MCMCThreads()` so the posterior
is well resolved and the cross-chain ``\hat R`` diagnostic is available; the
slightly raised target acceptance rate keeps the sampler stable on the
hierarchical innovation scale. We differentiate with
[Mooncake](https://chalk-lab.github.io/Mooncake.jl/), the recommended backend for
this package (see [Automatic differentiation backend](@ref ad-backend)).

```@example renewal
posterior = as_turing_model(model, y_obs, n)
chain = sample(
    posterior, NUTS(0.9; adtype = AutoMooncake(; config = nothing)),
    MCMCThreads(), 250, 2; progress = false)
nothing # hide
```

Sampling returns a chain whose parameters are namespaced by the component slot
that samples them, so a prior's inner variables never collide across the model.
`sample` returns a
[FlexiChains](https://github.com/penelopeysm/FlexiChains.jl) chain, which
`summarystats` summarises directly — no conversion step — giving point estimates
*and* their uncertainty alongside the effective sample size and ``\hat{R}``
convergence diagnostic. The autoregressive damping ``\rho``
(`damp_AR.θ[1]`), the innovation scale ``\sigma`` (`std`), and the
observation cluster factor ``\sqrt{1/\phi}`` (`cluster_factor`) are all
identified from the observed South Korean series:

```@example renewal
using MCMCChains
summarystats(chain)
```

## Prior versus posterior

Before reading the trajectories it is worth asking what the data taught us.
Sampling the *same* model with [`Prior`](https://turinglang.org/) — ignoring the
observations — gives a prior draw over the same parameters, and overlaying it on
the posterior shows which parameters moved. We load a
[Makie](https://docs.makie.org) backend and
[PairPlots.jl](https://sefffal.github.io/PairPlots.jl/); the FlexiChains PairPlots
extension turns a chain (subset to a few keys with `chain[[...]]`) into a
`PairPlots.Series`, so prior and posterior overlay on one corner plot.

```@example renewal
using CairoMakie, PairPlots

prior_chain = sample(posterior, Prior(), 1000; progress = false)
pp_keys = [@varname(damp_AR.θ), @varname(std),
    @varname(cluster_factor), @varname(init_incidence.θ)]
pairplot(
    PairPlots.Series(chain[pp_keys]; label = "posterior"),
    PairPlots.Series(prior_chain[pp_keys]; label = "prior"))
```

The innovation scale ``\sigma`` (`std`) is sharply updated away from
its prior — the data are informative about how much ``\log R_t`` wiggles — while
the autoregressive damping ``\rho`` (`damp_AR.θ`), the cluster factor and the
initial infections stay closer to their priors on this short window.

## Posterior trajectories

The reproduction number ``R_t = \exp(Z_t)`` is a *generated quantity* rather than
a sampled parameter. [`generated_observables`](@ref) re-runs the fitted model
over the chain to recover the latent ``Z_t`` and infection ``I_t`` trajectories
per draw. The reported counts ``y_t`` are scored element-wise, so their posterior
*predictive* distribution — fresh counts drawn under each posterior parameter set
— comes from `predict` on the same model with the observations set to `missing`.

A couple of small helpers reduce the per-draw trajectories to credible bands and
draw a median line with 50% and 95% ribbons.

```@setup renewal
using Statistics

const CI_QS = [0.025, 0.25, 0.5, 0.75, 0.975]

# time × 5 credible bands from a time × draws matrix
function credible_bands(mat; qs = CI_QS)
    reduce(hcat, (map(eachrow(mat)) do row
        vals = collect(skipmissing(row))
        isempty(vals) ? missing : quantile(vals, q)
    end for q in qs))
end

# median line with 50% and 95% ribbons
function ci_ribbon!(ax, ts, bands; color, label)
    keep = findall(!ismissing, view(bands, :, 3))
    x, b = ts[keep], Float64.(bands[keep, :])
    band!(ax, x, b[:, 1], b[:, 5]; color = (color, 0.15))
    band!(ax, x, b[:, 2], b[:, 4]; color = (color, 0.3))
    lines!(ax, x, b[:, 3]; color = color, linewidth = 2, label = label)
end

# posterior-predictive y_t bands from a `predict` chain; any leading
# indices a reporting delay leaves unscored are filled with `missing`
function predictive_bands(pred, n)
    ndraws = length(vec(pred[@varname(y_t[n])]))
    rows = map(1:n) do i
        try
            permutedims(vec(pred[@varname(y_t[i])]))
        catch
            fill(missing, 1, ndraws)
        end
    end
    credible_bands(reduce(vcat, rows))
end
```

Stack the per-draw ``Z_t`` into an ``R_t`` band, draw the posterior-predictive
``y_t`` from the unconditioned model, and plot both against the observed series:

```@example renewal
gens = vec(generated_observables(posterior, y_obs, chain).generated)
Rt = credible_bands(reduce(hcat, (exp.(g.Z_t) for g in gens)))

pred = predict(as_turing_model(model, fill(missing, n), n), chain)
yt = predictive_bands(pred, n)

fig = Figure(; size = (760, 620))
ax1 = Axis(fig[1, 1]; ylabel = "Reproduction number Rₜ")
ci_ribbon!(ax1, 1:size(Rt, 1), Rt; color = :purple, label = "posterior median")
hlines!(ax1, [1.0]; color = :grey, linestyle = :dash)
axislegend(ax1; position = :rt)
ax2 = Axis(fig[2, 1]; xlabel = "Day", ylabel = "Reported cases")
ci_ribbon!(ax2, 1:size(yt, 1), yt; color = :teal,
    label = "posterior predictive")
scatter!(ax2, 1:n, y_obs; color = :black, markersize = 7, label = "observed")
axislegend(ax2; position = :lt)
fig
```

The posterior-predictive band tracks the observed South Korean series closely,
and the ``R_t`` path recovers the first-wave turn-over: an early rise well above
one, a fall through ``R_t = 1`` as the wave peaks, and a decline below one as
cases drop.

## Swap a component

Because the parts share one interface, an alternative observation assumption is
a one-line change. Swapping the negative-binomial reporting for a
[`PoissonError`](@ref) leaves the renewal infection process — and its latent
``R_t`` process — untouched:

```@example renewal
poisson_model = IDModel(renewal, PoissonError())
length(rand(as_turing_model(poisson_model, fill(missing, n), n)))
```

## References

```@bibliography
Pages = ["renewal-negbin.md"]
Canonical = false
```
