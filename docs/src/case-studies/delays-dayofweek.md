# [Reporting delays and day-of-week effects](@id case-study-delays)

Real surveillance data is rarely a clean count of infections on the day they
occur. Cases are reported after a delay — an incubation period followed by a
reporting lag — and the number reported depends on the day of the week. Tools
for real-time estimation such as those of [abbott2020estimating](@citet) build
these features into the observation model so that the latent infection signal is
estimated free of reporting artefacts.

This case study keeps the renewal infection core of the
[previous example](@ref case-study-renewal) but replaces the simple observation
model with a layered one: infections are convolved through two delay
distributions and then modulated by a day-of-week reporting pattern. It also
shows the latent process as an ARIMA-style differenced process broadcast to a
weekly timescale, and assembles everything with [`IDProblem`](@ref). The model
follows the configuration of the [EpiNow2](https://epiforecasts.io/EpiNow2/)
package [abbott2020estimating](@citep) and is fit to daily confirmed COVID-19
cases from Italy's first wave in 2020.

## The model

```math
\begin{aligned}
Z_w &\sim \text{ARIMA}(2,1,1), & R_t &= \exp\!\big(Z_{\lfloor t/7 \rfloor}\big), \\
I_t &= R_t \sum_{s \ge 1} g_s\, I_{t-s}, \\
S_t &= \sum_{s \ge 1} \eta_s\, I_{t-s}, & D_t &= \sum_{s \ge 1} \xi_s\, S_{t-s}, \\
y_t &\sim \mathrm{NegBinomial}\big(\omega_{t \bmod 7}\, D_t, \phi\big),
\end{aligned}
```

where ``\eta`` is the incubation-period pmf, ``\xi`` the reporting-delay pmf,
and ``\omega`` a day-of-week reporting multiplier.

## A weekly latent process

The latent process is an ARIMA(2,1,1): an [`AR`](@ref)/[`MA`](@ref) combination
([`arma`](@ref)) wrapped in a [`DiffLatentModel`](@ref) to difference it once.
Differencing makes the level a random walk rather than mean-reverting, which
suits a reproduction number that can drift.

```@example delays
using ComposableTuringIDModels, Distributions, Random, Turing, Mooncake
using ADTypes: AutoMooncake
Random.seed!(20240601)

arma21 = arma(
    init = [Normal(0, 0.2), Normal(0, 0.2)],
    damp = [truncated(Normal(0.1, 0.2), 0, 1), truncated(Normal(0.1, 0.05), 0, 1)],
    θ = [truncated(Normal(0.0, 0.2), -1, 1)],
    ϵ_t = HierarchicalNormal(std = HalfNormal(0.1)))

arima211 = DiffLatentModel(; model = arma21, init = [Normal(0.3, 0.3)])
nothing # hide
```

[`broadcast_weekly`](@ref) makes the process piecewise-constant by week: a new
value is drawn each week and held for seven days. This models ``R_t`` as
changing weekly rather than daily, which both regularises the estimate and cuts
the number of latent parameters.

```@example delays
weekly_latent = broadcast_weekly(arima211)
nothing # hide
```

## The infection process

As before, a [`Renewal`](@ref) process driven by a discretised generation
interval. Here we use a ``\mathrm{Gamma}(1.4, 1/0.38)`` generation time. The
weekly ``\log R_t`` process built above is folded into the renewal model's `rt`
slot.

```@example delays
data = IDData(gen_distribution = Gamma(1.4, 1 / 0.38))
renewal = Renewal(data;
    rt = weekly_latent, initialisation = Normal(log(1.0), 1.0))
nothing # hide
```

## A layered observation model

We start from the [`NegativeBinomialError`](@ref) link and build outward.
[`ascertainment_dayofweek`](@ref) wraps it with a partially pooled day-of-week
multiplier, so reporting can be systematically higher or lower on particular
weekdays.

```@example delays
negbin = NegativeBinomialError(cluster_factor = HalfNormal(0.1))
dayofweek_negbin = ascertainment_dayofweek(
    negbin; latent_model = HierarchicalNormal(std = HalfNormal(1.0)))
nothing # hide
```

[`LatentDelay`](@ref) convolves the expected observations with a delay
distribution (discretised by double interval censoring). Two layers compose
sequentially: a fixed incubation period from infection to symptom onset, then a
reporting delay from onset to report whose parameters are *inferred*. The
reporting delay is an [`UncertainDelay`](@ref): its `LogNormal` log-scale mean and
standard deviation carry priors, so the delay is rediscretised each draw and
estimated jointly with the reproduction number rather than fixed from external
data.

```@example delays
incubation = LogNormal(1.6, 0.42)   # infection -> symptom onset (fixed)
reporting = UncertainDelay(         # symptom onset -> report (inferred)
    LogNormal, [Normal(0.58, 0.3), truncated(Normal(0.47, 0.2), 0, Inf)];
    D = 8.0)

observation = LatentDelay(LatentDelay(dayofweek_negbin, incubation), reporting)
nothing # hide
```

That single `observation` object now carries, from the inside out: a negative
binomial link, a day-of-week ascertainment modifier, a fixed incubation-delay
convolution, and an inferred reporting-delay convolution — assembled entirely by
composition. The reporting-delay parameters flow through the same priors seam as
every other parameter, so inferring the delay needs no change to the rest of the
model.

## The data

We fit the model to the daily confirmed COVID-19 cases from Italy's first wave
(the example series shipped with the EpiNow2 package), stored with the docs.

```@example delays
using CSV, DataFrames
datapath = joinpath(pkgdir(ComposableTuringIDModels),
    "docs", "src", "case-studies", "data", "italy_data.csv")
italy = CSV.read(datapath, DataFrame)
n = 42
y_obs = italy.confirm[1:n]
(n = n, total_cases = sum(y_obs), from = italy.date[1], to = italy.date[n])
```

## Assemble and fit

[`IDProblem`](@ref) ties the latent, infection, and observation models to a
time span. Its [`as_turing_model`](@ref) method takes data as a named tuple with
a `y_t` field (passing `missing` values would instead simulate from the prior).

```@example delays
problem = IDProblem(
    infection = renewal,
    observation_model = observation,
    tspan = (1, n))
nothing # hide
```

Fitting conditions on the observed reports, differentiating with the recommended
[Mooncake](https://chalk-lab.github.io/Mooncake.jl/) backend (see
[Automatic differentiation backend](@ref ad-backend)). We draw two chains in
parallel with `MCMCThreads()`, which gives a cross-chain ``\hat R``:

```@example delays
posterior = as_turing_model(problem, (y_t = y_obs,))
chain = sample(
    posterior, NUTS(0.9; adtype = AutoMooncake(; config = nothing)),
    MCMCThreads(), 250, 2; progress = false)
nothing # hide
```

`sample` returns a [FlexiChains](https://github.com/penelopeysm/FlexiChains.jl)
chain, which `summarystats` summarises directly — no conversion step. The
day-of-week scale (`DayofWeek.std`), the negative-binomial overdispersion
(`cluster_factor`) and the inferred reporting-delay parameters (`delay.θ`, the
`LogNormal` log-mean and log-sd) appear alongside the latent-process parameters:

```@example delays
using MCMCChains
summarystats(chain)
```

`DayofWeek.std` is the scale of the partially
pooled weekday multipliers (its own block, namespaced because the ascertainment
modifier introduces a named sub-process); `cluster_factor` is the
negative-binomial overdispersion; `delay.θ` are the inferred reporting-delay
parameters. The day-of-week effect, the reporting delay, and the weekly
reproduction number were all estimated jointly — and any of them can be swapped,
fixed, or removed by editing one line of the composition above.

## Prior versus posterior

Sampling the same model with [`Prior`](https://turinglang.org/) gives a prior
draw over the same parameters. Overlaying it on the posterior with
[PairPlots.jl](https://sefffal.github.io/PairPlots.jl/) — the FlexiChains
extension turns each chain, subset to a few keys, into a `PairPlots.Series` —
shows which parameters the six weeks of Italian data moved.

```@example delays
using CairoMakie, PairPlots

prior_chain = sample(posterior, Prior(), 1000; progress = false)
pp_keys = [@varname(damp_AR.θ), @varname(θ.θ),
    @varname(std), @varname(cluster_factor)]
pairplot(
    PairPlots.Series(chain[pp_keys]; label = "posterior"),
    PairPlots.Series(prior_chain[pp_keys]; label = "prior"))
```

The innovation scale ``\sigma`` (`std`)
and the negative-binomial overdispersion (`cluster_factor`) tighten under the
data, while the autoregressive damping
(`damp_AR.θ`) and moving-average
(`θ.θ`) coefficients of the ARIMA process stay
close to their weakly informative priors.

## Posterior trajectories

``R_t = \exp(Z_t)`` and the infections ``I_t`` are generated quantities recovered
per draw with [`generated_observables`](@ref); the reports ``y_t`` are scored
element-wise, so their posterior-predictive distribution comes from `predict` on
the model with the observations set to `missing`. Two small helpers reduce the
per-draw trajectories to credible bands.

```@setup delays
using Statistics

const CI_QS = [0.025, 0.25, 0.5, 0.75, 0.975]

function credible_bands(mat; qs = CI_QS)
    reduce(hcat, (map(eachrow(mat)) do row
        vals = collect(skipmissing(row))
        isempty(vals) ? missing : quantile(vals, q)
    end for q in qs))
end

function ci_ribbon!(ax, ts, bands; color, label)
    keep = findall(!ismissing, view(bands, :, 3))
    x, b = ts[keep], Float64.(bands[keep, :])
    band!(ax, x, b[:, 1], b[:, 5]; color = (color, 0.15))
    band!(ax, x, b[:, 2], b[:, 4]; color = (color, 0.3))
    lines!(ax, x, b[:, 3]; color = color, linewidth = 2, label = label)
end

# the two delay convolutions leave the first few reference days unscored,
# so those predictive entries are filled with `missing` and skipped
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

```@example delays
gens = vec(generated_observables(posterior, (y_t = y_obs,), chain).generated)
Rt = credible_bands(reduce(hcat, (exp.(g.Z_t) for g in gens)))

pred = predict(as_turing_model(problem, (y_t = fill(missing, n),)), chain)
yt = predictive_bands(pred, n)

fig = Figure(; size = (760, 620))
ax1 = Axis(fig[1, 1]; ylabel = "Reproduction number Rₜ")
ci_ribbon!(ax1, 1:size(Rt, 1), Rt; color = :purple, label = "posterior median")
hlines!(ax1, [1.0]; color = :grey, linestyle = :dash)
axislegend(ax1; position = :rt)
ax2 = Axis(fig[2, 1]; xlabel = "Day", ylabel = "Confirmed cases")
ci_ribbon!(ax2, 1:size(yt, 1), yt; color = :teal,
    label = "posterior predictive")
scatter!(ax2, 1:n, y_obs; color = :black, markersize = 7, label = "observed")
axislegend(ax2; position = :lt)
fig
```

The weekly ``R_t`` is piecewise-constant by construction, stepping down through
one as the first wave turns over. The posterior-predictive band starts partway
into the series — the two delay convolutions leave the earliest reference days
without a fully supported expected count — and from there tracks the observed
Italian reports, the layered observation model having absorbed the reporting
pattern rather than the infection signal.

## A time-varying reporting pattern

The day-of-week multiplier above is *static*: one weekly profile held fixed
across the series. Reporting behaviour can itself drift — testing capacity
changes, weekend effects strengthen or weaken — and the same composition
expresses that. Because the ascertainment modifier takes any latent model,
replacing the pooled [`HierarchicalNormal`](@ref) weekday effect with a
[`BroadcastLatentModel`](@ref) over a process that evolves week to week turns the
fixed profile into a time-varying one, at the cost of more latent parameters. The
structural change is again local to the observation model; the infection and
latent ``R_t`` parts are untouched. We keep the static pattern here — it is
identifiable from six weeks of data, where a fully time-varying weekday process
would not be — and flag the richer variant rather than fit it.

## References

```@bibliography
Pages = ["delays-dayofweek.md"]
Canonical = false
```
