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
weekly timescale, and assembles everything with [`EpiProblem`](@ref).

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
using EpiAwarePrototype, Distributions, Random, Turing, Mooncake
using ADTypes: AutoMooncake
Random.seed!(20240601)

arma21 = arma(
    init = [Normal(0, 0.2), Normal(0, 0.2)],
    damp = [truncated(Normal(0.1, 0.2), 0, 1), truncated(Normal(0.1, 0.05), 0, 1)],
    θ = [truncated(Normal(0.0, 0.2), -1, 1)],
    ϵ_t = HierarchicalNormal(std_prior = HalfNormal(0.1)))

arima211 = DiffLatentModel(; model = arma21, init_priors = [Normal(0.3, 0.3)])
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
data = EpiData(gen_distribution = Gamma(1.4, 1 / 0.38))
renewal = Renewal(data;
    rt = weekly_latent, initialisation_prior = Normal(log(1.0), 1.0))
nothing # hide
```

## A layered observation model

We start from the [`NegativeBinomialError`](@ref) link and build outward.
[`ascertainment_dayofweek`](@ref) wraps it with a partially pooled day-of-week
multiplier, so reporting can be systematically higher or lower on particular
weekdays.

```@example delays
negbin = NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1))
dayofweek_negbin = ascertainment_dayofweek(
    negbin; latent_model = HierarchicalNormal(std_prior = HalfNormal(1.0)))
nothing # hide
```

[`LatentDelay`](@ref) convolves the expected observations with a delay
distribution (discretised by double interval censoring). Two layers compose
sequentially: an incubation period from infection to symptom onset, then a
reporting delay from onset to report.

```@example delays
incubation = LogNormal(1.6, 0.42)   # infection -> symptom onset
reporting = LogNormal(0.58, 0.47)   # symptom onset -> report

observation = LatentDelay(LatentDelay(dayofweek_negbin, incubation), reporting)
nothing # hide
```

That single `observation` object now carries, from the inside out: a negative
binomial link, a day-of-week ascertainment modifier, an incubation-delay
convolution, and a reporting-delay convolution — assembled entirely by
composition.

## Assemble, simulate, fit

[`EpiProblem`](@ref) ties the latent, infection, and observation models to a
time span. Its [`as_turing_model`](@ref) method takes data as a named tuple with
a `y_t` field; `missing` values simulate.

```@example delays
n = 42
problem = EpiProblem(
    epi_model = renewal,
    observation_model = observation,
    tspan = (1, n))

sim = as_turing_model(problem, (y_t = fill(missing, n),))()
y_obs = sim.generated_y_t
first(y_obs, 14)
```

Fitting conditions on the simulated reports (short run for the docs build),
differentiating with the recommended [Mooncake](https://chalk-lab.github.io/Mooncake.jl/)
backend (see [Automatic differentiation backend](@ref ad-backend)):

```@example delays
chain = sample(
    as_turing_model(problem, (y_t = y_obs,)),
    NUTS(; adtype = AutoMooncake(; config = nothing)), 50; progress = false)
nothing # hide
```

```@example delays
using MCMCChains, Statistics
mc = MCMCChains.Chains(chain)
summarystats(mc[[:cluster_factor, Symbol("DayofWeek.std")]])
```

`DayofWeek.std` is the scale of the partially pooled weekday multipliers (its
own block, prefixed because the ascertainment modifier introduces a named
sub-process); `cluster_factor` is the negative-binomial overdispersion. The
day-of-week effect, the two delay kernels, and the weekly reproduction number
were all estimated jointly — and any of them can be swapped or removed by
editing one line of the composition above.

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
