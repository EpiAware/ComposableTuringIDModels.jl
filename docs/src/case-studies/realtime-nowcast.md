# [Real-time nowcasting: correcting right-truncation](@id case-study-nowcast)

In real time, the most recent days of a surveillance series are **incomplete**:
a case with reference day ``t`` is only reported after a reporting delay, so on
day ``\text{now}`` we have seen just a fraction of the cases that will eventually
be attributed to recent days. Fitting a renewal model to this *right-truncated*
tail without correction produces the familiar artefact of real-time estimation —
an apparent late down-turn in ``R_t`` that is really an artefact of not-yet-reported
counts [abbott2020estimating](@citep).

This case study takes a real, fully-reported case series — the daily confirmed
COVID-19 cases from Italy's first wave used in the
[delays example](@ref case-study-delays) — treats it as the *eventual* totals,
and truncates its recent tail to mimic the real-time snapshot an analyst would
have seen mid-outbreak. It then contrasts two fits: a **naive** one that treats
the truncated counts as complete, and one that wraps the observation model in
[`RightTruncate`](@ref) to scale each reference day's expected count by the
fraction of its eventual total reported so far. The correction is the
EpiNow2-style CDF-scaling nowcast [abbott2020estimating](@citep), expressed here
as a one-line observation modifier.

## The idea

The infection pipeline produces ``Y_t = \mu_t``, the expected *eventual* total
for reference day ``t``. At time ``\text{now}`` a reference day of age
``a = \text{now} - t`` has only had a fraction ``F[a+1]`` of its eventual total
reported, where ``F`` is the reporting-delay CDF. The expected observed-so-far is
therefore ``\mu_t \cdot F[a+1]``. [`RightTruncate`](@ref) conditions the
observation error on that down-weighted mean, so the model's ``Y_t`` stays the
eventual total and the nowcast is just ``Y_t`` read back out.

```math
\begin{aligned}
R_t &= \exp(Z_t), & Z_t &= \rho\, Z_{t-1} + \epsilon_t, \\
I_t &= R_t \sum_{s \ge 1} g_s\, I_{t-s}, \\
y_t &\sim \mathrm{NegBinomial}\!\big(I_t \cdot F[(\text{now}-t)+1],\ \phi\big).
\end{aligned}
```

A naive fit drops the ``F[\cdot]`` factor (equivalently assumes ``F \equiv 1``).

## The full-data model

We build the same composed renewal model as the
[renewal case study](@ref case-study-renewal): an autoregressive ``\log R_t``
process folded into a [`Renewal`](@ref) infection process, observed with a
[`NegativeBinomialError`](@ref).

```@example nowcast
using ComposableTuringIDModels, Distributions, Random, Turing
using CSV, DataFrames
Random.seed!(20240625)

latent = AR(
    damp = [truncated(Normal(0.8, 0.05), 0, 1)],
    init = [Normal(0.0, 0.25)],
    ϵ_t = HierarchicalNormal(std = HalfNormal(0.1)))
data = IDData(gen_distribution = Gamma(1.4, 1 / 0.38))
renewal = Renewal(data; rt = latent, initialisation = Normal(log(1.0), 1.0))
error = NegativeBinomialError(cluster_factor = HalfNormal(0.1))
nothing # hide
```

## Take a real series and truncate it

We use the fully-reported Italy confirmed-case series as the *eventual* totals,
then impose a reporting delay and truncate at ``\text{now}`` to reconstruct the
partially-reported tail a real-time analyst would have seen. The reporting delay
is a continuous distribution discretised to a CDF with the same released
[CensoredDistributions.jl](https://github.com/EpiAware/CensoredDistributions.jl)
path the rest of the package uses; [`ReportingCDF`](@ref) builds the completeness
curve ``F``.

```@example nowcast
datapath = joinpath(pkgdir(ComposableTuringIDModels),
    "docs", "src", "case-studies", "data", "italy_data.csv")
italy = CSV.read(datapath, DataFrame)
n = 50
eventual = italy.confirm[1:n]                # the eventual (complete) totals

reporting_delay = LogNormal(1.6, 0.5)        # mean ≈ 5.6 days
cdf_curve = ReportingCDF(reporting_delay)
nothing # hide
```

Now form the right-truncated snapshot: thin each reference day's eventual total
to the fraction reported by ``\text{now} = n``. The completeness by reference day
is ``F`` reversed onto the time axis (the most recent day is least complete).

```@example nowcast
completeness = as_turing_model(cdf_curve, n)()    # F by age, F[1] = age 0
scale = reverse(completeness)                     # by reference day t = 1..n
observed_so_far = @. rand(Binomial(eventual, scale))
(complete_tail = eventual[(end - 4):end], truncated_tail = observed_so_far[(end - 4):end])
```

The last few days are visibly thinned: the most recent day shows only a fraction
of its eventual count.

## Fit 1 — naive (no truncation correction)

Condition the plain renewal model on the truncated counts as though they were
complete.

```@example nowcast
naive_model = IDModel(renewal, error)
naive_post = as_turing_model(naive_model, observed_so_far, n)
naive_chain = sample(
    naive_post, NUTS(0.9), MCMCThreads(), 300, 2; progress = false)
nothing # hide
```

## Fit 2 — corrected with [`RightTruncate`](@ref)

Wrap the *same* error model in [`RightTruncate`](@ref) with the *same* reporting
delay. Nothing else changes — the infection process and its latent ``R_t`` are
untouched; only how the expected counts are compared to the truncated data.

```@example nowcast
corrected_obs = RightTruncate(error, reporting_delay)
corrected_model = IDModel(renewal, corrected_obs)
corrected_post = as_turing_model(corrected_model, observed_so_far, n)
corrected_chain = sample(
    corrected_post, NUTS(0.9), MCMCThreads(), 300, 2; progress = false)
nothing # hide
```

## The right-truncation bias and its correction

The reproduction number ``R_t = \exp(Z_t)`` is a generated quantity of the model,
not a sampled parameter. Re-running the fitted model over the chain with
[`returned`](https://turinglang.org/Turing.jl/) (re-exported by `Turing`) recovers
the latent ``Z_t`` trajectory per draw, from which we take the posterior-mean
``R_t`` and average it over the most recent window. As a reference we also fit the
plain model to the **complete** (untruncated) series — what the analyst would
eventually see. Right-truncation biases the naive fit's recent ``R_t`` *downward*
(the not-yet-reported counts look like a decline); the [`RightTruncate`](@ref) fit
removes that bias.

```@example nowcast
using Statistics

function recent_Rt(posterior, chain; window = 7)
    # Per-draw R_t = exp(Z_t); average over draws, then over the recent window.
    gens = vec(returned(posterior, chain))
    Rt_mean = mean(exp.(g.Z_t) for g in gens)
    mean(Rt_mean[(end - window + 1):end])
end

complete_post = as_turing_model(naive_model, eventual, n)
complete_chain = sample(
    complete_post, NUTS(0.9), MCMCThreads(), 300, 2; progress = false)

R_complete_recent = recent_Rt(complete_post, complete_chain)
R_naive_recent = recent_Rt(naive_post, naive_chain)
R_corrected_recent = recent_Rt(corrected_post, corrected_chain)

(complete = round(R_complete_recent, digits = 2),
    naive = round(R_naive_recent, digits = 2),
    corrected = round(R_corrected_recent, digits = 2))
```

The naive recent-``R_t`` sits **below** the complete-data estimate — the
artefactual late down-turn produced by treating the not-yet-reported tail as
complete — whereas the [`RightTruncate`](@ref) fit, which knows the recent days
are incomplete, pulls the recent ``R_t`` back **up** off that spurious decline
towards the complete-data value. (There is still Monte Carlo noise in the exact
values, but the robust, repeatable signal is the *direction* — the naive fit
under-estimates recent ``R_t``, and the correction removes that downward pull.)
The nowcast of the eventual totals is the corrected model's ``Y_t``, recovered the
same way from the per-draw generated quantities; the figures below make the
correction visible.

## Prior versus posterior

Sampling the corrected model with [`Prior`](https://turinglang.org/) gives a prior
draw over the shared renewal parameters. Overlaying it on the posterior with
[PairPlots.jl](https://sefffal.github.io/PairPlots.jl/) confirms the truncation
correction still identifies them from the thinned tail.

```@example nowcast
using CairoMakie, PairPlots

prior_chain = sample(corrected_post, Prior(), 1000; progress = false)
pp_keys = [@varname(Z_t.damp_AR.θ), @varname(Z_t.ϵ_t.std.θ),
    @varname(cluster_factor.θ), @varname(init_incidence.θ)]
pairplot(
    PairPlots.Series(corrected_chain[pp_keys]; label = "posterior"),
    PairPlots.Series(prior_chain[pp_keys]; label = "prior"))
```

## The correction in a figure

``R_t = \exp(Z_t)`` is a generated quantity recovered per draw with
[`generated_observables`](@ref). Plotting the posterior median and 95% band for
all three fits over time — the complete-data reference, the naive truncated fit,
and the [`RightTruncate`](@ref)-corrected fit — shows the right-truncation
artefact and its removal in the recent window (shaded).

```@setup nowcast
using Statistics

const CI_QS = [0.025, 0.25, 0.5, 0.75, 0.975]

function credible_bands(mat; qs = CI_QS)
    reduce(hcat, (map(eachrow(mat)) do row
        vals = collect(skipmissing(row))
        isempty(vals) ? missing : quantile(vals, q)
    end for q in qs))
end

# posterior R_t = exp(Z_t) credible bands for one fit
function rt_bands(post, data, chn)
    gens = vec(generated_observables(post, data, chn).generated)
    credible_bands(reduce(hcat, (exp.(g.Z_t) for g in gens)))
end

# median line with a 95% ribbon
function rt_line!(ax, ts, bands; color, label)
    band!(ax, ts, Float64.(bands[:, 1]), Float64.(bands[:, 5]);
        color = (color, 0.12))
    lines!(ax, ts, Float64.(bands[:, 3]); color = color, linewidth = 2,
        label = label)
end
```

```@example nowcast
Rt_complete = rt_bands(complete_post, eventual, complete_chain)
Rt_naive = rt_bands(naive_post, observed_so_far, naive_chain)
Rt_corrected = rt_bands(corrected_post, observed_so_far, corrected_chain)

ts = 1:n
fig = Figure(; size = (760, 420))
ax = Axis(fig[1, 1]; xlabel = "Reference day",
    ylabel = "Reproduction number Rₜ")
vspan!(ax, n - 6, n; color = (:grey, 0.15))
rt_line!(ax, ts, Rt_complete; color = :black, label = "complete (reference)")
rt_line!(ax, ts, Rt_naive; color = :crimson, label = "naive (truncated)")
rt_line!(ax, ts, Rt_corrected; color = :seagreen, label = "corrected")
hlines!(ax, [1.0]; color = :grey, linestyle = :dash)
axislegend(ax; position = :lb)
fig
```

In the shaded recent window the naive fit (red) dips below the complete-data
reference (black) — the artefactual late down-turn — while the
[`RightTruncate`](@ref) correction (green) pulls the recent ``R_t`` back up
towards the reference, having accounted for the not-yet-reported counts.

## Reading the shared parameters

Wrapping the error model in [`RightTruncate`](@ref) does not touch the renewal
process, so the corrected fit recovers the *same* shared parameters: the
autoregressive damping ``\rho`` (`Z_t.damp_AR.θ[1]`), the innovation scale
``\sigma`` (`Z_t.ϵ_t.std.θ`), the observation overdispersion
(`cluster_factor.θ`), and the initial infections (`init_incidence.θ`). Each is
namespaced by the component slot that samples it, so a prior's inner variables
never collide across the model.

`sample` returns a [FlexiChains](https://github.com/penelopeysm/FlexiChains.jl)
chain, which `summarystats` summarises directly — no conversion step — giving
point estimates *and* their uncertainty alongside the effective sample size and
``\hat{R}`` convergence diagnostic:

```@example nowcast
using MCMCChains
summarystats(corrected_chain)
```

## From the marginal to the joint

[`RightTruncate`](@ref) corrects right-truncation by conditioning on each
reference day's observed-so-far *total* — the **marginal** of the full
reference-day × reporting-delay structure. When the delay structure itself is of
interest (e.g. reporting that drifts over the outbreak), the
[`ReportTriangle`](@ref) observation model keeps the full 2D reporting triangle
and scores it cell by cell; its observed row-sums reconcile with this CDF-scaling
to machine precision. The marginal correction here is the cheaper, released-code
first step.

## References

```@bibliography
Pages = ["realtime-nowcast.md"]
Canonical = false
```
