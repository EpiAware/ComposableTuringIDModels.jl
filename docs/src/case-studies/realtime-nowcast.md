# [Real-time nowcasting: right-truncation vs the reporting triangle](@id case-study-nowcast)

In real time the most recent days of a surveillance series are **incomplete**.
A case with reference day ``t`` is only counted after a reporting delay, so on
day ``\text{now}`` we have seen a fraction of the cases that will eventually be
attributed to recent days.
Fitting a renewal model to this *right-truncated* tail without correction
produces the familiar artefact of real-time estimation, an apparent late
down-turn in ``R_t`` that is really the not-yet-reported counts
[abbott2020estimating](@citep).

This package corrects right-truncation two ways, and this case study fits and
compares both on the same simulated data.
[`RightTruncate`](@ref) is the **marginal** correction.
It keeps only each reference day's observed-so-far *total* and scales the
expected count by the reporting-delay CDF, the EpiNow2-style CDF-scaling nowcast.
[`ReportTriangle`](@ref) is the **joint** correction.
It keeps the full reference-day × reporting-delay triangle and scores it cell by
cell, the epinowcast-style nowcast.
The two share one reporting-delay kernel, so the triangle's observed row-sums are
exactly the observed-so-far totals the marginal correction conditions on.
Both leave the renewal infection process and its latent ``R_t`` untouched, and
both keep the model's ``Y_t`` as the eventual total, so the nowcast is just
``Y_t`` read back out.

## The two corrections

The infection pipeline produces ``Y_t = \mu_t``, the expected *eventual* total for
reference day ``t``.
Write ``p`` for the reporting-delay PMF (delay ``d = 0, 1, \dots``) and
``F[a+1] = \sum_{d=0}^{a} p[d+1]`` for its CDF.
A reference day of age ``a = \text{now} - t`` has reported a fraction ``F[a+1]``
of its eventual total.

```math
\begin{aligned}
R_t &= \exp(Z_t), & Z_t &= \rho\, Z_{t-1} + \epsilon_t, \\
I_t &= R_t \sum_{s \ge 1} g_s\, I_{t-s}, & Y_t &= I_t, \\[2pt]
\text{marginal:} \quad y_t &\sim \mathrm{NegBinomial}\!\big(Y_t\, F[(\text{now}-t)+1],\ \phi\big), \\
\text{joint:} \quad N_{t,d} &\sim \mathrm{NegBinomial}\!\big(Y_t\, p[d+1],\ \phi\big), \quad t + d \le \text{now}.
\end{aligned}
```

Summing the joint model's observed cells over ``d`` recovers the marginal mean
``Y_t\, F[(\text{now}-t)+1]``, so the two are consistent by construction.
A naive fit drops the delay factor (equivalently ``F \equiv 1``).

## The renewal model

Both corrections wrap the same composed renewal model as the
[renewal case study](@ref case-study-renewal), an autoregressive ``\log R_t``
process folded into a [`Renewal`](@ref) infection process observed with a
[`NegativeBinomialError`](@ref).

```@example nowcast
using ComposableTuringIDModels, Distributions, Random, Turing, Mooncake
using ADTypes: AutoMooncake
using CSV, DataFrames
Random.seed!(20240625)

adt = AutoMooncake(; config = nothing)

latent = AR(
    damp = [truncated(Normal(0.8, 0.05), 0, 1)],
    init = [Normal(0.0, 0.25)],
    ϵ_t = HierarchicalNormal(std = HalfNormal(0.1)))
renewal = Renewal(; generation_time = Gamma(1.4, 1 / 0.38),
    rt = latent, initialisation = Normal(log(1.0), 1.0))
error = NegativeBinomialError(cluster_factor = HalfNormal(0.1))
nothing # hide
```

## Simulate a reporting triangle

We take the fully-reported Italy confirmed-case series as the *eventual* totals,
split each day's total across reporting delays with the delay PMF, then mask the
cells not yet reported at ``\text{now} = n``.
[`ReportingPMF`](@ref) discretises the reporting delay with the same
[CensoredDistributions.jl](https://github.com/EpiAware/CensoredDistributions.jl)
path the rest of the package uses.

```@example nowcast
datapath = joinpath(pkgdir(ComposableTuringIDModels),
    "docs", "src", "case-studies", "data", "italy_data.csv")
italy = CSV.read(datapath, DataFrame)
n = 45
eventual = italy.confirm[1:n]                 # eventual (complete) totals

reporting_delay = LogNormal(1.6, 0.5)         # mean ≈ 5.6 days
delay_pmf = ReportingPMF(reporting_delay; D = 10)
pmf = delay_pmf.pmf                            # delays 0 … Dmax
Dmax = length(pmf) - 1
nothing # hide
```

Split each reference day's eventual total across delays with a multinomial draw,
then keep only the cells reported by ``\text{now} = n`` (``t + d \le n``).
The observed-so-far total per reference day is the row-sum of the reported cells,
the marginal the CDF-scaling correction sees.

```@example nowcast
full_triangle = reduce(vcat,
    (permutedims(rand(Multinomial(eventual[t], pmf))) for t in 1:n))
mask = [t + d <= n for t in 1:n, d in 0:Dmax]
reported_triangle = full_triangle .* mask
observed_so_far = vec(sum(reported_triangle, dims = 2))
(complete_tail = eventual[(end - 4):end],
    truncated_tail = observed_so_far[(end - 4):end])
```

The recent days are visibly thinned, the most recent day showing only a fraction
of its eventual count.

## Visualise the data

The reporting triangle is the native data object.
Plotting it with [AlgebraOfGraphics](https://aog.makie.org), the not-yet-reported
cells blanked, shows the staircase of missing counts in the recent corner.

```@example nowcast
<<<<<<< ours
naive_model = IDModel(renewal, error)
naive_post = as_turing_model(naive_model, observed_so_far, n)
naive_chain = sample(
    naive_post, NUTS(0.9), MCMCThreads(), 250, 2; progress = false)
nothing # hide
=======
using AlgebraOfGraphics, CairoMakie, DataFrames

tri_df = DataFrame(
    reference = repeat(1:n, outer = Dmax + 1),
    delay = repeat(0:Dmax, inner = n),
    count = vec(full_triangle),
    reported = vec(mask))
tri_df.shown = ifelse.(tri_df.reported, Float64.(tri_df.count), NaN)

draw(data(tri_df) *
     mapping(:reference, :delay, :shown => "Reported cases") * visual(Heatmap);
    axis = (xlabel = "Reference day", ylabel = "Reporting delay (days)"))
>>>>>>> theirs
```

The same truncation reads as a shortfall in the tail when the observed-so-far
row-sums are drawn against the eventual totals.

```@example nowcast
<<<<<<< ours
corrected_obs = RightTruncate(error, reporting_delay)
corrected_model = IDModel(renewal, corrected_obs)
corrected_post = as_turing_model(corrected_model, observed_so_far, n)
corrected_chain = sample(
    corrected_post, NUTS(0.9), MCMCThreads(), 250, 2; progress = false)
nothing # hide
=======
comp_df = DataFrame(
    reference = repeat(1:n, 2),
    count = vcat(eventual, observed_so_far),
    series = repeat(["eventual total", "observed so far"], inner = n))

draw(data(comp_df) * mapping(:reference, :count, color = :series) *
     visual(Lines); axis = (xlabel = "Reference day", ylabel = "Cases"))
>>>>>>> theirs
```

## Three fits

We fit the plain renewal model three ways with NUTS, drawing **two chains** with
**1000 warmup** iterations each so the cross-chain ``\hat R`` diagnostic is
available, differentiating with
[Mooncake](https://chalk-lab.github.io/Mooncake.jl/) (see
[Automatic differentiation backend](@ref ad-backend)).
The naive fit treats the truncated totals as complete, and motivates the problem.
[`RightTruncate`](@ref) applies the marginal correction to the same observed-so-far
totals.
[`ReportTriangle`](@ref) applies the joint correction to the reporting triangle,
built through the shared [`define_y_t`](@ref) hook.

```@example nowcast
naive_model = IDModel(renewal, error)
naive_post = as_turing_model(naive_model, observed_so_far, n)
naive_chain = sample(
    naive_post, NUTS(1000, 0.9; adtype = adt), MCMCThreads(), 250, 2;
    progress = false)

rt_obs = RightTruncate(error, ReportingCDF(reporting_delay; D = 10))
rt_model = IDModel(renewal, rt_obs)
rt_post = as_turing_model(rt_model, observed_so_far, n)
rt_chain = sample(
    rt_post, NUTS(1000, 0.9; adtype = adt), MCMCThreads(), 250, 2;
    progress = false)

tri_obs = ReportTriangle(error, delay_pmf)
tri_data = define_y_t(tri_obs, reported_triangle, eventual)
tri_model = IDModel(renewal, tri_obs)
tri_post = as_turing_model(tri_model, tri_data, n)
tri_chain = sample(
    tri_post, NUTS(1000, 0.9; adtype = adt), MCMCThreads(), 250, 2;
    progress = false)
nothing # hide
```

As a reference we also fit the plain model to the **complete** (untruncated)
series, what the analyst would eventually see.

```@example nowcast
complete_post = as_turing_model(naive_model, eventual, n)
complete_chain = sample(
<<<<<<< ours
    complete_post, NUTS(0.9), MCMCThreads(), 250, 2; progress = false)

R_complete_recent = recent_Rt(complete_post, complete_chain)
R_naive_recent = recent_Rt(naive_post, naive_chain)
R_corrected_recent = recent_Rt(corrected_post, corrected_chain)

(complete = round(R_complete_recent, digits = 2),
    naive = round(R_naive_recent, digits = 2),
    corrected = round(R_corrected_recent, digits = 2))
=======
    complete_post, NUTS(1000, 0.9; adtype = adt), MCMCThreads(), 250, 2;
    progress = false)
nothing # hide
>>>>>>> theirs
```

## Recent Rt

``R_t = \exp(Z_t)`` is a generated quantity, recovered per draw by re-running the
fitted model over the chain with [`generated_observables`](@ref).
Averaging it over the most recent window shows the bias and its removal.
Right-truncation biases the naive recent ``R_t`` *downward*, so the naive value
sits below the two corrections, which lift it back up towards the complete-data
reference.
The exact numbers carry Monte Carlo noise on this short run, but the direction is
the robust, repeatable signal.

```@example nowcast
using Statistics

function recent_Rt(post, chain; window = 7)
    gens = vec(generated_observables(post, nothing, chain).generated)
    Rt_mean = mean(exp.(g.Z_t) for g in gens)
    round(mean(Rt_mean[(end - window + 1):end]), digits = 2)
end

(complete = recent_Rt(complete_post, complete_chain),
    naive = recent_Rt(naive_post, naive_chain),
    right_truncate = recent_Rt(rt_post, rt_chain),
    report_triangle = recent_Rt(tri_post, tri_chain))
```

## Posterior prediction, nowcast, and Rt

Three views of the fits share the reference-day axis with the recent window
shaded.
The **``R_t``** panel plots the ``\exp(Z_t)`` bands of all four fits.
The **nowcast** panel plots the reconstructed eventual totals ``Y_t = I_t`` for
the two corrections against the true eventual totals and the truncated
observed-so-far.
The **posterior prediction** panel plots the [`RightTruncate`](@ref) model's
posterior-predictive observed-so-far, recovered with `predict` on the same model
with the observations set to `missing`, against the data it was fit to.

```@setup nowcast
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

# generated-quantity bands: stack f(draw) over the chain
function gen_bands(post, chain, f)
    gens = vec(generated_observables(post, nothing, chain).generated)
    credible_bands(reduce(hcat, (f(g) for g in gens)))
end

# posterior-predictive y_t bands from a `predict` chain
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

```@example nowcast
ts = 1:n
Rt_complete = gen_bands(complete_post, complete_chain, g -> exp.(g.Z_t))
Rt_naive = gen_bands(naive_post, naive_chain, g -> exp.(g.Z_t))
Rt_rt = gen_bands(rt_post, rt_chain, g -> exp.(g.Z_t))
Rt_tri = gen_bands(tri_post, tri_chain, g -> exp.(g.Z_t))

now_rt = gen_bands(rt_post, rt_chain, g -> g.I_t)
now_tri = gen_bands(tri_post, tri_chain, g -> g.I_t)

rt_pred = predict(as_turing_model(rt_model, fill(missing, n), n), rt_chain)
yt = predictive_bands(rt_pred, n)

fig = Figure(; size = (780, 820))
ax1 = Axis(fig[1, 1]; ylabel = "Reproduction number Rₜ")
vspan!(ax1, n - 6, n; color = (:grey, 0.15))
ci_ribbon!(ax1, ts, Rt_complete; color = :black, label = "complete (reference)")
ci_ribbon!(ax1, ts, Rt_naive; color = :crimson, label = "naive")
ci_ribbon!(ax1, ts, Rt_rt; color = :seagreen, label = "right-truncate")
ci_ribbon!(ax1, ts, Rt_tri; color = :steelblue, label = "report-triangle")
hlines!(ax1, [1.0]; color = :grey, linestyle = :dash)
axislegend(ax1; position = :lb, nbanks = 2)

ax2 = Axis(fig[2, 1]; ylabel = "Eventual cases (nowcast)")
vspan!(ax2, n - 6, n; color = (:grey, 0.15))
ci_ribbon!(ax2, ts, now_rt; color = :seagreen, label = "right-truncate")
ci_ribbon!(ax2, ts, now_tri; color = :steelblue, label = "report-triangle")
lines!(ax2, ts, Float64.(eventual); color = :black, linewidth = 2,
    label = "true eventual")
scatter!(ax2, ts, observed_so_far; color = :grey, markersize = 6,
    label = "observed so far")
axislegend(ax2; position = :lt)

ax3 = Axis(fig[3, 1]; xlabel = "Reference day", ylabel = "Observed-so-far")
vspan!(ax3, n - 6, n; color = (:grey, 0.15))
ci_ribbon!(ax3, ts, yt; color = :teal, label = "posterior predictive")
scatter!(ax3, ts, observed_so_far; color = :black, markersize = 6,
    label = "observed")
axislegend(ax3; position = :lt)
fig
```

<<<<<<< ours
In the shaded recent window the naive fit (red) dips below the complete-data
reference (black) — the artefactual late down-turn — while the
[`RightTruncate`](@ref) correction (green) pulls the recent ``R_t`` back up
towards the reference, having accounted for the not-yet-reported counts.

## Reading the shared parameters

Wrapping the error model in [`RightTruncate`](@ref) does not touch the renewal
process, so the corrected fit recovers the *same* shared parameters: the
autoregressive damping ``\rho`` (`damp_AR`), the innovation scale
``\sigma`` (`std`), the observation overdispersion
(`cluster_factor`), and the initial infections (`init_incidence`). Each is
namespaced by the component slot that samples it, so a prior's inner variables
never collide across the model.

=======
In the shaded recent window the naive ``R_t`` (crimson) dips below the
complete-data reference (black), while the two corrections lift the recent
``R_t`` back up off that spurious decline.
The nowcast panel shows the same story on the count scale, the corrections lifting
the reconstructed eventual totals above the truncated observed-so-far towards the
true eventual line.
The posterior-predictive panel confirms the [`RightTruncate`](@ref) fit reproduces
the observed-so-far series it saw.

## Shared parameters

Neither correction touches the renewal process, so both recover the *same* shared
parameters, the autoregressive damping ``\rho`` (`damp_AR[1]`), the innovation
scale ``\sigma`` (`std`), the observation overdispersion (`cluster_factor`) and
the initial infections (`init_incidence`).
>>>>>>> theirs
`sample` returns a [FlexiChains](https://github.com/penelopeysm/FlexiChains.jl)
chain that `summarystats` summarises directly, giving point estimates and their
uncertainty alongside the effective sample size and ``\hat R``.

```@example nowcast
using MCMCChains
summarystats(tri_chain)
```

## Prior versus posterior

Sampling the reporting-triangle model with [`Prior`](https://turinglang.org/) and
overlaying it on the posterior with
[PairPlots.jl](https://sefffal.github.io/PairPlots.jl/) confirms the joint
correction still identifies the shared parameters from the thinned triangle.

```@example nowcast
using PairPlots

prior_chain = sample(tri_post, Prior(), 1000; progress = false)
pp_keys = [@varname(damp_AR), @varname(std),
    @varname(cluster_factor), @varname(init_incidence)]
pairplot(
    PairPlots.Series(tri_chain[pp_keys]; label = "posterior"),
    PairPlots.Series(prior_chain[pp_keys]; label = "prior"))
```

## References

```@bibliography
Pages = ["realtime-nowcast.md"]
Canonical = false
```
