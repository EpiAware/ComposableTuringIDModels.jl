# [Multiple observation streams: cases, deaths, and strata](@id tutorial-split)

Real-time surveillance rarely watches an epidemic through a single lens.
The same infections surface as reported cases, hospital admissions, deaths, and
often each of these split by age, region, or variant.
These streams share one underlying infection process but differ in their
reporting delay, ascertainment, and noise [sherratt2021surveillance](@citep).
Fitting them jointly — one infection trajectory, several observation streams —
propagates uncertainty correctly and lets a sparse stream (deaths) borrow
strength from a dense one (cases).

This tutorial uses one construct, [`Split`](@ref), for every multi-stream
shape.
`Split` observes the expected series arriving at the point where it sits in the
pipeline through several named streams, so *where you place it* chooses the
composition:

  - **parallel** — placed high, on infections: every stream observes the *same*
    ``I_t`` (cases and deaths each a delayed, ascertained fraction of ``I_t``);
  - **cascade** — placed low, after a shared layer: a later stream is observed
    *downstream* of an earlier one (deaths as a delayed fraction of the
    *expected reported cases*);
  - **strata** — one stream per data-defined group (an age band).

## How `Split` threads streams

Every observation model in the package returns the uniform pair
`(; y_t, expected)`: the sampled observations `y_t` and the pre-error `expected`
series the error was scored against.
Exposing `expected` is what lets `Split` do all three shapes with one mechanism.
`Split` feeds each stream the `expected` series reaching it, and — because
`Split` is itself an observation model — a shared modifier can run *before* it.
`Split((cases = …, deaths = …))` on its own splits infections (parallel), while
`LatentDelay(Split((cases = …, deaths = …)), pmf)` applies a common delay first
and then splits, so a stream nested inside another stream's pipeline sits
downstream of it (cascade).

!!! note "The threaded quantity is the expected, not the realised, series"
    A downstream stream reads its upstream stream's **expected** (pre-error)
    series, never its realised, sampled counts.
    So a cascade threads the *mean* reported cases into deaths, not a noisy draw.
    The case where an observation depends on another stream's *realised*
    (error-corrupted) observation — feeding sampled cases, not expected cases,
    into deaths — is not covered here and is out of scope for now.

`Split` also prefixes each stream's sampled variables with the stream name
automatically, so the streams stay distinct without any manual prefix layer.

## Parallel: cases and deaths from shared infections

We drive the streams with a renewal infection process, exactly as in the
[renewal tutorial](@ref tutorial-renewal), and observe it through two
pipelines.
Cases are a short-delay, high-ascertainment negative-binomial stream.
Deaths are a long-delay stream whose ascertainment — the infection-fatality
ratio — is itself *estimated*: each stream is a full observation model, so its
ascertainment can be a fixed fraction or, as here, a latent [`Intercept`](@ref)
model with a prior.

```@example split
using ComposableTuringIDModels, Distributions, Random, Turing, Mooncake
using ADTypes: AutoMooncake
Random.seed!(1234)

latent = AR(
    damp = [truncated(Normal(0.8, 0.05), 0, 1),
        truncated(Normal(0.1, 0.05), 0, 1)],
    init = [Normal(0.0, 0.2), Normal(0.0, 0.2)],
    ϵ_t = HierarchicalNormal(std = HalfNormal(0.1)))
renewal = Renewal(; generation_time = Gamma(6.5, 0.62),
    rt = latent, initialisation = Normal(log(100.0), 0.1))

cases = LatentDelay(
    Ascertainment(NegativeBinomialError(cluster_factor = HalfNormal(0.1)),
        FixedIntercept(log(0.6))),                     # ~60% case ascertainment
    LogNormal(1.6, 0.5))                                # short infection→report delay
deaths = LatentDelay(
    Ascertainment(NegativeBinomialError(cluster_factor = HalfNormal(0.1)),
        Intercept(Normal(log(0.015), 0.25))),          # estimated ~1.5% IFR
    LogNormal(2.8, 0.4))                                # long infection→death delay

parallel = Split((cases = cases, deaths = deaths))
```

The composed model assembles the renewal infection process and the two-stream
observation model exactly like a single-stream study.

```@example split
model = IDModel(renewal, parallel)
```

Passing `missing` data simulates a synthetic outbreak.
The per-stream data contract is a `NamedTuple` keyed by stream name, and the
returned `generated_y_t` is a `NamedTuple` of the two simulated series.

```@example split
n = 70
sim = as_turing_model(model, (cases = missing, deaths = missing), n)()
y = sim.generated_y_t
(total_cases = sum(skipmissing(y.cases)), total_deaths = sum(skipmissing(y.deaths)))
```

Fitting conditions on both streams at once.
We draw two chains in parallel with `MCMCThreads()`, matching the other
tutorials, and differentiate with
[Mooncake](https://chalk-lab.github.io/Mooncake.jl/), the recommended
backend for this package (see
[Automatic differentiation backend](@ref ad-backends)).

```@example split
ydata = (cases = y.cases, deaths = y.deaths)
posterior = as_turing_model(model, ydata, n)
chain = sample(
    posterior, NUTS(0.95; adtype = AutoMooncake(; config = nothing)),
    MCMCThreads(), 250, 2; progress = false)
nothing # hide
```

The two streams keep their own overdispersion parameters — `Split` prefixes them
`cases.cluster_factor` and `deaths.cluster_factor` — while sharing the one
infection trajectory, and the deaths stream's estimated IFR intercept
(`deaths.Ascertainment.intercept`) is recovered alongside them.
The dense case stream pins the shared ``R_t`` process; the sparse death stream is
observed jointly rather than fit in isolation.

```@example split
using MCMCChains
summarystats(chain)
```

Summary statistics confirm the parameters converged, but they do not show
whether the fit actually tracks the two simulated series.
Posterior-predictive draws from the same fitted chain, plotted against the
simulated counts, close that loop.
The helpers below turn a `time × draws` matrix into 50% and 95% credible
bands and draw a median line with ribbons.
`predictive_bands` walks a named stream's per-index sampled `y_t`, filling
any index a stream's delay leaves unscored with `missing`.

```@setup split
using CairoMakie, Statistics

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

# posterior-predictive bands for one stream's sampled y_t, indexed through a
# varname-building closure (e.g. `i -> @varname(cases.y_t[i])`)
function predictive_bands(pred, n, vn)
    ndraws = length(vec(pred[vn(n)]))
    rows = map(1:n) do i
        try
            permutedims(vec(pred[vn(i)]))
        catch
            fill(missing, 1, ndraws)
        end
    end
    credible_bands(reduce(vcat, rows))
end
```

`Split` prefixes each stream's `y_t`, so `predictive_bands` reads
`cases.y_t[i]` and `deaths.y_t[i]` off the `predict` chain, the prefixed
equivalent of the bare `y_t[i]` the [renewal tutorial](@ref tutorial-renewal)
reads from a single-stream model:

```@example split
missmodel = as_turing_model(model, (cases = missing, deaths = missing), n)
pred = predict(missmodel, chain)
cases_bands = predictive_bands(pred, n, i -> @varname(cases.y_t[i]))
deaths_bands = predictive_bands(pred, n, i -> @varname(deaths.y_t[i]))

fig = Figure(; size = (760, 620))
ax1 = Axis(fig[1, 1]; ylabel = "Cases")
ci_ribbon!(ax1, 1:n, cases_bands; color = :teal,
    label = "posterior predictive")
scatter!(ax1, 1:n, y.cases; color = :black, markersize = 6,
    label = "simulated")
axislegend(ax1; position = :lt)
ax2 = Axis(fig[2, 1]; xlabel = "Day", ylabel = "Deaths")
ci_ribbon!(ax2, 1:n, deaths_bands; color = :firebrick,
    label = "posterior predictive")
scatter!(ax2, 1:n, y.deaths; color = :black, markersize = 6,
    label = "simulated")
axislegend(ax2; position = :lt)
fig
```

The band covers the simulated series on almost every day for both streams,
and both medians track the outbreak's rise and fall rather than sitting flat
at the mean.
The sparser death series (82 simulated deaths against 6576 cases over the
same 70 days) still recovers: its band is visibly wider, but it moves with
the same underlying trajectory rather than needing its own signal to do so.

## Cascade: deaths downstream of reported cases

In the parallel model, cases and deaths both branch off infections, so a
reporting artefact in the case series (a weekend dip, an ascertainment change)
does *not* touch deaths.
Sometimes we want the opposite: deaths modelled as a delayed fraction of the
*reported cases*, so whatever is reflected in cases propagates into deaths.
That is a cascade ``I_t \to \text{cases} \to \text{deaths}``, and it needs no new
construct and no mode flag — it is the same [`Split`](@ref) placed *lower* in the
stack.
Share the infection→case-report delay, then split: the cases stream applies its
error to the delayed expectation, and the deaths stream sits downstream, delayed
again by the case-report→death interval and scaled by the fatality fraction.

```@example split
cascade = LatentDelay(                                   # infection→case delay
    Split((
        cases = NegativeBinomialError(cluster_factor = HalfNormal(0.1)),
        deaths = LatentDelay(                            # case→death delay
            Ascertainment(NegativeBinomialError(cluster_factor = HalfNormal(0.1)),
                FixedIntercept(log(0.02))),
            LogNormal(2.2, 0.3)))),
    LogNormal(1.6, 0.5))
cascade_model = IDModel(renewal, cascade)
cas = as_turing_model(cascade_model, (cases = missing, deaths = missing), n)()
```

The `Split` sits *after* the shared case delay and *before* the error leaves, so
the deaths stream's expected input is the delayed-and-ascertained *expected
cases*, not the raw infections: it is both scaled by the fatality fraction and
shortened by the case delay.

```@example split
(cases_expected_length = length(cas.expected_y_t.cases),
    deaths_expected_length = length(cas.expected_y_t.deaths),
    deaths_are_a_fraction_of_cases =
        sum(cas.expected_y_t.deaths) < sum(cas.expected_y_t.cases))
```

## Strata: one stream per age band

A stratified stream — one observation series per age band, region, or variant —
is again the same construct, here composed with the renewal infection process
and observed through one named stream per band.
Each band is a full observation model, so its delay and ascertainment can differ,
and its parameters are namespaced by the band name.

```@example split
strata_obs = Split((
    young = LatentDelay(
        Ascertainment(NegativeBinomialError(cluster_factor = HalfNormal(0.1)),
            FixedIntercept(log(0.7))), LogNormal(1.5, 0.4)),
    old = LatentDelay(
        Ascertainment(NegativeBinomialError(cluster_factor = HalfNormal(0.1)),
            FixedIntercept(log(0.4))), LogNormal(1.8, 0.4))))
strata_model = IDModel(renewal, strata_obs)
strata_sim = as_turing_model(
    strata_model, (young = missing, old = missing), n)().generated_y_t
map(s -> sum(skipmissing(s)), strata_sim)                # totals per band
```

The streams above each observe the *same* infections. When the streams instead
draw on a **weighted mix** of infections — one band, another band, and a summed
total — the same `Split` carries an `observation-strata × infection-strata` weight
matrix, and a single **template** model is replicated once per data stream.
`Split(template, W)` projects the infection series reaching it through `W`, so it
composes inside an `IDModel` like any other observation model: the infections come
from the modelled process, not a hand-built series.
One weight matrix covers the one-to-one (an identity map), many-to-one (an
aggregation row summing infection strata into one stream), and many-to-many
(a general matrix) infection → observation cases.

Here the renewal process supplies one infection stratum, and `W` maps it onto a
`young` band, an `old` band, and their `total`:

```@example split
W = reshape([0.7, 0.3, 1.0], 3, 1)                  # young, old, and their total
weighted = Split(LatentDelay(PoissonError(), LogNormal(1.6, 0.5)), W)
weighted_model = IDModel(renewal, weighted)
age = as_turing_model(
    weighted_model, (young = missing, old = missing, total = missing), n)()
map(s -> sum(skipmissing(s)), age.generated_y_t)         # simulated total per band
```

The aggregate `total` stream sees the summed expected infections of both bands —
its expected series is exactly `young .+ old`.

Simulating checks that the forward map runs.
It does not check that a many-to-one `W` is actually **recoverable** from
data.
Fitting `weighted_model` to its own simulated streams answers that: the fit
conditions on `young`, `old`, and the aggregate `total` together, exactly as
[`Split`](@ref) conditions on any other named streams.

```@example split
weighted_data = (young = age.generated_y_t.young, old = age.generated_y_t.old,
    total = age.generated_y_t.total)
weighted_posterior = as_turing_model(weighted_model, weighted_data, n)
weighted_chain = sample(
    weighted_posterior, NUTS(0.95; adtype = AutoMooncake(; config = nothing)),
    MCMCThreads(), 250, 2; progress = false)
nothing # hide
```

`young`, `old`, and `total` are not three independent counts: `total` is
exactly `young + old`, so all three read the same one-stratum ``R_t`` path
through fixed, unequal weights rather than each pinning it independently.
That collinearity makes the ``\epsilon_t`` innovations mix more slowly than
the two-stream parallel fit above — worth knowing before trusting any one
tutorial's diagnostics at face value:

```@example split
summarystats(weighted_chain)
```

Posterior-predictive bands per stream, plotted against the simulated counts,
show whether the shared infection process and the `W` weights together
recover each stratum — including `total`, which is nowhere in the infection
process itself, only assembled from it by `W`:

```@example split
weighted_pred = predict(as_turing_model(
        weighted_model, (young = missing, old = missing, total = missing), n),
    weighted_chain)
young_bands = predictive_bands(weighted_pred, n, i -> @varname(young.y_t[i]))
old_bands = predictive_bands(weighted_pred, n, i -> @varname(old.y_t[i]))
total_bands = predictive_bands(weighted_pred, n, i -> @varname(total.y_t[i]))

fig2 = Figure(; size = (760, 780))
ax_young = Axis(fig2[1, 1]; ylabel = "Young")
ax_old = Axis(fig2[2, 1]; ylabel = "Old")
ax_total = Axis(fig2[3, 1]; ylabel = "Total (young + old)", xlabel = "Day")
for (ax, bands, obs, color) in (
        (ax_young, young_bands, age.generated_y_t.young, :teal),
        (ax_old, old_bands, age.generated_y_t.old, :firebrick),
        (ax_total, total_bands, age.generated_y_t.total, :purple))
    ci_ribbon!(ax, 1:n, bands; color = color, label = "posterior predictive")
    scatter!(ax, 1:n, obs; color = :black, markersize = 6, label = "simulated")
end
axislegend(ax_young; position = :lt)
fig2
```

Despite the slower mixing, the 95% band covers the simulated counts on almost
every day for all three streams, so the many-to-one map is recovered, not
merely simulated.
A single shared ``R_t`` path, read through three collinear weighted views of
it, is enough to pin that path: `young` and `old` need no independent signal
of their own, and `total` — a stream the infection process never draws
directly — still lands on its simulated series because `W` ties it to the
same shared path.
`Split(template, W)` many-to-one aggregation is not just forward-simulated
here, it is fit and recovered.

Here the single renewal process supplied one infection stratum, broadcast
through `W`. When the strata are genuinely **separate infection processes** —
several distinct regions, say, each with its own latent — swap the single
infection model for [`CombineInfections`](@ref): it draws each process
independently and stacks the results into the same `infection-strata × time`
matrix `Split`/`StrataMap` already expect, so `IDModel(CombineInfections([...]),
Split(template, W))` maps several distinct infection processes onto streams
end-to-end. For one process carried across a strata axis instead, with
partially pooled per-stratum deviations, see [`Stratify`](@ref) and [Partial
pooling across groups](@ref tutorial-hierarchy). It also composes with
[`Renewal`](@ref)'s `mixing` slot, see [Coupled patch models](@ref
tutorial-patches).

## References

```@bibliography
Pages = ["split-observations.md"]
Canonical = false
```
