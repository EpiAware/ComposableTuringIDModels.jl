# [Multiple observation streams: cases, deaths, and strata](@id case-study-split)

Real-time surveillance rarely watches an epidemic through a single lens.
The same infections surface as reported cases, hospital admissions, deaths, and
often each of these split by age, region, or variant.
These streams share one underlying infection process but differ in their
reporting delay, ascertainment, and noise [sherratt2021surveillance](@citep).
Fitting them jointly — one infection trajectory, several observation streams —
propagates uncertainty correctly and lets a sparse stream (deaths) borrow
strength from a dense one (cases).

This case study uses one construct, [`Split`](@ref), for every multi-stream
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
[renewal case study](@ref case-study-renewal), and observe it through two
pipelines.
Cases are a short-delay, high-ascertainment negative-binomial stream.
Deaths are a long-delay stream whose ascertainment — the infection-fatality
ratio — is itself *estimated*: each stream is a full observation model, so its
ascertainment can be a fixed fraction or, as here, a latent [`Intercept`](@ref)
model with a prior.

```@example split
using EpiAwarePrototype, Distributions, Random, Turing, Mooncake
using ADTypes: AutoMooncake
Random.seed!(1234)

data = EpiData(gen_distribution = Gamma(6.5, 0.62))
latent = AR(
    damp_priors = [truncated(Normal(0.8, 0.05), 0, 1),
        truncated(Normal(0.1, 0.05), 0, 1)],
    init_priors = [Normal(0.0, 0.2), Normal(0.0, 0.2)],
    ϵ_t = HierarchicalNormal(std_prior = HalfNormal(0.1)))
renewal = Renewal(data; rt = latent, initialisation_prior = Normal(log(100.0), 0.1))

cases = LatentDelay(
    Ascertainment(NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1)),
        FixedIntercept(log(0.6))),                     # ~60% case ascertainment
    LogNormal(1.6, 0.5))                                # short infection→report delay
deaths = LatentDelay(
    Ascertainment(NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1)),
        Intercept(Normal(log(0.015), 0.25))),          # estimated ~1.5% IFR
    LogNormal(2.8, 0.4))                                # long infection→death delay

parallel = Split((cases = cases, deaths = deaths))
```

The composed model assembles the renewal infection process and the two-stream
observation model exactly like a single-stream study.

```@example split
model = EpiAwareModel(renewal, parallel)
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
We draw a full chain with NUTS, matching the other case studies, and
differentiate with [Mooncake](https://chalk-lab.github.io/Mooncake.jl/), the
recommended backend for this package (see
[Automatic differentiation backend](@ref ad-backend)).

```@example split
ydata = (cases = y.cases, deaths = y.deaths)
posterior = as_turing_model(model, ydata, n)
chain = sample(
    posterior, NUTS(0.9; adtype = AutoMooncake(; config = nothing)), 1000;
    progress = false)
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
        cases = NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1)),
        deaths = LatentDelay(                            # case→death delay
            Ascertainment(NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1)),
                FixedIntercept(log(0.02))),
            LogNormal(2.2, 0.3)))),
    LogNormal(1.6, 0.5))
cascade_model = EpiAwareModel(renewal, cascade)
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
        Ascertainment(NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1)),
            FixedIntercept(log(0.7))), LogNormal(1.5, 0.4)),
    old = LatentDelay(
        Ascertainment(NegativeBinomialError(cluster_factor_prior = HalfNormal(0.1)),
            FixedIntercept(log(0.4))), LogNormal(1.8, 0.4))))
strata_model = EpiAwareModel(renewal, strata_obs)
strata_sim = as_turing_model(
    strata_model, (young = missing, old = missing), n)().generated_y_t
map(s -> sum(skipmissing(s)), strata_sim)                # totals per band
```

When the bands come from *different infection strata* rather than the same shared
infections, a [`StrataMap`](@ref) supplies the expected series as an
`infection-strata × time` matrix together with an
`observation-strata × infection-strata` weight matrix, and a single **template**
model is replicated once per data stream.
The mapping covers the one-to-one, many-to-one, and many-to-many infection →
observation cases with one weight matrix: an identity map is one-to-one, an
aggregation row sums several infection strata into one stream (many-to-one), and
a general matrix is many-to-many.

```@example split
template = Split(LatentDelay(PoissonError(), LogNormal(1.6, 0.5)))
young_inf = fill(200.0, n)
old_inf = fill(80.0, n)
inf_strata = permutedims(hcat(young_inf, old_inf))       # 2 × n (inf strata × time)
W = [1.0 0.0; 0.0 1.0; 1.0 1.0]                          # bands and their total
agemap = StrataMap(inf_strata, W)
age = as_turing_model(
    template, (young = missing, old = missing, total = missing), agemap)()
map(length, age.y_t)                                     # one series per stream
```

The aggregate `total` stream sees the summed expected infections of both bands,
and swapping the identity/aggregation rows of `W` for estimated weights is the
seam a partially-observed or cross-classified reporting structure grows from.

## References

```@bibliography
Pages = ["split-observations.md"]
Canonical = false
```
