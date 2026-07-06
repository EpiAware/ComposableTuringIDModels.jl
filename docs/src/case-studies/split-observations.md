# [Multiple observation streams: cases, deaths, and strata](@id case-study-split)

Real-time surveillance rarely watches an epidemic through a single lens.
The same infections surface as reported cases, hospital admissions, deaths, and
often each of these split by age, region, or variant.
These streams share one underlying infection process but differ in their
reporting delay, ascertainment, and noise, and their biases differ too
[sherratt2021surveillance](@citep).
Fitting them jointly — one infection trajectory, several observation streams —
propagates uncertainty correctly and lets a sparse stream (deaths) borrow
strength from a dense one (cases).

This case study uses one construct, [`Split`](@ref), for every multi-stream
shape.
`Split` observes one expected series through several named streams, and the only
thing that changes between the shapes is where each stream's expected input comes
from:

  - **parallel** — every stream observes the *same* infections (cases and deaths
    each a delayed, ascertained fraction of ``I_t``);
  - **sequential** — a stream is observed *downstream* of another (deaths as a
    delayed fraction of the *expected reported cases*, so a case-reporting
    artefact flows into deaths);
  - **strata** — one stream per data-defined group (an age band), with a mapping
    from infection strata to observation streams.

## How `Split` threads streams

Every observation model in the package returns the uniform pair
`(; y_t, expected)`: the sampled observations `y_t` and the pre-error `expected`
series the error was scored against.
Exposing `expected` is what lets `Split` do all three shapes with one mechanism.
A parallel stream reads the incoming infections; a sequential stream reads an
earlier stream's `expected` output; a strata stream reads its slice of a
multi-stratum expected.
Because `Split` is itself an observation model, it can be placed **high** (on
infections) or **low** (after a shared delay): `LatentDelay(Split(...), pmf)`
applies a common delay first and then splits, so a pipeline can be split at any
point.

## Parallel: cases and deaths from shared infections

We drive the streams with a renewal infection process, exactly as in the
[renewal case study](@ref case-study-renewal), and observe it through two
pipelines.
Cases are a short-delay, high-ascertainment negative-binomial stream; deaths are
a long-delay, low-ascertainment (an infection-fatality-ratio) stream.
Each stream is a full observation model in its own right; [`Split`](@ref) keys
them by name and prefixes their parameters so they stay distinct.

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
        FixedIntercept(log(0.015))),                   # ~1.5% infection-fatality ratio
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
We differentiate with [Mooncake](https://chalk-lab.github.io/Mooncake.jl/), the
recommended backend for this package (see
[Automatic differentiation backend](@ref ad-backend)).

```@example split
ydata = (cases = y.cases, deaths = y.deaths)
posterior = as_turing_model(model, ydata, n)
chain = sample(
    posterior, NUTS(0.9; adtype = AutoMooncake(; config = nothing)), 100;
    progress = false)
nothing # hide
```

The two streams keep their own overdispersion parameters — the stack prefixes
them `cases.cluster_factor` and `deaths.cluster_factor` — while sharing the one
infection trajectory.
The dense case stream pins the shared ``R_t`` process; the sparse death stream is
observed jointly rather than fit in isolation.

```@example split
using MCMCChains
summarystats(chain)
```

## Sequential: deaths downstream of reported cases

In the parallel model, cases and deaths both branch off infections, so a
reporting artefact in the case series (a weekend dip, an ascertainment change)
does *not* touch deaths.
Sometimes we want the opposite: deaths modelled as a delayed fraction of the
*reported cases*, so whatever is reflected in cases propagates into deaths.
That is a **sequential** cascade ``I_t \to \text{cases} \to \text{deaths}``, and
it is the same [`Split`](@ref) with `sequential = true`.

```@example split
sequential = Split((
        cases = LatentDelay(PoissonError(), LogNormal(1.6, 0.5)),
        deaths = LatentDelay(
            Ascertainment(PoissonError(), FixedIntercept(log(0.02))),
            LogNormal(2.2, 0.3))); sequential = true)
seq_model = EpiAwareModel(renewal, sequential)
seq = as_turing_model(seq_model, (cases = missing, deaths = missing), n)()
nothing # hide
```

The threaded quantity is the earlier stream's **expected** (pre-error) series,
never its noisy sampled output.
Reading the returned `expected_y_t` shows the deaths stream's expected input is
the delayed-and-ascertained *expected cases*, not the raw infections: it is both
scaled (by the fatality fraction) and shortened by the case delay.

```@example split
(cases_expected_length = length(seq.expected_y_t.cases),
    deaths_expected_length = length(seq.expected_y_t.deaths))
```

## Strata: one stream per age band

A stratified stream — one observation series per age band, region, or variant —
is the same construct built from a single **template** observation model instead
of named streams.
The number of streams and their names come from the *data* at model-build time
(one stream per entry of the `y_t` `NamedTuple`), so the strata count is not
hard-coded on the model.

A [`StrataMap`](@ref) supplies the expected series as an
`infection-strata × time` matrix together with an
`observation-strata × infection-strata` weight matrix, and maps infections onto
observation streams.
The mapping covers the one-to-one, many-to-one, and many-to-many infection →
observation cases with one weight matrix: an identity map is one-to-one, an
aggregation row sums several infection strata into one stream (many-to-one), and
a general matrix is many-to-many.

```@example split
template = Split(LatentDelay(PoissonError(), LogNormal(1.6, 0.5)))

# Two infection age strata over 40 days (here just illustrative constant means).
young = fill(200.0, 40)
old = fill(80.0, 40)
strata = permutedims(hcat(young, old))          # 2 × 40 (inf strata × time)

# Three observation streams: the two age bands and their aggregate total.
W = [1.0 0.0; 0.0 1.0; 1.0 1.0]                  # (obs strata × inf strata)
agemap = StrataMap(strata, W)

age_data = (young = missing, old = missing, total = missing)
age = as_turing_model(template, age_data, agemap)()
map(length, age.y_t)                             # one series per stream
```

Each stream is a prefixed copy of the template, so an age band's parameters are
namespaced (`young.`, `old.`, `total.`), and the aggregate stream sees the summed
expected infections of both bands.
Swapping the identity/aggregation rows of `W` for estimated weights is the seam a
partially-observed or cross-classified reporting structure grows from.

## References

```@bibliography
Pages = ["split-observations.md"]
Canonical = false
```
