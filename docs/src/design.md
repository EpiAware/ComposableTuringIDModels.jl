# Composable design

`ComposableTuringIDModels` treats an epidemiological model as a composition of
independent parts. A full model has two top-level parts:

  - an **infection model** generates unobserved infections ``I_t``. It *owns* a
    **latent model** internally — an unobserved process ``Z_t`` (e.g. a log
    reproduction number or growth rate) that the infection model maps to ``I_t``;
  - an **observation model** maps infections to the observed data ``y_t``.

A **latent model** describes an unobserved process ``Z_t`` over time. It is not
a mandatory top-level component: the latent (e.g. ``\log R_t``) is not always
the estimand, so it is folded into the infection model that consumes it. This
gives more flexibility — you choose an infection model and hand it whatever
latent process you want to drive it — and decouples the generation interval so
that only the [`Renewal`](@ref) model carries one.

Every part is a plain struct that implements a single method of the generic
constructor [`as_turing_model`](@ref). There is no deep type hierarchy: a part
is identified by the method it implements, not by its place in a tree.

## One constructor, composed by submodels

`as_turing_model(component, args...)` returns a `DynamicPPL.Model`. A component
that contains another component builds the inner model and samples it as a
submodel:

```julia
z ~ as_turing_submodel(inner_model, n)
```

`as_turing_submodel` disables automatic variable prefixing by default, so
parameter names stay flat (pass `prefix = true` to namespace a slot). Because every component speaks the same `as_turing_model` protocol,
components nest freely: an [`AR`](@ref) process can carry a
[`HierarchicalNormal`](@ref) error model, a [`DiffLatentModel`](@ref) can wrap
that `AR` to produce an ARIMA-style process, and that whole latent process can be
folded into a [`DirectInfections`](@ref) model observed with a
[`NegativeBinomialError`](@ref).

## The latent is folded into the infection model

The latent process is supplied to the infection model rather than to the
composer. An infection model takes a latent slot — `Z` for [`DirectInfections`](@ref),
`rt` for [`ExpGrowthRate`](@ref) and [`Renewal`](@ref) — and generates that
process internally before mapping it to infections. So `as_turing_model` for an
infection model takes only a series length and returns `(; I_t, Z_t)`: the
infection path and the internal latent draw, kept accessible as a generated
quantity. Only [`Renewal`](@ref) needs a generation interval, so it alone takes
one; the others take a `transformation` directly.

## Swap-in, swap-out

Because the parts share an interface, you change a modelling assumption by
swapping one struct for another, leaving the rest untouched:

```@example design
using ComposableTuringIDModels, Distributions

# An ARIMA-style latent process: a differenced AR.
latent = DiffLatentModel(; model = AR(), init = [Normal(), Normal()])

# Fold the latent into a direct-infections process, then swap the observation
# model without touching the rest.
poisson_model = IDModel(
    DirectInfections(; Z = latent, initialisation = Normal()),
    PoissonError())

negbin_model = IDModel(
    DirectInfections(; Z = latent, initialisation = Normal()),
    NegativeBinomialError())
nothing # hide
```

## Composing accumulation steps

The recurrences that drive the time series — a random walk, an autoregression, a
renewal process — are expressed as [`accumulate_scan`](@ref) steps. Within the
renewal family these steps compose: a [`RenewalStep`](@ref) is a
force-of-infection core plus an ordered tuple of modifiers that share one
incidence window. The first modifier, [`SusceptibleDepletion`](@ref), scales the
proposed incidence by the available susceptible fraction and depletes the pool,
turning the renewal process into one with a fixed population. [`Renewal`](@ref) is
a step-composing helper — pass the modifier and it is composed onto the step:

```@example design
gen_int = [0.2, 0.3, 0.5]

# A renewal process with a fixed population of 1000 and susceptible depletion.
depleting = Renewal(gen_int, SusceptibleDepletion(1000.0); rt = RandomWalk())
nothing # hide
```

A scan step is a deterministic function, so a modifier that needs *sampled*
parameters — an importation rate, say — draws them before the scan runs. That is
the modifier's optional `as_turing_model(mod, n)` method: it samples through the
same [`as_turing_submodel`](@ref) seam as everything else and returns the
modifier the scan uses. Modifiers that sample nothing return themselves, so the
step resolves its whole modifier tuple through one call and nothing in the
renewal model tests what a modifier is. [`ImportedCases`](@ref) is the worked
example, and modifiers apply in the order given:

```@example design
# Susceptible depletion, then importation added on top of the depleted
# incidence with its rate estimated. The importation prior is on the
# unconstrained scale, mapped onto a positive rate by the modifier.
seeded = Renewal(gen_int, SusceptibleDepletion(1000.0),
    ImportedCases(Normal(0.0, 1.0)); rt = RandomWalk())
nothing # hide
```

See [Renewal modifiers](@ref renewal-modifiers) for what each contributes to a fitted model.

## Inference

A composed model is an ordinary Turing model. Pass observed data instead of
`missing` to condition it, then sample. We set the automatic-differentiation
backend explicitly with `NUTS(; adtype = ...)`:
[Mooncake](https://chalk-lab.github.io/Mooncake.jl/) is the recommended default
for this package (see [Automatic differentiation backend](@ref ad-backends)).

```julia
using Turing, Mooncake
using ADTypes: AutoMooncake
y = as_turing_model(poisson_model, fill(missing, 30), 30)().generated_y_t
posterior = as_turing_model(poisson_model, y, 30)
chain = sample(posterior, NUTS(; adtype = AutoMooncake(; config = nothing)), 1_000)
```

The standard Turing tools — `rand` for prior draws, `fix` to pin parameters,
`condition` (or `|`) to condition on values, and `sample` for inference — all
apply unchanged.

## [How much of the data an observation chain scores](@id lead-in)

The `n` passed to `as_turing_model` is the length of the **infection** series,
not the number of observations, and the two are not always the same.
Every [`LatentDelay`](@ref) in an observation chain convolves the expected
series with a delay PMF and returns a series shorter by `length(pmf) - 1`: the
head of a convolution is only partially observed, so it is dropped rather than
fitted.
The observation-error model then right-aligns the data against what is left.

So a chain with delays scores the **last** `n - lead_in` observations and
silently ignores any earlier ones, where `lead_in` is the sum of
`length(pmf) - 1` over the chain's delays.
[`observation_lead_in`](@ref) reads that number off an assembled model, and
[`observation_coverage`](@ref) reports what a given `n` scores:

```@example design
delayed = LatentDelay(
    LatentDelay(PoissonError(), fill(1 / 15, 15)), fill(1 / 30, 30))
y = fill(10, 60)

observation_lead_in(delayed), observation_coverage(delayed, y, length(y))
```

Add the lead-in to the series length to score every observation:

```@example design
n = length(y) + observation_lead_in(delayed)
observation_coverage(delayed, y, n)
```

The same idiom sets an [`IDProblem`](@ref)'s time span,
`tspan = (1, length(y) + observation_lead_in(observation_model))`.
Length-preserving modifiers ([`Ascertainment`](@ref), [`RightTruncate`](@ref),
[`ReportTriangle`](@ref), a [`Split`](@ref)'s streams) add nothing of their own,
so the lead-in comes from the delays alone.

What counts as an observation follows the model that scores it, not the shape of
the data.
A [`BinomialError`](@ref) takes `(y = successes, N = trials)` and is counted by
its `y` field, a [`ReportTriangle`](@ref) by the reference days of its triangle,
and a [`Split`](@ref) one stream at a time.
Only the per-time-point error families right-align, so only there does a
non-zero `n_unscored` mean observations quietly dropped; a `ReportTriangle`
asserts its reference days instead, and the report warns before the assertion
fires.

The extended series length then has to travel with the model.
[`forecast`](@ref) rebuilds the model over the horizon, and defaults to the
length of `y`, so a chain fitted with the lead-in added back needs the same `n`:

```julia
n = length(y) + observation_lead_in(model)
chain = sample(as_turing_model(model, y, n), NUTS(), 1000)
fc = forecast(model, y, chain, 14; n = n)
```

An [`IDProblem`](@ref) already records that length in its `tspan`, so
`forecast(problem, y, chain, 14)` needs no keyword.
Forecasting from the shorter default asks for a shorter latent stream than the
chain holds, which `forecast` refuses.

## Infection↔observation mappings

An infection model does not have to generate a single curve. Two components
widen `I_t` to an `inf_strata x time` matrix, one row per stratum:

  - [`Stratify`](@ref) puts a stratum axis on **one shared** process,
    combining a shared path with a per-stratum deviation — a partially pooled
    panel, the "same many" case.
  - [`CombineInfections`](@ref) draws **several distinct** infection processes
    independently and stacks them — the "different many" case (a region-level
    epidemic, say, where each region's curve is genuinely its own process).

On the observation side, [`Split`](@ref) (with an optional [`StrataMap`](@ref)
weight matrix) reads that matrix and projects it onto observation streams
through an `obs_strata x inf_strata` weight, so one mechanism covers every
mapping cardinality: `map = I` is one-to-one, an aggregation row is
many-to-one, and a general matrix is many-to-many.
A panel model composes the same way as a single-group one.
Build the infection side with [`Stratify`](@ref) and pass it straight to the
plain `IDModel(infection_model, observation_model)` constructor.
No separate panel constructor is needed.
See [Partial pooling across groups](@ref tutorial-hierarchy) and [Multiple
observation streams](@ref tutorial-split) for worked examples.

## [Inspecting and updating an observation chain](@id obs-traversal)

An observation model is a nesting of modifiers around an error model, each
holding what it wraps in a field, so a chain with a few modifiers is several
levels deep.
Displaying one prints the whole tree:

```@example traversal
using ComposableTuringIDModels, Distributions
obs = LatentDelay(
    Ascertainment(
        Split((cases = PoissonError(), deaths = PoissonError())),
        FixedIntercept(log(0.1))
    ),
    [0.5, 0.3, 0.2]
)
```

Reading a component out of that structure, or putting a different one in, goes
through two functions rather than through field paths.
`ComposableTuringIDModels.wrapped_models` reports what a single component wraps,
and `ComposableTuringIDModels.observation_components` is the walk built on it —
every component of a chain, outermost first, with a [`Split`](@ref)'s branches
followed in stream order.
Both are public but not exported, so they are reached through the module name.

Locating a component is a `filter` over the walk instead of a field path:

```@example traversal
using ComposableTuringIDModels: observation_components
delays = filter(x -> x isa LatentDelay, observation_components(obs))
length(delays), first(delays).delay
```

A quantity accumulated over a chain is a `sum` over the same walk.
A quantity that does not accumulate reads the branching structure through
`wrapped_models` instead, because a `Split`'s streams run in parallel off one
expected series rather than nesting.

`ComposableTuringIDModels.rewrap` is the other direction — the same wrapper
rebuilt around new wrapped models, with everything else it holds carried across.
Swapping every component of a given type is those two together:

```@example traversal
using ComposableTuringIDModels: wrapped_models, rewrap
swap(f, m) = f(rewrap(m, map(x -> swap(f, x), wrapped_models(m))))
swapped = swap(x -> x isa PoissonError ? NegativeBinomialError() : x, obs)
```

Several components transform or derive a field on construction.
[`LatentDelay`](@ref) stores its delay PMF reversed for the convolution,
[`Ascertainment`](@ref) stores its prior already wrapped in a
[`PrefixLatentModel`](@ref), and [`Aggregate`](@ref) derives its presence mask
from its window lengths.
Each declares a `ConstructionBase.constructorof` that takes its fields as
stored, so `rewrap` and `Accessors.@set` both rebuild them correctly.

A modifier defined outside the package inherits both functions as long as it
holds what it wraps in its own field.
One that holds it somewhere the field walk cannot see, such as inside a
container, must define `wrapped_models` and `rewrap` for itself, and one whose
constructor cannot accept its own stored fields back must define
`ConstructionBase.constructorof`.
The docstrings state both contracts.
