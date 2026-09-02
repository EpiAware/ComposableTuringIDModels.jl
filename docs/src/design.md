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
parameter names stay flat (pass `prefix = true` to namespace a slot).
Because every component speaks the same `as_turing_model` protocol,
components nest freely: an [`AR`](@ref) process can carry a
[`HierarchicalNormal`](@ref) error model, a [`DiffLatentModel`](@ref) can wrap
that `AR` to produce an ARIMA-style process, and that whole latent process can be
folded into a [`DirectInfections`](@ref) model observed with a
[`NegativeBinomialError`](@ref).

## Parameter names, and prefixing where they clash

A component names its own parameters for what they are, not for which
component owns them: an initial value is `init`, a damping coefficient is
`damp`, a marginal standard deviation is `σ`.
Names stay short because a component knows nothing about what it will be
composed with, and a self-describing name would be wrong as soon as the same
component appeared twice.

The cost is that two components can want the same name.
A Gaussian process draws a marginal standard deviation `σ`, and so does
[`NormalError`](@ref).
Composed bare, both reach the top level of the same model and land on one
variable:

```@example prefixing
using ComposableTuringIDModels, Distributions
using DynamicPPL: DebugUtils, VarInfo

renewal(rt) = Renewal(;
    generation_time = [0.3, 0.4, 0.3], rt = rt, initialisation = Normal())
bare = as_turing_model(
    IDModel(renewal(HilbertSpaceGP()), NormalError()), fill(10.0, 20), 20)
keys(VarInfo(bare))
```

There are four names for five parameters, and the model says so:

```@example prefixing
DebugUtils.check_model(bare; error_on_failure = false)
```

That is a hard stop, not a silent error.
`sample` runs `check_model` first and refuses to run a model that fails it.

The fix is a prefix applied where the conflict is, rather than a longer name
everywhere.
[`PrefixLatentModel`](@ref) wraps a latent process so its variables are
namespaced, and [`PrefixObservationModel`](@ref) does the same for an
observation model:

```@example prefixing
prefixed = as_turing_model(
    IDModel(
        renewal(PrefixLatentModel(; model = HilbertSpaceGP(), prefix = "gp")),
        NormalError()
    ), fill(10.0, 20), 20)
keys(VarInfo(prefixed))
```

```@example prefixing
DebugUtils.check_model(prefixed; error_on_failure = false)
```

The Gaussian process's hyperparameters are now `gp.ℓ` and `gp.σ`, and the
error model keeps its own `σ`.
A chain reads each by its own path, `chain[@varname(gp.σ)]`, and a value is
pinned the same way with `fix(prefixed, (gp = (σ = 0.2,),))`.

Components that compose several children of the same kind prefix them
already, so nothing needs adding there.
[`CombineLatentModels`](@ref) names its components `Combine.1`, `Combine.2`
and so on, [`ConcatLatentModels`](@ref) uses `Concat.1`, `Concat.2`,
[`DiffLatentModel`](@ref) namespaces the process it differences under `diff`,
and [`Split`](@ref) prefixes each observation stream by its name.
A Gaussian process reached through one of those never collides.

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

## [An effect confined to a window](@id windowed-effect)

[`broadcast_window`](@ref) builds a latent model that equals its inner model on a declared index window and zero everywhere else.
Summing it onto a base path with [`CombineLatentModels`](@ref) gives an effect that is active only over that window.

```@example design
using ComposableTuringIDModels, Distributions

base = RandomWalk(init = Normal(0, 0.1), ϵ_t = IID(Normal(0, 0.05)))
windowed = CombineLatentModels(
    [base, broadcast_window(Normal(1, 0.2), 20:30)], ["", "Window"]
)
keys(rand(as_turing_model(windowed, 40)))
```

The inner slot is a PATH slot, so one helper covers a known effect and an estimated one.
A [`FixedIntercept`](@ref) gives a constant.
A bare `Distribution` gives one estimated level held across the window.
A process gives an effect that varies within the window, and it is generated over the window rather than over the whole series.

```@example design
known = broadcast_window(FixedIntercept(2.0), 20:30)
varying = broadcast_window(RandomWalk(), 20:30)
nothing # hide
```

Windows stack by adding members, and each keeps its own prefix.

```@example design
stacked = CombineLatentModels(
    [
        base,
        broadcast_window(Normal(1.5, 0.1), 1:5),
        broadcast_window(Normal(2, 0.1), 20:30),
    ], ["", "Early", "Late"]
)
keys(rand(as_turing_model(stacked, 40)))
```

A window touching either end of the series is an ordinary case.
The window is part of the broadcast rule, so it shows in the model tree.

```@example design
broadcast_window(Normal(1, 0.2), 1:10)
```

A latent process and a parameter prior share one role, so a windowed effect drops into a per-step parameter slot.
Here an [`AR`](@ref) is more strongly damped over the window and constant outside it.

```@example design
shifted = AR(
    damp = CombineLatentModels(
        [FixedIntercept(0.4), broadcast_window(Normal(0.4, 0.01), 20:30)],
        ["", "Shift"]
    ), ϵ_t = IID(Normal(0, 0.1))
)
keys(rand(as_turing_model(shifted, 40)))
```

### How it is built

[`broadcast_window`](@ref) is a [`BroadcastLatentModel`](@ref) under an [`InWindow`](@ref) rule, in the same way that [`broadcast_dayofweek`](@ref) is one under [`RepeatEach`](@ref).
The rule asks for a series as long as the window and places it, so nothing is drawn for the indices outside.

Without the helper the same effect is a transform that reads the index, summed on by [`CombineLatentModels`](@ref).
This generates the effect over the whole series and then masks it, so it is the longer route for a process, but it shows what the helper is doing.

```@example design
manual = CombineLatentModels(
    [
        base,
        TransformLatentModel(
            Intercept(Normal(1, 0.2)),
            x -> x .* in.(eachindex(x), Ref(20:30))
        ),
    ], ["", "Window"]
)
keys(rand(as_turing_model(manual, 40)))
```

[`CombineLatentModels`](@ref) sums, so a multiplicative window effect is the same composition in log space under an `exp` transform.
The window multiplies the path by `exp` of the drawn effect, and leaves it unchanged elsewhere.

```@example design
multiplied = TransformLatentModel(
    CombineLatentModels(
        [base, broadcast_window(Normal(log(2), 0.1), 20:30)], ["", "Window"]
    ), x -> exp.(x)
)
keys(rand(as_turing_model(multiplied, 40)))
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

The generation time is a slot like any other, so a modifier composes onto a
continuous generation time the constructor discretises for you.
The discretisation keywords stay available alongside the modifiers, which keeps
one discretisation path rather than a second one written outside the package.

```@example design
# The same modifier on a Gamma generation time discretised to 15 days.
discretised = Renewal(Gamma(2, 1.5), SusceptibleDepletion(1000.0);
    D_gen = 15.0, rt = RandomWalk())
nothing # hide
```

See [Renewal modifiers](@ref renewal-modifiers) for what each contributes to a fitted model.

## The seeding window

A renewal process needs infections before the modelled window starts.
By default `initialisation` is a level at ``t_0``: one value, decaying backwards
at the growth rate implied by ``\mathcal R_1``.
Wrapping a process in a [`SeedingPath`](@ref) estimates that run-up instead, so
the data inform its shape rather than ``\mathcal R_1`` fixing it.

```@example design
seeding = Renewal(; generation_time = gen_int, rt = RandomWalk(),
    initialisation = SeedingPath(RandomWalk(; init = Normal(log(50), 0.5))))
as_turing_model(seeding, 20)().I_seed
```

The window is returned as `I_seed`, whether it was drawn or decayed, and is not
part of `I_t`.

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
chain = sample(posterior, NUTS(; adtype = AutoMooncake(; config = nothing)),
    MCMCThreads(), 1_000, 2)
```

The standard Turing tools — `rand` for prior draws, `fix` to pin parameters,
`condition` (or `|`) to condition on values, and `sample` for inference — all
apply unchanged.

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
