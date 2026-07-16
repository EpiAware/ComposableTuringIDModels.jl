# Composable design

`ComposableTuringIDModels` treats an epidemiological model as a composition of
independent parts. A full model has two top-level parts:

  - an **infection model** generates unobserved infections ``I_t``. It *owns* a
    **latent model** internally — an unobserved process ``Z_t`` (e.g. a log
    reproduction number or growth rate) that the infection model maps to ``I_t``;
  - an **observation model** maps infections to the observed data ``y_t``.

A **latent model** describes an unobserved process ``Z_t`` over time. It is no
longer a mandatory top-level component: the latent (e.g. ``\log R_t``) is not
always the estimand, so it is folded into the infection model that consumes it.
This gives more flexibility — you choose an infection model and hand it whatever
latent process you want to drive it — and decouples the generation interval so
that only the [`Renewal`](@ref) model carries one.

Every part is a plain struct that implements a single method of the generic
constructor [`as_turing_model`](@ref). There is no deep type hierarchy: a part
is identified by the method it implements, not by its place in a tree. This is
the central design change from the package this one is adapted from, which
used separate `generate_latent`, `generate_latent_infs`, and
`generate_observations` functions dispatched over a layered abstract hierarchy.

## One constructor, composed by submodels

`as_turing_model(component, args...)` returns a `DynamicPPL.Model`. A component
that contains another component builds the inner model and samples it as a
submodel:

```julia
z ~ to_submodel(as_turing_model(inner_model, n), false)
```

The trailing `false` disables automatic variable prefixing, so parameter names
stay flat. Because every component speaks the same `as_turing_model` protocol,
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
quantity. Only [`Renewal`](@ref) needs a generation interval, so it alone carries
an [`IDData`](@ref); the others take a `transformation` directly.

## Swap-in, swap-out

Because the parts share an interface, you change a modelling assumption by
swapping one struct for another, leaving the rest untouched:

```@example design
using ComposableTuringIDModels, Distributions

# An ARIMA-style latent process: a differenced AR.
latent = DiffLatentModel(; model = AR(), init_priors = [Normal(), Normal()])

# Fold the latent into a direct-infections process, then swap the observation
# model without touching the rest.
poisson_model = IDModel(
    DirectInfections(; Z = latent, initialisation_prior = Normal()),
    PoissonError())

negbin_model = IDModel(
    DirectInfections(; Z = latent, initialisation_prior = Normal()),
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
data = IDData([0.2, 0.3, 0.5], exp)

# A renewal process with a fixed population of 1000 and susceptible depletion.
depleting = Renewal(data, SusceptibleDepletion(1000.0); rt = RandomWalk())
nothing # hide
```

Depletion is not baked into the renewal recurrence: it is a separate modifier
composed onto the core, so further renewal-family mechanisms can be added the same
way without rewriting the model.

## Inference

A composed model is an ordinary Turing model. Pass observed data instead of
`missing` to condition it, then sample. We set the automatic-differentiation
backend explicitly with `NUTS(; adtype = ...)`:
[Mooncake](https://chalk-lab.github.io/Mooncake.jl/) is the recommended default
for this package (see [Automatic differentiation backend](@ref ad-backend)).

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
