# Composable design

`EpiAwarePrototype` treats an epidemiological model as a composition of
independent parts. Each part answers one question:

  - a **latent model** describes an unobserved process ``Z_t`` over time;
  - an **infection model** maps that latent process to unobserved infections
    ``I_t``;
  - an **observation model** maps infections to the observed data ``y_t``.

Every part is a plain struct that implements a single method of the generic
constructor [`as_turing_model`](@ref). There is no deep type hierarchy: a part
is identified by the method it implements, not by its place in a tree. This is
the central design change from the package this prototype is adapted from, which
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
that `AR` to produce an ARIMA-style process, and the whole thing can drive a
[`DirectInfections`](@ref) process observed with a [`NegativeBinomialError`](@ref).

## Swap-in, swap-out

Because the parts share an interface, you change a modelling assumption by
swapping one struct for another, leaving the rest untouched:

```@example design
using EpiAwarePrototype, Distributions
data = EpiData([0.2, 0.3, 0.5], exp)

# An ARIMA-style latent process: a differenced AR.
latent = DiffLatentModel(; model = AR(), init_priors = [Normal(), Normal()])

# Swap the observation model without touching the rest.
poisson_model = EpiAwareModel(latent,
    DirectInfections(; data = data, initialisation_prior = Normal()),
    PoissonError())

negbin_model = EpiAwareModel(latent,
    DirectInfections(; data = data, initialisation_prior = Normal()),
    NegativeBinomialError())
nothing # hide
```

## Inference

A composed model is an ordinary Turing model. Pass observed data instead of
`missing` to condition it, then sample:

```julia
using Turing
y = rand(as_turing_model(poisson_model, missing, 30)).generated_y_t
posterior = as_turing_model(poisson_model, y, 30)
chain = sample(posterior, NUTS(), 1_000)
```

The standard Turing tools — `rand` for prior draws, `fix` to pin parameters,
`condition` (or `|`) to condition on values, and `sample` for inference — all
apply unchanged.
