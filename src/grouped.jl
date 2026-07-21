# Grouped / panel composition: one shared infection process observed by several
# groups, each reporting it at its own partially-pooled level. This lifts the
# hand-orchestrated panel pattern (a shared infection, a per-group effect over the
# grouping axis, and a per-group observation) into a first-class composed model so
# a multi-group epidemic is expressed by composition rather than a bespoke
# `@model` (issue #45).

@doc raw"
A composed grouped/panel epidemiological model: one shared infection process seen
by several groups, each reporting it at its own group effect.

Every group sees the *same* infection curve ``I_t`` (drawn once from
`infection_model`) but observes it at its own level ``\ell_g``. The per-group
effects come from a single `group_effect` prior over the grouping axis — a
[`Hierarchy`](@ref) for partial pooling, an [`IID`](@ref) for independent levels,
a [`RandomWalk`](@ref) for correlated ordered strata, or any other prior/latent
model. `combiner` maps the shared curve and a group's effect to that group's
expected series (default: a multiplicative effect on the exponential scale,
``\text{expected}_g = e^{\ell_g}\, I_t``), and each group is then observed through
`observation_model`.

The grouping dimension is **not** a field: it is read from the data at build time.
Sampling [`as_turing_model(model, Y)`](@ref as_turing_model) reads `n_time` and
`n_groups` from `size(Y)` (rows are time, columns are groups), draws the shared
infection once, draws the `n_groups` effects, and loops the groups applying the
combiner and the observation model. The group effect is namespaced through the
prior-slot prefix convention so its own innovations cannot collide with the
infection process's latent, and each group's observation variables are prefixed
with `group<g>`, so a full panel composes without any hand-written orchestration.

The returned generated quantities are
`(; I_t, Z_t, group_levels, generated_y_t, expected_y_t, y)`: the shared infection
path `I_t`, its internal latent `Z_t`, the per-group effects `group_levels`, the
`n_time × n_groups` matrices of sampled and pre-error expected counts, and `y`,
the per-group vector of `(; y_t, expected)` observation returns. Pass an
all-`missing` `Y` to simulate from the prior, or a data matrix to condition.

## Fields

  - `infection_model`: the shared infection process generating ``I_t`` (and its
    internal latent ``Z_t``), drawn once for all groups.
  - `group_effect`: the prior over the grouping axis generating the per-group
    effects (a [`Hierarchy`](@ref), a latent process, or a bare `Distribution`).
  - `observation_model`: the observation model mapping each group's expected
    series to its observed counts.
  - `combiner`: the function `(I_t, level_g)` mapping the shared infection curve
    and a group's scalar effect to that group's expected series.

# Examples
```@example GroupedIDModel
using ComposableTuringIDModels, Distributions
model = GroupedIDModel(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2)),
    Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5))),
    PoissonError())
Ymiss = Matrix{Union{Missing, Float64}}(missing, 12, 3)   # 12 time steps, 3 groups
sim = as_turing_model(model, Ymiss)()
size(sim.generated_y_t)
```
"
struct GroupedIDModel{
    I <: AbstractInfectionModel, G <: PriorLike,
    O <: AbstractObservationModel, F <: Function} <: AbstractComposableModel
    "Shared infection process generating ``I_t`` (and its internal latent ``Z_t``)."
    infection_model::I
    "Prior over the grouping axis generating the per-group effects."
    group_effect::G
    "Observation model mapping each group's expected series to observed counts."
    observation_model::O
    "Function `(I_t, level_g)` mapping the shared curve and a group's effect to that group's expected series."
    combiner::F

    function GroupedIDModel(infection_model::I, group_effect::G,
            observation_model::O, combiner::F) where {
            I <: AbstractInfectionModel, G <: PriorLike,
            O <: AbstractObservationModel, F <: Function}
        @assert hasmethod(combiner, Tuple{AbstractVector, Real}) "combiner must have a method for (AbstractVector, Real)"
        return new{I, G, O, F}(
            infection_model, group_effect, observation_model, combiner)
    end
end

# Default combiner: a multiplicative group effect on the exponential scale, so a
# group effect `level` scales the shared infection curve by `exp(level)`.
_grouped_combiner(I_t, level) = exp(level) .* I_t

function GroupedIDModel(infection_model::AbstractInfectionModel, group_effect,
        observation_model::AbstractObservationModel;
        combiner = _grouped_combiner)
    return GroupedIDModel(infection_model, group_effect, observation_model, combiner)
end

# Lift an existing `IDModel` (a shared infection + observation) to a panel by
# adding a per-group effect over the grouping axis.
function GroupedIDModel(idmodel::IDModel, group_effect; combiner = _grouped_combiner)
    return GroupedIDModel(idmodel.infection_model, group_effect,
        idmodel.observation_model, combiner)
end

@model function as_turing_model(model::GroupedIDModel, Y)
    n_time, n_groups = size(Y)
    # Per-group effects over the grouping axis. Prefixed through the prior-slot
    # convention (the `group_levels` left-hand name namespaces the whole submodel)
    # so the group prior's own innovations can never collide with the infection
    # process's latent — the group-requirement threading, handled internally.
    group_levels ~ as_turing_submodel(model.group_effect, n_groups; prefix = true)
    # One shared infection process, drawn once and seen by every group.
    infections ~ as_turing_submodel(model.infection_model, n_time)
    I_t = infections.I_t
    Z_t = infections.Z_t
    ys = Vector{Any}(undef, n_groups)
    exps = Vector{Any}(undef, n_groups)
    for g in 1:n_groups
        expected_g = model.combiner(I_t, group_levels[g])
        # Namespace each group's observation variables so the per-group error
        # draws stay distinct across the panel.
        og = PrefixObservationModel(model.observation_model, "group$g")
        obs_g ~ as_turing_submodel(og, Y[:, g], expected_g)
        ys[g] = obs_g.y_t
        exps[g] = obs_g.expected
    end
    generated_y_t = reduce(hcat, ys)
    expected_y_t = reduce(hcat, exps)
    y = [(; y_t = ys[g], expected = exps[g]) for g in 1:n_groups]
    return (; I_t, Z_t, group_levels, generated_y_t, expected_y_t, y)
end
