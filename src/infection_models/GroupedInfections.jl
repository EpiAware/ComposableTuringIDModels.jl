# One shared infection process replicated across a grouping axis by a
# per-group effect, into one groups x time `I_t` matrix — the "same many"
# many-to-many infection-side mapping. Composes with `IDModel` and `Split`
# like every other infection model.

@doc raw"
Replicate one shared infection process across a grouping axis by a per-group
effect, into an `n_groups x n_time` `I_t` matrix.

Every group sees the *same* infection curve ``I_t`` (drawn once from
`infection_model`) but at its own level ``\ell_g``. The per-group effects come
from a single `group_effect` prior over the grouping axis — a
[`Hierarchy`](@ref) for partial pooling, an [`IID`](@ref) for independent
levels, a [`RandomWalk`](@ref) for correlated ordered strata, or any other
prior/latent model. `combiner` maps the shared curve and a group's effect to
that group's row of the returned matrix (default: a multiplicative effect on
the exponential scale, ``\text{row}_g = e^{\ell_g}\, I_t``).

The grouping dimension is **not** a field: it is supplied at build time via
`n`, exactly the way a series length is passed to `as_turing_model(latent,
n)`. Because a group count sits alongside the series length here,
`as_turing_model(model, n)` takes `n` as a `(n_time, n_groups)` NamedTuple
rather than a bare `Int` — the infection role's `n` argument is duck-typed to
whatever shape a particular infection model needs.

The returned `Z_t` is `(; Z_t, group_levels)`: the shared infection model's
own internal latent (kept under the same name so a plain infection model and
a `GroupedInfections` both expose their latent under `.Z_t`), plus the
per-group effects, so both stay recoverable as generated quantities through
[`IDModel`](@ref)'s unchanged, generic body.

The group prior is namespaced through the prior-slot prefix convention
(`prefix = true`, matching every other prior slot in the package) so its own
innovations can never collide with the shared infection process's latent; the
shared infection process itself stays flat (prefix off), the default for
every infection model.

Combine with [`Split`](@ref) on the observation side — `Split(template)`
observes every group through the same observation model; `IDModel`'s outer
constructors (`IDModel(infection_model, group_effect, observation_model)`)
build this pair directly.

## Fields

  - `infection_model`: the shared infection process generating ``I_t`` (and
    its internal latent ``Z_t``), drawn once for all groups.
  - `group_effect`: the prior over the grouping axis generating the per-group
    effects (a [`Hierarchy`](@ref), a latent process, or a bare
    `Distribution`).
  - `combiner`: the function `(I_t, level_g)` mapping the shared infection
    curve and a group's scalar effect to that group's row.

# Examples
```@example GroupedInfections
using ComposableTuringIDModels, Distributions
model = GroupedInfections(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2)),
    Hierarchy(; mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5))))
sim = as_turing_model(model, (n_time = 12, n_groups = 3))()
size(sim.I_t)   # 3 groups x 12 time points
```
"
struct GroupedInfections{I <: AbstractInfectionModel, G <: PriorLike, F <: Function} <:
       AbstractInfectionModel
    "Shared infection process generating ``I_t`` (and its internal latent ``Z_t``)."
    infection_model::I
    "Prior over the grouping axis generating the per-group effects."
    group_effect::G
    "Function `(I_t, level_g)` mapping the shared curve and a group's effect to that group's row."
    combiner::F

    function GroupedInfections(infection_model::I, group_effect::G,
            combiner::F) where {
            I <: AbstractInfectionModel, G <: PriorLike, F <: Function}
        @assert hasmethod(combiner, Tuple{AbstractVector, Real}) "combiner must have a method for (AbstractVector, Real)"
        return new{I, G, F}(infection_model, group_effect, combiner)
    end
end

# Default combiner: a multiplicative group effect on the exponential scale, so
# a group effect `level` scales the shared infection curve by `exp(level)`.
_grouped_combiner(I_t, level) = exp(level) .* I_t

function GroupedInfections(infection_model::AbstractInfectionModel, group_effect;
        combiner = _grouped_combiner)
    return GroupedInfections(infection_model, group_effect, combiner)
end

@model function as_turing_model(
        model::GroupedInfections, n::NamedTuple{(:n_time, :n_groups)})
    n_time, n_groups = n.n_time, n.n_groups
    # The shared infection process is drawn once, flat (prefix off) — the
    # default for every infection model.
    infections ~ as_turing_submodel(model.infection_model, n_time)
    I_t = infections.I_t
    # Per-group effects, namespaced through the prior-slot prefix convention so
    # the group prior's own innovations can never collide with the infection
    # process's own latent.
    group_levels ~ as_turing_submodel(model.group_effect, n_groups; prefix = true)
    rows = [model.combiner(I_t, group_levels[g]) for g in 1:n_groups]
    I_mat = permutedims(reduce(hcat, rows))
    Z_t = (; Z_t = infections.Z_t, group_levels)
    return (; I_t = I_mat, Z_t)
end
