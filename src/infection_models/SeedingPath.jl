# The seeding window of a renewal process, estimated as a process rather than
# decayed deterministically from a single initial level.

@doc raw"
Mark a latent process as the **seeding path** of a [`Renewal`](@ref): the whole
pre-window incidence is estimated, rather than decaying from a single initial
level.

[`Renewal`](@ref)'s `initialisation` slot is a *level*: one value (or one per
stratum) transformed to ``I_0``, with the pre-window incidence decaying at the
growth rate implied by ``\mathcal R_1`` (see [`renewal_init_window`](@ref)).
Wrapping a process in `SeedingPath` changes what that slot means. The process is
drawn at the **window shape** — `len_gen_int` values for a single series,
`n_strata × len_gen_int` for a stratified renewal (so the process needs a strata
axis, e.g. a [`Stratify`](@ref)) — transformed by the renewal's
`transformation`, and used directly as the initial incidence window. Nothing is
decayed and no growth rate is implied, so the data inform the shape of the
run-up.

The window runs oldest to newest, so its last value is the incidence at ``t_0``
and the scan's first step convolves the drawn window. It is returned as `I_seed`
in the renewal's generated quantities, whether it was drawn or decayed, so the
seeding phase can be inspected.

A geometric random walk over the seeding phase, with the default `exp`
transformation, is `SeedingPath(RandomWalk(; init = Normal(log(50), 0.5)))`.

## Fields

  - `model`: the latent process generating the (unconstrained) seeding path.

# Examples
```@example SeedingPath
using ComposableTuringIDModels, Distributions
renewal = Renewal(;
    generation_time = [0.2, 0.3, 0.5], rt = RandomWalk(),
    initialisation = SeedingPath(RandomWalk(; init = Normal(log(50), 0.5))))
as_turing_model(renewal, 10)().I_seed
```
"
struct SeedingPath{M <: AbstractPriorModel} <: AbstractPriorModel
    "The latent process generating the (unconstrained) seeding path."
    model::M
end

# A bare `Distribution` widens to a constant path, as every other PATH slot
# does, giving a flat but estimated seeding window.
SeedingPath(model::Distribution) = SeedingPath(path_prior(model))

# Only `Renewal` knows what a seeding window is, so a `SeedingPath` anywhere else
# is an error rather than a silently unwrapped process.
# One method per `ModelShape` member, so neither is ambiguous with the `Dims{2}`
# path-model guard.
function _seeding_path_misuse(m::SeedingPath)
    return error(
        "`SeedingPath` marks the seeding window of a `Renewal` and is only " *
            "meaningful in its `initialisation` slot. Pass " *
            "$(nameof(typeof(m.model))) itself in other prior slots."
    )
end

as_turing_model(m::SeedingPath, n::Int) = _seeding_path_misuse(m)
as_turing_model(m::SeedingPath, n::Dims{2}) = _seeding_path_misuse(m)

# The shape a seeding path is drawn at: the incidence window the renewal scan
# starts from, which is `len_gen_int` values per stratum.
_seeding_shape(::Int, len_gen_int) = len_gen_int
_seeding_shape(n::Dims{2}, len_gen_int) = (n[1], len_gen_int)
