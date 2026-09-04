# A composed model paired with the data it is fitted to, so the shape of the
# infection process comes from the data rather than being stored on the model.

@doc raw"
A composed model and the data it is fitted to, held together.

An [`IDModel`](@ref) says what the process is. The data says how long it runs
and how many streams it has. `IDProblem` is the pair, so the two travel
together and cannot drift apart. From it,
[`as_turing_model(problem)`](@ref as_turing_model) builds the
`DynamicPPL.Model` and [`data_requirements(problem)`](@ref data_requirements)
reports what the model asks of the data, neither of which needs the length
restating.

The shape of the infection process is read from the observation model and the
data at build time (see [`infection_strata`](@ref)), not stored. A plain vector
gives a single-series infection process, while a `strata x time` matrix or a
`NamedTuple` of streams gives a stratified one.

Printing an `IDProblem` shows the component tree and a summary of the data.
That is the reason to hold one rather than a conditioned `DynamicPPL.Model`,
which renders as its full nested parametric type with the observations dumped
inline.

Construction rejects data the model cannot fit, but the check is narrow and
worth knowing the shape of. It cannot see a length disagreement on a single
series, nor across streams that imply their own strata, because in both the
observation count is read from the same data it would be checked against. Ask
properly with [`data_requirements`](@ref).

There is no data-free `IDProblem`, because the pairing is what the type is.
[`IDModel`](@ref) is the object for a model on its own. A problem whose
observations have not arrived yet is one over a blank series,
`Vector{Missing}(missing, n)`, which fixes the shape and carries the length the
way observations do. The data is attached the way anything else in a
composition is respecified.

```julia
using Accessors
problem = IDProblem(model, Vector{Missing}(missing, 20))  # shape fixed, none observed
fitted = @set problem.data = y                            # observations attached
```

# Arguments

  - `model`: the composed [`IDModel`](@ref), or an infection model and an
    observation model to compose into one.
  - `data`: the observations the model is fitted to. A vector, a
    `strata x time` matrix, or a `NamedTuple` of streams.

# Examples
```@example IDProblem
using ComposableTuringIDModels, Distributions
problem = IDProblem(
    DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
    PoissonError(),
    Vector{Missing}(missing, 20))
```

```@example IDProblem
rand(as_turing_model(problem))
```

## Fields

  - `model`: the composed [`IDModel`](@ref).
  - `data`: the observations.
"
struct IDProblem{M <: IDModel, D}
    "The composed model."
    model::M
    "The observations the model is fitted to."
    data::D

    function IDProblem(model::M, data::D) where {M <: IDModel, D}
        _assert_data_fits(model, data)
        return new{M, D}(model, data)
    end
end

# An inner constructor, so `@set problem.data` is checked on the same terms.
# The guard fires only when the streams disagree with each other, a stream
# longer than the one the observation count is read from.
# It never rejects data `as_turing_model` would accept, because `data_fits` is
# the weaker of the two conditions.
# The build check adds an exact-length rule and the legacy-`n` rule on top.
function _assert_data_fits(model::IDModel, data)
    required = data_requirements(model, data, _data_shape(model, data))
    data_fits(required) && return nothing
    return throw(
        ArgumentError(
            "the data does not fit the model:\n" *
                sprint(show, MIME"text/plain"(), required)
        )
    )
end

function IDProblem(
        infection::AbstractInfectionModel, observation_model::AbstractObservationModel,
        data
    )
    return IDProblem(IDModel(infection, observation_model), data)
end

# The problem's lead-in and observation chain are its model's, so a requirements
# report reads the same either way.
observation_lead_in(problem::IDProblem) = observation_lead_in(problem.model)
_observation_chain(problem::IDProblem) = _observation_chain(problem.model)

# There is deliberately no `data_fits(::IDProblem)`.
# The count would be read from the same data it was checked against, so for a
# single series it is true by construction and would be false reassurance.
# `data_fits(model, y_t, n)` takes the two independently.
_problem_shape(problem::IDProblem) = _data_shape(problem.model, problem.data)

data_requirements(problem::IDProblem) = data_requirements(
    problem.model, problem.data, _problem_shape(problem)
)

# The generic method would report as if nothing had been supplied, discarding
# the data the problem exists to hold, so a bare length is refused by name.
function data_requirements(::IDProblem, ::ModelShape)
    return throw(
        ArgumentError(
            "an `IDProblem` holds its own data, so `data_requirements(problem)` " *
                "takes no length; to report over a different length or dataset " *
                "use `data_requirements(problem.model, y_t, n)`"
        )
    )
end

@doc raw"
Build the `DynamicPPL.Model` for an [`IDProblem`](@ref).

The problem already holds the data, so there is no length to restate.
`as_turing_model(problem)` is `as_turing_model(problem.model, problem.data)`.

There is deliberately no method taking data alongside the problem. An
`IDProblem` *is* a model and its data, so fitting the same model to a different
series is a different problem, not a different call on this one. Build it with
`Accessors`, which is how every other part of this package is respecified:

```julia
using Accessors
refit = @set problem.data = y_new
```

# Arguments

  - `problem`: the [`IDProblem`](@ref).
"
as_turing_model(problem::IDProblem) = as_turing_model(problem.model, problem.data)

# --- printing ---------------------------------------------------------------

# The model's components hang directly off the problem, with the data as a final
# sibling, so the pairing reads as one tree.
function Base.show(io::IO, ::MIME"text/plain", problem::IDProblem)
    print(io, "IDProblem")
    _print_component_tree(
        io, _component_children(problem.model), "";
        trailing = (string("data: ", _data_summary(problem.data)),)
    )
    return nothing
end

Base.show(io::IO, ::IDProblem) = print(io, "IDProblem")

# Data with nothing observed in it reads as a simulation whatever shape it
# takes, because a count of values that are all blank reads as if something had
# been supplied.
_blank_summary(n) = "none, $n time points (simulating from the prior)"
function _data_summary(y::AbstractVector)
    n = length(y)
    blank = count(ismissing, y)
    blank == n && return _blank_summary(n)
    counted = string(n, " observation", n == 1 ? "" : "s")
    blank == 0 && return string(counted, " (", eltype(y), ")")
    return string(counted, ", ", blank, " missing")
end

function _data_summary(y::AbstractMatrix)
    strata, n = size(y)
    all(ismissing, y) && return _blank_summary(n)
    return string(
        strata, " strat", strata == 1 ? "um" : "a", " x ", n,
        " observation", n == 1 ? "" : "s"
    )
end

# A `NamedTuple` is a stream per entry for a `Split`, and a field per entry for
# a model that reads more than the counts (a `BinomialError` takes `(y, N)`), so
# the summary names the entries without claiming which.
# `data_requirements` is what distinguishes them.
function _data_summary(y::NamedTuple)
    n = _series_time_length(y)
    all(v -> all(ismissing, v), y) && return _blank_summary(n)
    return string(n, " observations in each of ", join(keys(y), ", "))
end

# Anything else names its type, since the summary has to stay one line.
_data_summary(y) = string(typeof(y))
