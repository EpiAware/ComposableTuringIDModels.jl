# `IDProblem`: a composed model paired with the data it is fitted to.

@doc raw"
A composed model and the data it is fitted to, held together.

An [`IDModel`](@ref) says what the process is; the data says how long it runs
and how many streams it has. `IDProblem` is the pair, so the two travel
together and cannot drift apart. From it,
[`as_turing_model(problem)`](@ref as_turing_model) builds the
`DynamicPPL.Model` and [`data_requirements(problem)`](@ref data_requirements)
reports what the model asks of the data, neither of which needs the length
restating.

The shape of the infection process is read from the observation model and the
data at build time (see [`infection_strata`](@ref)), not stored: a plain vector
gives a single-series infection process, while a `strata x time` matrix or a
`NamedTuple` of streams gives a stratified one. Simulating from the prior is a
problem over a blank series, `Vector{Missing}(missing, n)`, which carries the
length the same way observations do.

Printing an `IDProblem` shows the component tree and a summary of the data,
which is the reason to hold one rather than a conditioned `DynamicPPL.Model`:
the latter renders as its full nested parametric type with the observations
dumped inline.

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

# The problem holds both halves, so it reports what the model asks of the data
# with no arguments at all.
# There is deliberately no `data_fits(::IDProblem)`: the observation count comes
# from the same data the report is checked against, so for a single series the
# answer is true by construction and would be false reassurance. Asking whether
# a dataset fits a model is `data_fits(model, y_t, n)`, where the two are
# independent.
_problem_shape(problem::IDProblem) = _data_shape(problem.model, problem.data)

data_requirements(problem::IDProblem) = data_requirements(
    problem.model, problem.data, _problem_shape(problem)
)

@doc raw"
Build the `DynamicPPL.Model` for an [`IDProblem`](@ref).

The problem already holds the data, so there is no length to restate:
`as_turing_model(problem)` is `as_turing_model(problem.model, problem.data)`.
Pass a second argument to refit the same model to different data, as in a
rolling refit over a growing series.

# Arguments

  - `problem`: the [`IDProblem`](@ref).
  - `data`: (optional) observations to use in place of the problem's own.
"
as_turing_model(problem::IDProblem) = as_turing_model(problem.model, problem.data)
as_turing_model(problem::IDProblem, data) = as_turing_model(problem.model, data)

# --- printing ---------------------------------------------------------------

# The model's components hang directly off the problem, with the data as a final
# sibling, so the pairing reads as one tree rather than as a model with a
# footnote.
function Base.show(io::IO, ::MIME"text/plain", problem::IDProblem)
    print(io, "IDProblem")
    _print_component_tree(
        io, _component_children(problem.model), "";
        trailing = (string("data: ", _data_summary(problem.data)),)
    )
    return nothing
end

Base.show(io::IO, ::IDProblem) = print(io, "IDProblem")

# One line describing what the problem holds, in the same vocabulary
# `data_requirements` uses: streams, observations and the shape they come in.
function _data_summary(y::AbstractVector)
    n = length(y)
    blank = count(ismissing, y)
    blank == n && return "none, $n time points (simulating from the prior)"
    counted = string(n, " observation", n == 1 ? "" : "s")
    blank == 0 && return string(counted, " (", eltype(y), ")")
    return string(counted, ", ", blank, " missing")
end

function _data_summary(y::AbstractMatrix)
    strata, n = size(y)
    return string(
        strata, " strat", strata == 1 ? "um" : "a", " x ", n,
        " observation", n == 1 ? "" : "s"
    )
end

# A `NamedTuple` is a stream per entry for a `Split`, and a field per entry for a
# model that reads more than the counts (a `BinomialError` takes `(y, N)`), so
# the summary names the entries without claiming which. `data_requirements` says
# which, and is where to look.
function _data_summary(y::NamedTuple)
    return string(
        _series_time_length(y), " observations in each of ", join(keys(y), ", ")
    )
end

# Anything else names its type rather than printing itself, since the point of
# the summary is to stay one line.
_data_summary(y) = string(typeof(y))
