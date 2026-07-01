# 2D reporting-triangle nowcasting (epinowcast-style): the reference-date ×
# reporting-delay triangle data structure and the per-cell observation model.

@doc raw"
A **reporting triangle**: the reference-date × reporting-delay count matrix with
the not-yet-reported cells masked off.

`counts[t, d + 1]` is the number of events with reference day `t` first reported
at delay `d = 0, 1, …, Dmax`. A cell is **observed** at the present time `now`
iff `t + d ≤ now`; the remaining lower-right cells have not yet been reported and
are masked out by `observed`. This is the native object of epinowcast-style
nowcasting: it keeps the full joint reference-day × delay structure rather than
collapsing it to a per-reference-day observed-so-far total (the marginal that
[`RightTruncate`](@ref) conditions on — see the consistency note below).

Build one with [`define_y_t`](@ref) from either a dense matrix or a long-form
table of `(reference, delay, count)` rows; pass it as the `y_t` data to a
[`ReportTriangle`](@ref) observation model.

## Consistency with right-truncation (CDF-scaling)

Summing the **observed** cells of reference day `t` (delays `d = 0 … now − t`)
gives `Σ_d counts[t, d+1]`, the observed-so-far total — which in expectation is
`μ_t · F[(now − t) + 1]`, exactly the CDF-scaled expected observed-so-far that
[`RightTruncate`](@ref) conditions on. So the observed row-sums of the triangle
are the marginal of the joint model; the triangle additionally models how that
total is split across delays.

## Fields

  - `counts`: the `N[t, d+1]` count matrix (reference day × delay). Unobserved
    cells are ignored by [`ReportTriangle`](@ref); they may hold `missing`, `0`,
    or any placeholder.
  - `observed`: the boolean mask of reported cells (`t + d ≤ now`).
  - `Dmax`: the maximum reporting delay (the number of delay columns is
    `Dmax + 1`, delays `0 … Dmax`).

# Examples
```@example ReportingTriangle
using EpiAwarePrototype
# A 4 × 3 matrix (reference days 1..4, delays 0..2); now = 4.
N = [10 5 2; 12 6 3; 14 7 4; 16 8 5]
rt = define_y_t(ReportTriangle(PoissonError(), [0.5, 0.3, 0.2]), N, fill(20.0, 4))
rt.observed
```
"
struct ReportingTriangle{A <: AbstractMatrix, M <: AbstractMatrix{Bool}}
    "The `N[t, d+1]` reference-day × delay count matrix."
    counts::A
    "The boolean mask of reported cells (`t + d ≤ now`)."
    observed::M
    "The maximum reporting delay (delay columns are `0 … Dmax`)."
    Dmax::Int
end

# Build the `t + d ≤ now` observation mask for an `n × (Dmax + 1)` triangle.
function _triangle_mask(n::Int, Dmax::Int, now::Int)
    return BitMatrix([t + d <= now for t in 1:n, d in 0:Dmax])
end

@doc raw"
A composable reporting-delay PMF component: a `Distribution` (or a precomputed
vector) turned into the reporting-delay PMF `p` (delays `0 … Dmax`), for use as
the delay submodel of [`ReportTriangle`](@ref).

Sampling `as_turing_model(c::ReportingPMF, n)` returns the length-`(Dmax + 1)` PMF
(the argument `n` is the role-interface series length and is ignored — the PMF is
indexed by delay, not reference day). The PMF is non-negative and sums to one;
`Dmax = length(pmf) - 1`.

It is the **fixed-delay** default used by [`ReportTriangle`](@ref): the PMF is
precomputed once and held constant. Because the delay is supplied to
`ReportTriangle` as a submodel — mirroring how [`RightTruncate`](@ref) takes a
[`ReportingCDF`](@ref) — a user can instead pass any latent component producing the
delay PMF, the seam an estimated / time-varying delay grows from.

## Constructors

  - `ReportingPMF(distribution; D, Δd)` — discretise a continuous reporting-delay
    distribution via double-interval censoring (CensoredDistributions.jl), exactly
    the released-CD path [`LatentDelay`](@ref) / [`EpiData`](@ref) use.
  - `ReportingPMF(pmf)` — from a precomputed delay PMF (non-negative, sums to 1).

# Examples
```@example ReportingPMF
using EpiAwarePrototype, Distributions
c = ReportingPMF(truncated(Normal(2.0, 1.0), 0.0, Inf))
as_turing_model(c, 10)()
```

## Fields

  - `pmf`: the reporting-delay PMF `p` (delays `0 … Dmax`, non-negative, sums to
    1); `Dmax = length(pmf) - 1`.
"
struct ReportingPMF{T <: AbstractVector{<:Real}} <: AbstractLatentModel
    "The reporting-delay PMF (delays `0 … Dmax`)."
    pmf::T

    function ReportingPMF(pmf::T) where {T <: AbstractVector{<:Real}}
        @assert all(>=(0), pmf) "Delay PMF must be non-negative"
        @assert isapprox(sum(pmf), 1) "Delay PMF must sum to 1"
        return new{T}(pmf)
    end
end

function ReportingPMF(distribution::C; D = nothing, Δd = 1.0) where {
        C <: ContinuousDistribution}
    return ReportingPMF(_discretised_pmf(distribution; Δd = Δd, D = D))
end

# The maximum reporting delay carried by the PMF component (delays `0 … Dmax`).
_pmf_Dmax(c::ReportingPMF) = length(c.pmf) - 1

@model function as_turing_model(c::ReportingPMF, n)
    return c.pmf
end

@doc raw"
A **2D reporting-triangle** observation model (epinowcast-style nowcasting), the
*joint* counterpart to the [`RightTruncate`](@ref) marginal.

`ReportTriangle` consumes the expected eventual totals `Y_t = μ_t` (per reference
day, the same quantity the infection pipeline produces) together with a
reporting-delay PMF `p` drawn from a **submodel**, expands them to per-cell
expected means `μ_{t,d} = μ_t · p[d + 1]`, and scores **only the observed cells**
of a [`ReportingTriangle`](@ref) (`t + d ≤ now`) under a per-cell count error
model. The not-yet-reported cells are never sampled. Because the model's `Y_t`
stays the eventual total, the nowcast of the eventual total is just `Y_t` read out
as a generated quantity (and the completed triangle is `μ_{t,d}` over *all* `d`).

The delay PMF is supplied as a **composable submodel** `delay_model` — mirroring
how [`RightTruncate`](@ref) takes a [`ReportingCDF`](@ref) — sampled with
`to_submodel(..., false)` inside the model. The default [`ReportingPMF`](@ref)
wraps a fixed PMF (the fixed-delay, independent-cell variant: each observed cell is
an independent Poisson / negative-binomial draw about its mean, via the per-cell
[`AbstractObservationErrorModel`](@ref) supplied as `error_model`). Because the
delay is a submodel, an estimated / time-varying delay — and the multinomial-split
parameterisation — are the seams this grows from (a later phase of the nowcasting
design).

The default PMF is built with the same released-CD `double_interval_censored` +
`pdf` discretisation path that [`LatentDelay`](@ref) / [`EpiData`](@ref) use, so
the triangle's per-cell means and the right-truncation nowcast share one delay
kernel.

## Constructors

  - `ReportTriangle(error_model, delay_model)` — from a delay submodel producing
    the PMF (e.g. a [`ReportingPMF`](@ref)); `Dmax` is read from it.
  - `ReportTriangle(error_model, pmf::AbstractVector)` — wrap a precomputed delay
    PMF (non-negative, sums to 1) in the default [`ReportingPMF`](@ref).
  - `ReportTriangle(error_model, distribution; D, Δd)` — discretise a continuous
    reporting-delay distribution via double-interval censoring
    (CensoredDistributions.jl) into a [`ReportingPMF`](@ref), exactly as
    [`LatentDelay`](@ref).

## The `y_t` data contract

The observation data is a [`ReportingTriangle`](@ref), built through the shared
[`define_y_t`](@ref) hook from either a matrix or a long-form table:

```julia
y_t = define_y_t(obs, N, Y_t)                       # from a count matrix
y_t = define_y_t(obs, reports, Y_t; now = now)      # from (reference, delay, count) rows
```

Pass `y_t = missing` to **simulate**: `ReportTriangle` builds a fully observed
triangle (`now = n + Dmax`) of `missing` cells and fills them predictively.

## Fields

  - `error_model`: the per-cell count-error model (e.g. [`PoissonError`](@ref),
    [`NegativeBinomialError`](@ref)).
  - `delay_model`: the delay submodel producing the reporting-delay PMF `p`
    (delays `0 … Dmax`); cell `(t, d)` has expected mean `Y_t[t] · p[d + 1]`. The
    default [`ReportingPMF`](@ref) holds a fixed PMF.

# Examples
```@example ReportTriangle
using EpiAwarePrototype, Distributions
obs = ReportTriangle(PoissonError(), truncated(Normal(2.0, 1.0), 0.0, Inf))
# Simulate a triangle for 15 reference days of expected total 50.
sim = as_turing_model(obs, missing, fill(50.0, 15))()
sim.observed
```
"
@kwdef struct ReportTriangle{E <: AbstractObservationErrorModel,
    D <: AbstractLatentModel} <: AbstractObservationModel
    "The per-cell count-error model."
    error_model::E
    "The delay submodel producing the reporting-delay PMF (delays `0 … Dmax`)."
    delay_model::D
end

function ReportTriangle(error_model::E,
        pmf::T) where {
        E <: AbstractObservationErrorModel, T <: AbstractVector{<:Real}}
    return ReportTriangle(error_model, ReportingPMF(pmf))
end

function ReportTriangle(error_model::E, distribution::C; D = nothing,
        Δd = 1.0) where {
        E <: AbstractObservationErrorModel, C <: ContinuousDistribution}
    return ReportTriangle(error_model, ReportingPMF(distribution; D = D, Δd = Δd))
end

# The maximum reporting delay carried by the model (delays `0 … Dmax`). Read
# statically from the delay submodel so `define_y_t` can size the triangle before
# the PMF is sampled.
_triangle_Dmax(o::ReportTriangle) = _pmf_Dmax(o.delay_model)

@doc raw"
Build the [`ReportingTriangle`](@ref) data a [`ReportTriangle`](@ref) scores.

Dispatched on [`ReportTriangle`](@ref) as the triangle method of the shared
[`define_y_t`](@ref) data-unpacking hook (the vector / `NamedTuple` methods serve
the per-time-point error families). It accepts:

  - a [`ReportingTriangle`](@ref) — returned unchanged (already built);
  - a dense matrix `N[t, d+1]` — the present time defaults to `now = size(N, 1)`
    (the last reference day), masking cell `(t, d)` observed iff `t + d ≤ now`;
  - a long-form table of `(reference, delay, count)` rows (any `Tables.jl` table,
    e.g. a `DataFrame`) — accumulated into the matrix, with `now` supplied as a
    keyword;
  - `missing` — a fully observed triangle of `missing` cells for predictive
    simulation (`now = n + Dmax`, so every delay of every reference day is
    reported).

# Arguments

  - `obs_model`: the [`ReportTriangle`](@ref) model (fixes `Dmax`).
  - `y_t`: the raw data — a [`ReportingTriangle`](@ref), a matrix, a long-form
    table, or `missing`.
  - `Y_t`: the expected eventual-total series (its length `n` sets the number of
    reference days when sizing a `missing` triangle).

# Keyword Arguments

  - `now`: the present time used for the `t + d ≤ now` mask. Defaults to the
    number of reference days for a matrix, and is required for a long-form table.
  - `reference`, `delay`, `count`: the column names in a long-form table
    (defaults `:reference`, `:delay`, `:count`).

# Examples
```@example define_y_t_triangle
using EpiAwarePrototype
obs = ReportTriangle(PoissonError(), [0.5, 0.3, 0.2])      # Dmax = 2
N = [10 5 2; 12 6 3; 14 7 4]                            # 3 reference days
define_y_t(obs, N, fill(20.0, 3)).observed
```
"
function define_y_t(obs_model::ReportTriangle, y_t::ReportingTriangle, Y_t)
    Dmax = _triangle_Dmax(obs_model)
    @assert y_t.Dmax==Dmax "The triangle's Dmax ($(y_t.Dmax)) must match the model's delay PMF (Dmax = $Dmax)"
    return y_t
end

function define_y_t(obs_model::ReportTriangle, y_t::Missing, Y_t)
    n = length(Y_t)
    Dmax = _triangle_Dmax(obs_model)
    counts = Matrix{Missing}(missing, n, Dmax + 1)
    # Simulating: every cell is reported, so take `now` past the last delay of the
    # last reference day (`now = n + Dmax`) — the full rectangle is observed.
    observed = _triangle_mask(n, Dmax, n + Dmax)
    return ReportingTriangle(counts, observed, Dmax)
end

function define_y_t(obs_model::ReportTriangle, y_t::AbstractMatrix, Y_t; now::Int = size(
        y_t, 1))
    n = size(y_t, 1)
    Dmax = _triangle_Dmax(obs_model)
    @assert size(y_t, 2)==Dmax + 1 "The count matrix has $(size(y_t, 2)) delay columns; the model expects Dmax + 1 = $(Dmax + 1)"
    observed = _triangle_mask(n, Dmax, now)
    return ReportingTriangle(y_t, observed, Dmax)
end

function define_y_t(obs_model::ReportTriangle, reports, Y_t; now::Int,
        reference::Symbol = :reference, delay::Symbol = :delay,
        count::Symbol = :count)
    # Long-form `(reference, delay, count)` rows → the dense `N[t, d+1]` matrix.
    n = length(Y_t)
    Dmax = _triangle_Dmax(obs_model)
    rows = rowtable(reports)
    N = zeros(Int, n, Dmax + 1)
    for r in rows
        t = getproperty(r, reference)
        d = getproperty(r, delay)
        (1 <= t <= n) || error("reference day $t out of range 1:$n")
        (0 <= d <= Dmax) || error("delay $d out of range 0:$Dmax")
        N[t, d + 1] += getproperty(r, count)
    end
    return define_y_t(obs_model, N, Y_t; now = now)
end

@model function as_turing_model(obs_model::ReportTriangle, y_t, Y_t)
    n = length(Y_t)

    # Draw the reporting-delay PMF from the delay submodel (the default
    # `ReportingPMF` returns a fixed PMF; a richer component could sample it).
    pmf ~ to_submodel(as_turing_model(obs_model.delay_model, n), false)

    # Per-cell error priors (e.g. the NegBin cluster factor) are shared across all
    # observed cells, sampled once from the inner error family.
    priors ~ to_submodel(
        generate_observation_error_priors(obs_model.error_model, missing, Y_t), false)

    # Build (or pass through) the reporting triangle, then read the mask before
    # rebinding `y_t` to the bare count matrix. The tilde below scores `y_t[…]`
    # directly — using the model *argument* name — so DynamicPPL treats concrete
    # cells as conditioned observations (not latent variables to link) and a
    # `missing` matrix as predictive draws to fill in (widening/rebinding `y_t`).
    tri = define_y_t(obs_model, y_t, Y_t)
    observed = tri.observed
    Dmax = tri.Dmax
    @assert size(observed, 1)==n "The triangle has $(size(observed, 1)) reference days; Y_t has $n"
    y_t = tri.counts

    for t in 1:n, d in 0:Dmax

        observed[t, d + 1] || continue
        # Expected cell mean: eventual total × reporting-delay PMF at this delay.
        μ = Y_t[t] * pmf[d + 1] + 1e-6
        y_t[t, d + 1] ~ observation_error(obs_model.error_model, μ, priors...)
    end
    return ReportingTriangle(y_t, observed, Dmax)
end
