#src MANAGED by EpiAwarePackageTools.scaffold — do not edit by hand.
#src The AD page: a cost report, plus the short backend-choice note that
#src carries the `ad-backends` anchor. Lives under `docs/src/benchmarks/`
#src with its own top-level "Benchmarks" nav group, alongside the
#src performance-over-time page (when the package has one) — not nested
#src inside Tutorials (#305, the shape EpiAwareADTools#28 asked for).
#src
#src There was a separate `ad-backends.jl` tutorial under Getting started
#src carrying a support table and a how-to-choose narrative. It is retired:
#src the support table duplicated the README's per-backend coverage badge
#src row and this page's own per-backend scenario coverage, and the generic
#src advice belongs to the ecosystem rather than to any one package
#src (epiaware.github.io#28). Package pages across the org link to its
#src `@ref ad-backends` anchor, so that anchor now lives on the "Choosing a
#src backend" section below and must stay there.
#src
#src The page body is re-applied on every update so it stays kit-current;
#src everything package-specific it reports (scenarios, backends,
#src broken/skip declarations) is read at docs-build time from the
#src package-owned `test/ADFixtures` registry, so declare a broken scenario
#src there, never here. If this page cannot execute for this package, park
#src it via `FORCE_STUB_TUTORIALS` in `docs/docs_config.jl` instead of
#src editing it.

md"""
# [AD backend comparison](@id ad-comparison)

What each AD backend costs on ComposableTuringIDModels.jl's shared AD scenario set, so a
decision between backends can be made from numbers rather than the general
pattern alone.

Which backends are supported is the per-backend badge row in the
[README](https://github.com/EpiAware/ComposableTuringIDModels.jl#readme), which reports each backend's
CI status and its coverage flag from the gradient suite.

## Packages used
"""

md"""
```@raw html
<details><summary>Show setup code</summary>
```
"""

using ComposableTuringIDModels
using ADFixtures
import DifferentiationInterfaceTest as DIT
## DIT 0.11 dropped its Chairmarks dependency; `benchmark_differentiation`
## needs it loaded explicitly to resolve `run_benchmark!`.
using Chairmarks
using DataFramesMeta
using Statistics
using CairoMakie

CairoMakie.activate!(type = "png", px_per_unit = 2)
set_theme!(theme_latexfonts(); fontsize = 14)

## A DataFrame is `showable` as `text/html`, and both Literate and
## DocumenterVitepress take that branch first — so returning one from a cell
## drops DataFrames' own styled `<table>` (inline styles, a `Row` index
## column, a column-type row, an `N×M DataFrame` caption) straight into the
## page as raw HTML, outside VitePress's table styling. Wrapping the text in
## a type that is showable ONLY as `text/markdown` makes both writers emit a
## plain pipe table instead, which VitePress renders as a native table.
struct MarkdownTable
    text::String
end
Base.show(io::IO, ::MIME"text/markdown", t::MarkdownTable) = print(io, t.text)

## Render `df` as a markdown pipe table: first column left-aligned (the
## label), the rest right-aligned (the numbers). A `|` inside a cell would
## otherwise split it into two columns, so escape it -- registry backend
## names are free text.
_cell(x) = replace(string(x), "|" => "\\|")

function markdown_table(df)
    cols = string.(names(df))
    io = IOBuffer()
    println(io, "| ", join(_cell.(cols), " | "), " |")
    println(io, "|:---|", repeat("---:|", max(length(cols) - 1, 0)))
    for row in eachrow(df)
        println(io, "| ", join((_cell(row[c]) for c in cols), " | "), " |")
    end
    return MarkdownTable(String(take!(io)))
end

backend_entries = ADFixtures.backends()
scenario_list = ADFixtures.scenarios()

## The registry's optional bookkeeping accessors (see the ADRegistry
## contract): a missing accessor means no broken or skipped scenarios.
function _optional(name, default)
    return isdefined(ADFixtures, name) ? getfield(ADFixtures, name)() : default
end
global_broken = Set(String.(_optional(:broken_scenario_names, String[])))
backend_broken = _optional(
    :backend_broken_scenarios, Dict{String, Set{String}}()
)
backend_skip = _optional(
    :backend_skip_scenarios, Dict{String, Set{String}}()
);

md"""
```@raw html
</details>
```
"""

md"""
## Benchmark

`DifferentiationInterfaceTest.benchmark_differentiation` runs every
(backend, scenario) pair the registry supports.
Combinations declared broken or skipped in the registry are excluded from
their backend's rows, so they show up as reduced scenario coverage in the
`Scenarios` column below, rather than as timings of gradients that are
wrong or crash.
The figures are the prepared per-call cost.
DifferentiationInterface prepares each backend once, recording a tape for
ReverseDiff and compiling a rule for Enzyme and Mooncake, and we time the
reused operator, so that one-off preparation is excluded.
This matches repeated use such as an MCMC run, where preparation is
amortised over many gradient calls.
Each backend's time and allocations are then divided by the ForwardDiff
value on the same scenario, so ForwardDiff sits at 1.0 by construction;
values below 1.0 are faster (or lighter), above 1.0 slower (or heavier).
Timings use short per-measurement budgets so the page stays cheap to
build; treat small differences as indicative rather than exact.
"""

md"""
### Summary

Geometric mean of the relative cost across the scenarios each backend can
handle. `Scenarios` reports coverage, since a partial backend averages
only over the scenarios it differentiates.
"""

md"""
```@raw html
<details><summary>Show benchmark code</summary>
```
"""

bench_parts = map(backend_entries) do entry
    excluded = union(
        global_broken,
        get(backend_broken, entry.name, Set{String}()),
        get(backend_skip, entry.name, Set{String}())
    )
    scens = filter(s -> !(s.name in excluded), scenario_list)
    part = DataFrame(
        DIT.benchmark_differentiation(
            [entry.backend], scens;
            logging = false,
            benchmark_test = false,
            benchmark_seconds = 0.5
        )
    )
    ## Label rows with the registry's backend name, which distinguishes
    ## configurations (e.g. Enzyme forward vs reverse) that share a package.
    part[!, :backend_label] .= entry.name
    part
end
raw_bench = vcat(bench_parts...)

bench_long = @chain raw_bench begin
    @rsubset :operator == ^(:gradient)
    @rtransform begin
        :backend = :backend_label
        :scenario = :scenario.name
        :time_us = :time * 1.0e6
        :bytes_kb = :bytes / 1024
    end
    @rsubset isfinite(:time_us) && isfinite(:bytes_kb)
    @select :backend :scenario :time_us :bytes_kb
end;

## The baseline every cost is divided by: ForwardDiff when the registry has
## it (the org standard), otherwise the registry's first backend.
baseline = any(e -> e.name == "ForwardDiff", backend_entries) ?
    "ForwardDiff" : first(backend_entries).name

ref = @chain bench_long begin
    @rsubset :backend == baseline
    @select :scenario :ref_time = :time_us :ref_bytes = :bytes_kb
end

rel = @chain bench_long begin
    leftjoin(ref, on = :scenario)
    @rsubset !ismissing(:ref_time) && !ismissing(:ref_bytes)
    @rtransform begin
        :rel_time = :time_us / :ref_time
        :rel_bytes = :bytes_kb / :ref_bytes
    end
end;

## Geometric mean over positive values; guards against a zero-allocation
## scenario sending `log` to -Inf.
function geomean(x)
    pos = filter(>(0), x)
    return isempty(pos) ? NaN : exp(mean(log.(pos)))
end

n_total = length(scenario_list)

summary_table = @chain rel begin
    @by :backend begin
        :rel_time = round(geomean(:rel_time); digits = 2)
        :rel_bytes = round(geomean(:rel_bytes); digits = 2)
        :scenarios = "$(length(:scenario))/$(n_total)"
    end
    @orderby :rel_time
    rename(
        :backend => "Backend",
        :rel_time => "Relative time",
        :rel_bytes => "Relative allocations",
        :scenarios => "Scenarios"
    )
end;

md"""
```@raw html
</details>
```
"""

markdown_table(summary_table)

md"""
### Spread across scenarios

Each box summarises a backend's relative cost across the scenario set, on
a log scale so speed-ups and slow-downs are symmetric around the baseline
at 1.0.
"""

md"""
```@raw html
<details><summary>Show plotting code</summary>
```
"""

plot_df = @chain rel begin
    stack(
        [:rel_time, :rel_bytes],
        variable_name = :metric, value_name = :value
    )
    @rsubset isfinite(:value) && :value > 0
    @rtransform begin
        :metric = :metric == "rel_time" ? "Relative time" :
            "Relative allocations"
        :family = first(split(:backend))
        :mode = occursin("reverse", lowercase(:backend)) ? "reverse" :
            "forward"
    end
end

## Facet order: time then allocations. Plain CairoMakie rather than
## AlgebraOfGraphics -- the grammar-of-graphics `mapping`/`visual` calls pull
## in DimensionalData via Makie, which conflicts with FlexiChains' compat
## range in any package that hard-deps both (kit#283).
metric_order = ["Relative time", "Relative allocations"]

fig_relative = Figure(size = (1200, 500))
for (col, metric) in enumerate(metric_order)
    sub = @rsubset plot_df :metric == metric
    backend_order = sort(unique(sub.backend))
    ax = Axis(
        fig_relative[1, col];
        title = metric,
        ylabel = col == 1 ? "Cost relative to $baseline" : "",
        yscale = log10,
        xticks = (1:length(backend_order), backend_order),
        xticklabelrotation = pi / 4
    )
    xs = [findfirst(==(b), backend_order) for b in sub.backend]
    boxplot!(ax, xs, sub.value)
end

md"""
```@raw html
</details>
```
"""

fig_relative

md"""
### Per scenario

The same data with one point per scenario, so individual outliers show
rather than being summarised.
Scenarios on the horizontal axis, relative cost on the vertical axis (log
scale), backends by colour, faceted by metric.
"""

md"""
```@raw html
<details><summary>Show plotting code</summary>
```
"""

families = sort(unique(plot_df.family))
modes = sort(unique(plot_df.mode))
palette = Makie.wong_colors()
marker_shapes = [:circle, :utriangle, :rect, :diamond, :star5]

## Axes built up front (one assignment per binding, not mutated in the loop
## below) so a top-level `@example` block -- which runs each statement in
## global scope -- can't hit Julia's soft-scope "ambiguous assignment in a
## for loop" trap.
scenario_orders = [
    sort(unique((@rsubset plot_df :metric == m).scenario))
        for m in metric_order
]
fig_scenarios = Figure(size = (1600, 800))
axes_scenarios = [
    Axis(
            fig_scenarios[1, col];
            title = metric_order[col],
            ylabel = col == 1 ? "Cost relative to $baseline" : "",
            yscale = log10,
            xticks = (
                1:length(scenario_orders[col]),
                scenario_orders[col],
            ),
            xticklabelrotation = pi / 4
        )
        for col in eachindex(metric_order)
]

for (col, metric) in enumerate(metric_order)
    sub = @rsubset plot_df :metric == metric
    scenario_order = scenario_orders[col]
    ax = axes_scenarios[col]
    for (fi, fam) in enumerate(families), (mi, mode) in enumerate(modes)
        grp = @rsubset sub :family == fam && :mode == mode
        isempty(grp) && continue
        xs = [findfirst(==(s), scenario_order) for s in grp.scenario]
        scatter!(
            ax, xs, grp.value;
            color = palette[mod1(fi, length(palette))],
            marker = marker_shapes[mod1(mi, length(marker_shapes))],
            markersize = 11,
            label = "$fam ($mode)"
        )
    end
end
Legend(
    fig_scenarios[1, length(metric_order) + 1], axes_scenarios[1];
    merge = true, unique = true, title = "Backend family / Mode"
);

md"""
```@raw html
</details>
```
"""

fig_scenarios

md"""
The full long-format result is available as `raw_bench` if you want GC
fraction, compile fraction, the `value_and_gradient` rows, or absolute
timings.

## [Choosing a backend](@id ad-backends)

The numbers above are this package's scenarios, but the shape of them is
general.
Forward mode (ForwardDiff, Enzyme forward, Mooncake forward) costs one pass
per parameter, so it wins when the parameter count is small.
Reverse mode (ReverseDiff, Enzyme reverse, Mooncake reverse) costs one pass
per output regardless of the parameter count, so it pays off once this
package's quantities sit inside a larger model with many latent parameters.
Turing's
[AD guidance](https://turinglang.org/docs/usage/automatic-differentiation/)
puts the crossover around 20 parameters.

ForwardDiff is the simplest fast default below that and needs no
configuration.
Above it, switch to a reverse-mode backend through the sampler's `adtype`,
for example `sample(model, NUTS(; adtype = AutoMooncake()), 1000)`.
The surest choice is to benchmark the backends on your own model.

Where the registry enables Enzyme, the standard configuration defers
per-value activity decisions to runtime:

```julia
using ADTypes, Enzyme
AutoEnzyme(mode = Enzyme.set_runtime_activity(Enzyme.Reverse))
```

Runtime activity is not free. On paths that do not need it, it can make
Enzyme several times slower, so where one Enzyme configuration is applied
to every scenario the rows for it above are conservative.

A backend can match the reference gradient to machine precision and still
kill the process once a sampler drives it.
The AD suite therefore also runs a short NUTS smoke test per backend, in
serial and threaded modes, over a composed model from the same registry as
the gradients above.

When a backend misbehaves, start with ForwardDiff: it fails with ordinary
Julia `MethodError`s that point at the offending call, where Enzyme and
Mooncake report at the compiled-IR level.
`test/ad/run_selected.jl` checks a single (backend, scenario) pair without
running the full suite:

```
julia --project=test/ad test/ad/run_selected.jl --backend enzyme --scenario AR
```

A combination that is genuinely broken is declared in the `ADFixtures`
registry (`backend_broken_scenarios`, or `backend_skip_scenarios` when it
cannot run at all), which excludes it here and marks it `@test_broken` in
the gradient tests rather than leaving the suite red.

## Reproducing this page

The numbers above are measured on the docs-build machine, so they reflect
that CPU.
To regenerate locally:

```
task docs
```

or, equivalently:

```
julia --project=docs docs/make.jl
```

## See also

- `test/ad/` holds the gradient tests as tagged `@testitem`s, validated
  against a ForwardDiff reference with
  `DifferentiationInterfaceTest.test_differentiation`. Pass a backend tag
  (e.g. `TAG=enzyme_reverse task test-ad-backend`) to run a single backend,
  as the per-backend CI does.
- `test/ADFixtures` is the package-owned registry this page renders from;
  scenarios, backends, and broken/skip declarations all live there.
- The shared harness and the `ADRegistry` contract live in
  [EpiAwarePackageTools.jl](https://github.com/EpiAware/EpiAwarePackageTools.jl).
"""
