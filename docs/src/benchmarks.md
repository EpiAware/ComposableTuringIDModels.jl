# [Benchmarks](@id benchmarks)

`ComposableTuringIDModels` tracks the performance of representative modelling
operations over time.
The suite is a prototype: it covers a small set of representative models rather
than exhaustively measuring every component, and the numbers are indicative
rather than a guarantee.

Benchmarking reuses the shared tooling in
`EpiAwarePackageTools.Benchmarks` rather than re-implementing a runner or a
comparison report.
The package owns the suite definition (`benchmark/benchmarks.jl`); the kit owns
running it and turning results into a legible pull-request comment.

## What is measured

The suite is a `BenchmarkTools.BenchmarkGroup` named `SUITE`, defined in
`benchmark/benchmarks.jl`, with three groups.

- **Model evaluation** — building and evaluating representative models.
  For each model the suite times a prior draw (`rand`, which samples every
  random variable) and the forward pass (`model()`, which returns the generated
  quantities).
  The models are two latent processes (`AR`, `RandomWalk`) and two composed
  `IDModel`s (`DirectInfections` with `PoissonError`, and `Renewal` with
  `NegativeBinomialError`), each turned into a Turing model via
  `as_turing_model`.
- **Sampling** — a short NUTS run (50 draws) on a composed
  `DirectInfections` + `PoissonError` model conditioned on data simulated from
  its own prior.
- **AD gradients** — the gradient of a representative log-density across
  automatic-differentiation backends (`ForwardDiff`, `ReverseDiff`, `Mooncake`,
  and `Enzyme` where supported).
  Results are keyed by scenario and backend so the comparison report folds them
  into a per-(scenario × backend) matrix.

## Running the suite locally

The `benchmark/` directory is its own Julia environment.
Run the whole suite and save the results with the managed runner, which calls
`EpiAwarePackageTools.Benchmarks.run_suite`:

```sh
julia --project=benchmark benchmark/run.jl results.json
```

`run_suite` uses a short per-benchmark time budget so a full run stays
affordable while the minimum-time estimator used in the comparison stays
stable.
To compare two result files and write a Markdown report, use the managed
comparison script, which calls
`EpiAwarePackageTools.Benchmarks.compare_comment`:

```sh
julia --project=benchmark benchmark/compare.jl pr.json base.json comment.md
```

## Continuous integration

Two workflows drive benchmarking in CI, both building on the shared kit.

- `benchmark.yaml` runs on pull requests.
  It benchmarks the pull-request head and the base branch in separate jobs,
  then posts (and updates) a single comparison comment: a bucketed summary plus
  collapsed per-benchmark tables split into evaluation and AD-gradient groups.
- `benchmark-history.yaml` runs on pushes to `main` and on tags.
  It benchmarks the recent tagged releases plus the current commit with
  AirspeedVelocity and publishes a timeline to the repository's `benchmarks`
  branch.
  The kit's `asv_comment` / `flatten_asv` helpers read the same
  AirspeedVelocity result format when a report is needed.
