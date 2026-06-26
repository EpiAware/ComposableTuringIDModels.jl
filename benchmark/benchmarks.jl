# PACKAGE-OWNED — scaffold writes this once and never overwrites it.
#
# Benchmark suite definition. Build a BenchmarkTools `BenchmarkGroup` named
# `SUITE`; the managed `run.jl` / `compare.jl` consume it. Put AD-gradient
# benchmarks under the `"AD gradients"` group so the comparison comment folds
# them into a compact per-(scenario x backend) matrix. Edit freely.

using BenchmarkTools
using EpiAwarePrototype

const SUITE = BenchmarkGroup()

# Example evaluation benchmark — replace with the package's own:
# SUITE["Evaluation"]["example"] = @benchmarkable sum(rand(100))

# Example AD-gradient group (folded into a matrix by `compare.jl`):
# SUITE["AD gradients"]["scenario"]["ForwardDiff"] = @benchmarkable ...
