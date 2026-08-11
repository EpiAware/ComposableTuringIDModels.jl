#!/usr/bin/env julia
# MANAGED by EpiAwarePackageTools.scaffold — do not edit by hand.
#
# Runic format check, run in this isolated environment so its JuliaSyntax pin
# does not clash with JET. Checks the standard source trees without modifying
# them; exits non-zero if any file is not formatted. Calls `Runic.main`
# directly rather than `julia -m Runic`: the latter needs Julia >= 1.12, but
# CI also runs the LTS floor.
#
#   julia --project=test/formatter test/formatter/runtests.jl

using Runic

# Project root is two levels up from test/formatter.
project_root = dirname(dirname(@__DIR__))
dirs = filter(
    isdir,
    [
        joinpath(project_root, d)
            for d in ("src", "test", "docs", "benchmark", "ext")
    ]
)

if isempty(dirs)
    println("no source directories found to format-check")
    exit(0)
end

exit(Runic.main(["--check", "--diff", dirs...]))
