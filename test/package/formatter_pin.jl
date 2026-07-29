# PACKAGE-OWNED — not part of the managed template.
#
# The formatting check runs from the shared test environment (the managed
# `quality.jl` calls `test_formatting(QA_CONFIG.mod)` with no `env`), so the
# style it asserts is whatever JuliaFormatter that environment resolves. The
# package-owned `test/Project.toml` therefore pins it exactly.
#
# The version everything else uses is single-sourced from the kit and reaches
# this repo through MANAGED files that template-sync regenerates. This file is
# not one of them, so a kit bump moves the managed pin and leaves the shared
# environment behind, and the formatting check then asserts a style nothing
# applies. This guard turns that into an explicit failure naming both versions.
#
# Delete this once the kit runs the formatting check through the isolated
# `test/formatter` environment (as `test_linting` already does for JET), which
# makes both pins the same file.

@testitem "Quality: formatter pin" tags=[:quality] begin
    using Pkg: TOML

    root = joinpath(@__DIR__, "..", "..")
    shared = TOML.parsefile(joinpath(root, "test", "Project.toml"))
    managed = TOML.parsefile(joinpath(
        root, "test", "formatter", "Project.toml"))

    shared_pin = get(get(shared, "compat", Dict()), "JuliaFormatter", nothing)
    managed_pin = get(get(managed, "compat", Dict()), "JuliaFormatter", nothing)

    @test shared_pin !== nothing
    @test managed_pin !== nothing
    if shared_pin != managed_pin
        @warn "JuliaFormatter pins have drifted. Set the `[compat]` entry " *
              "in test/Project.toml to match test/formatter/Project.toml, " *
              "then reformat the tree with `task format`." shared_pin managed_pin
    end
    @test shared_pin == managed_pin
end
