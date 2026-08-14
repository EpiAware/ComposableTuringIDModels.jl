# MANAGED by EpiAwarePackageTools.scaffold — do not edit by hand.
# Re-run `scaffold(pkgdir(MyPackage))` to update; the scheduled template-sync
# regenerates this file. Package-specific inputs (ignore lists, extension
# names, broken quarantines) live in the package-owned `qa_config.jl`.
#
# Standard package-quality testset. Routes every generic QA check through the
# shared EpiAwarePackageTools helpers over the package's own module, supplied
# by `qa_config.jl`'s `QA_CONFIG` NamedTuple (see the template for fields).

@testitem "Quality: Aqua" tags = [:quality] begin
    using EpiAwarePackageTools
    include(joinpath(@__DIR__, "qa_config.jl"))
    test_aqua(QA_CONFIG.mod; QA_CONFIG.aqua...)
end

@testitem "Quality: ExplicitImports" tags = [:quality] begin
    using EpiAwarePackageTools
    include(joinpath(@__DIR__, "qa_config.jl"))
    test_explicit_imports(QA_CONFIG.mod; ignore = QA_CONFIG.ei_ignore)
end

@testitem "Quality: import centralisation" tags = [:quality] begin
    using EpiAwarePackageTools
    include(joinpath(@__DIR__, "qa_config.jl"))
    test_import_centralisation(QA_CONFIG.mod)
end

@testitem "Quality: docstring format" tags = [:quality] begin
    using EpiAwarePackageTools
    include(joinpath(@__DIR__, "qa_config.jl"))
    test_docstring_format(
        QA_CONFIG.mod;
        crossref_ignore = QA_CONFIG.crossref_ignore,
        QA_CONFIG.docstring...
    )
end

@testitem "Quality: README sections" tags = [:quality] begin
    using EpiAwarePackageTools
    include(joinpath(@__DIR__, "qa_config.jl"))
    # `readme` is a newer `QA_CONFIG` field; an adopter predating it has none.
    # Default to the repo-root README with standard sections (#163) rather
    # than erroring; warn so a typoed key doesn't silently revert (#188).
    cfg = if hasproperty(QA_CONFIG, :readme)
        QA_CONFIG.readme
    else
        @warn "QA_CONFIG has no `readme` field; checking the repo-root " *
            "README with the standard sections. Add one to qa_config.jl " *
            "to configure (or confirm) this."
        (; path = joinpath(@__DIR__, "..", ".."))
    end
    test_readme_sections(
        cfg.path;
        (k => v for (k, v) in pairs(cfg) if k !== :path)...
    )
end

@testitem "Quality: doctest" tags = [:quality] begin
    using EpiAwarePackageTools
    include(joinpath(@__DIR__, "qa_config.jl"))
    test_doctest(QA_CONFIG.mod)
end

@testitem "Quality: formatting" tags = [:quality] begin
    using EpiAwarePackageTools
    include(joinpath(@__DIR__, "qa_config.jl"))
    # `formatter_env` is a newer `QA_CONFIG` field; an adopter predating it
    # has none. Fall back to the in-process check, which floats with the
    # shared test environment's resolved Runic, rather than erroring; warn
    # so a typoed key doesn't silently revert (#188, #321).
    env = if hasproperty(QA_CONFIG, :formatter_env)
        QA_CONFIG.formatter_env
    else
        @warn "QA_CONFIG has no `formatter_env` field; checking formatting " *
            "in-process against the shared test environment, whose " *
            "Runic version floats with the CI Julia in use. Add one to " *
            "qa_config.jl to pin it via the isolated formatter env."
        nothing
    end
    test_formatting(QA_CONFIG.mod; env = env)
end

@testitem "Quality: linting (JET)" tags = [:quality] begin
    using EpiAwarePackageTools
    include(joinpath(@__DIR__, "qa_config.jl"))
    test_linting(QA_CONFIG.mod; env = QA_CONFIG.jet_env)
end

@testitem "Quality: extension ambiguities" tags = [:quality] begin
    using EpiAwarePackageTools
    include(joinpath(@__DIR__, "qa_config.jl"))
    for ext in QA_CONFIG.extensions
        # Load the trigger packages, then check the extension's surface.
        for trigger in ext.triggers
            Base.require(Main, Symbol(trigger))
        end
        test_ext_ambiguities(
            QA_CONFIG.mod, ext.name;
            prefixes = ext.prefixes,
            expect_phantoms = get(ext, :expect_phantoms, false),
            broken = get(ext, :broken, false)
        )
    end
end
