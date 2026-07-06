# PACKAGE-OWNED — scaffold writes this once and never overwrites it.
#
# QA configuration values the managed `quality.jl` testset reads. Fill in the
# package-specific inputs the shared helpers need; the standard testset logic
# stays in `quality.jl` (managed). Edit freely.

using ComposableTuringIDModels

const QA_CONFIG = (
    # The module under test.
    mod = ComposableTuringIDModels,

    # Path to the isolated JET environment (see test/jet/Project.toml).
    jet_env = joinpath(@__DIR__, "..", "jet"),

    # Per-check Aqua relaxations, e.g. (; ambiguities = false). Empty = all on.
    # `undocumented_names` is off because the package reexports all of
    # Distributions and Turing (`@reexport`), pulling their public names into
    # this module's namespace; those names are documented in their home
    # packages, not here, so the in-namespace undocumented-names check is not
    # meaningful for a reexporting package.
    aqua = (; undocumented_names = false),

    # ExplicitImports `ignore`: symbols an extension legitimately imports
    # non-publicly. Tuple of Symbols, e.g. (:_internal_helper,).
    ei_ignore = (),

    # Docstring `crossref_ignore`: upstream names docstrings link to via
    # `[`name`](@ref)`, e.g. (:pdf, :cdf, :logpdf).
    crossref_ignore = (),

    # Extra docstring-format options, e.g.
    # (; exported_only_examples = true, require_field_docs = true).
    docstring = (;),

    # Package extensions to ambiguity-check. Each entry:
    #   (; name = :MyPkgSomeTriggerExt,
    #      triggers = ("SomeTrigger",),       # packages to load first
    #      prefixes = ("MyPkg", "SomeTrigger"),
    #      expect_phantoms = false,    # true if a third party adds phantoms
    #      broken = false)             # true to quarantine a known ambiguity
    extensions = ()
)
