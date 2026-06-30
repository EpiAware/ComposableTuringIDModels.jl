using Pkg: Pkg
Pkg.instantiate()

using Documenter
using DocumenterVitepress
using DocumenterCitations
using EpiAwarePrototype

include("pages.jl")

# Doctest / example setup shared across docstrings.
DocMeta.setdocmeta!(
    EpiAwarePrototype, :DocTestSetup,
    :(using EpiAwarePrototype, Distributions); recursive = true)

# --- generate the API reference pages from the module's documented bindings ---
# Each binding is listed ONCE in a `@docs` block so a function with several
# documented methods appears under a single heading with a single index entry.
# Public vs internal mirrors Documenter's `@autodocs` Public/Private split.
function _is_public(mod::Module, sym::Symbol)
    return @static if isdefined(Base, :ispublic)
        Base.ispublic(mod, sym)
    else
        Base.isexported(mod, sym)
    end
end

function api_bindings(mod::Module)
    meta = Base.Docs.meta(mod)
    vars = sort!([b.var for b in keys(meta)]; by = string)
    public = Symbol[]
    private = Symbol[]
    for v in vars
        v === nameof(mod) && continue  # the module's own docstring lives on the home page
        push!(_is_public(mod, v) ? public : private, v)
    end
    return public, private
end

function write_api_page(path, title, anchor, page, intro, mod, names)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, anchor === nothing ? "# $title" : "# [$title](@id $anchor)")
        println(io)
        println(io, intro)
        println(io)
        println(io, "## Index")
        println(io)
        println(io, "```@index")
        println(io, "Pages = [\"$page\"]")
        println(io, "```")
        println(io)
        println(io, "```@docs")
        for name in names
            println(io, string(mod, ".", name))
        end
        println(io, "```")
    end
end

let (public, private) = api_bindings(EpiAwarePrototype)
    lib_dir = joinpath(@__DIR__, "src", "lib")
    write_api_page(
        joinpath(lib_dir, "public.md"),
        "Public API", "public-api", "public.md",
        "Documentation for the exported, public interface of `EpiAwarePrototype`.",
        EpiAwarePrototype, public)
    write_api_page(
        joinpath(lib_dir, "internals.md"),
        "Internal API", "internal-api", "internals.md",
        "Documentation for `EpiAwarePrototype`'s unexported internal helpers and " *
        "supertypes. These are not part of the stable public API; they are " *
        "documented because the public docstrings cross-reference them.",
        EpiAwarePrototype, private)
    println("Generated API pages: $(length(public)) public, $(length(private)) internal")
end

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :numeric)

makedocs(;
    modules = [EpiAwarePrototype],
    authors = "Sam Abbott and contributors",
    sitename = "EpiAwarePrototype.jl",
    # A prototype: keep the build resilient rather than strict. Missing/internal
    # docstrings and "see also" @refs to undocumented helpers warn rather than
    # fail the build.
    warnonly = [:docs_block, :missing_docs, :autodocs_block, :cross_references],
    doctest = true,
    pages = pages,
    plugins = [bib],
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "github.com/EpiAware/EpiAwarePrototype.jl",
        devbranch = "main",
        devurl = "dev",
        deploy_url = "epiawareprototype.epiaware.org"))

DocumenterVitepress.deploydocs(;
    repo = "github.com/EpiAware/EpiAwarePrototype.jl",
    target = joinpath(@__DIR__, "build"),
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true)
