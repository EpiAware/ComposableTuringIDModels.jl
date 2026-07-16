# Dependency-free pretty-printing for model components (replaces upstream's
# PrettyPrinting.jl tree display).
#
# A composed model is a tree of components: an `IDModel` holds an infection and
# an observation model, an infection model owns a latent process, a modifier
# wraps an inner model, and so on. The default Julia `show` renders such a value
# as its full nested parametric type (`RandomWalk{Normal{...}, ...}` — a
# screenful), which is unreadable. Here we render the component *tree* instead:
# each node is the role it plays plus its concrete component name, and we recurse
# only through the slots that are themselves components, leaving distributions,
# data and step structs as leaves. This is display-only; it changes nothing about
# how a model is constructed or sampled.

# Turn a field name into the role label shown in the tree, dropping a trailing
# `_model`/`_models` so `infection_model` reads as `infection` and `latent_model`
# as `latent`, while a plain name such as `rt` or `models` is kept as-is.
function _role_label(f::Symbol)
    s = String(f)
    for suffix in ("_models", "_model")
        if endswith(s, suffix) && length(s) > length(suffix)
            return s[1:(end - length(suffix))]
        end
    end
    return s
end

# Whether a value is a component the tree should recurse into. Components are
# `AbstractComposableModel`s. A raw prior slot holds a bare `Distribution` (or a
# vector of them), which is not a component and so is a leaf automatically; a
# richer prior (a latent model used as a prior, e.g. `RandomWalk`) is itself a
# component and is still recursed into.
_is_tree_component(v) = v isa AbstractComposableModel

# Collect the component children of `model` as `(role, child)` pairs. A field is a
# child when its value is a tree component (see [`_is_tree_component`](@ref)), or a
# vector/tuple containing components (as combined/stacked models hold); leaf fields
# (bare-distribution priors, generation intervals, step structs, functions) are
# skipped so the tree stays compact.
function _component_children(model::AbstractComposableModel)
    children = Tuple{String, AbstractComposableModel}[]
    for f in fieldnames(typeof(model))
        v = getfield(model, f)
        role = _role_label(f)
        if _is_tree_component(v)
            push!(children, (role, v))
        elseif v isa TimeVarying && _is_tree_component(v.prior)
            # A `TimeVarying`-marked slot is a leaf marker wrapping a prior
            # process; recurse into that process so the tree still shows it.
            push!(children, (role, v.prior))
        elseif v isa Union{AbstractVector, Tuple}
            for (i, el) in enumerate(v)
                _is_tree_component(el) &&
                    push!(children, (string(role, "[", i, "]"), el))
            end
        end
    end
    return children
end

# Recursively print the component children beneath an already-printed node, using
# box-drawing connectors and an accumulated `prefix` for indentation.
function _print_component_tree(io::IO, children, prefix::AbstractString)
    n = length(children)
    for (i, (role, child)) in enumerate(children)
        is_last = i == n
        print(io, '\n', prefix, is_last ? "└─ " : "├─ ",
            role, ": ", nameof(typeof(child)))
        child_prefix = string(prefix, is_last ? "   " : "│  ")
        _print_component_tree(io, _component_children(child), child_prefix)
    end
    return nothing
end

# Rich (REPL / `display`) rendering: the concrete component name followed by the
# recursive, indented component tree.
function Base.show(io::IO, ::MIME"text/plain", model::AbstractComposableModel)
    print(io, nameof(typeof(model)))
    _print_component_tree(io, _component_children(model), "")
    return nothing
end

# Compact (nested / `print` / `repr`) rendering: just the concrete component name,
# never the nested parametric type. Keeps a model legible inside arrays, tuples
# and error messages without dumping its whole type signature.
Base.show(io::IO, model::AbstractComposableModel) = print(io, nameof(typeof(model)))
