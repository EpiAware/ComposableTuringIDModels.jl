# The default Julia `show` renders a composed model as its full nested
# parametric type, which is a screenful.
# This renders the component tree instead, recursing only through slots that are
# themselves components and leaving distributions, data and step structs as
# leaves.
# It is display-only.

# `infection_model` reads as `infection`, while a plain name such as `rt` or
# `models` is kept as it is.
function _role_label(f::Symbol)
    s = String(f)
    for suffix in ("_models", "_model")
        if endswith(s, suffix) && length(s) > length(suffix)
            return s[1:(end - length(suffix))]
        end
    end
    return s
end

# A raw prior slot holds a bare `Distribution`, which is not a component and so
# is a leaf automatically.
# A latent model used as a prior is a component and is still recursed into.
_is_tree_component(v) = v isa AbstractComposableModel

# A field is a child when its value is a tree component, or a vector or tuple
# containing components.
# Leaf fields are skipped so the tree stays compact.
function _component_children(model::AbstractComposableModel)
    children = Tuple{String, AbstractComposableModel}[]
    for f in fieldnames(typeof(model))
        v = getfield(model, f)
        role = _role_label(f)
        if _is_tree_component(v)
            push!(children, (role, v))
        elseif v isa Union{AbstractVector, Tuple}
            for (i, el) in enumerate(v)
                _is_tree_component(el) &&
                    push!(children, (string(role, "[", i, "]"), el))
            end
        end
    end
    return children
end

# The label a node shows in the tree.
# A component whose behaviour is set by a non-component field rather than by its
# own type overrides this, because such fields are otherwise leaves the tree
# never reaches.
_node_label(model) = string(nameof(typeof(model)))

# Recursively print the component children beneath an already-printed node, using
# box-drawing connectors and an accumulated `prefix` for indentation.
function _print_component_tree(io::IO, children, prefix::AbstractString)
    n = length(children)
    for (i, (role, child)) in enumerate(children)
        is_last = i == n
        print(
            io, '\n', prefix, is_last ? "└─ " : "├─ ",
            role, ": ", _node_label(child)
        )
        child_prefix = string(prefix, is_last ? "   " : "│  ")
        _print_component_tree(io, _component_children(child), child_prefix)
    end
    return nothing
end

# Rich (REPL / `display`) rendering: the concrete component name followed by the
# recursive, indented component tree.
function Base.show(io::IO, ::MIME"text/plain", model::AbstractComposableModel)
    print(io, _node_label(model))
    _print_component_tree(io, _component_children(model), "")
    return nothing
end

# Compact rendering, for a model nested inside an array, a tuple or an error
# message.
Base.show(io::IO, model::AbstractComposableModel) = print(io, nameof(typeof(model)))
