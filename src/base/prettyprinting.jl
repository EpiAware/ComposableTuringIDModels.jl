# Dependency-free pretty-printing for model components (replaces upstream's
# PrettyPrinting.jl tree display).

# Pretty-printing shared by every model component. The original package leant on
# PrettyPrinting.jl; here we keep a dependency-free `show` that lists the public
# fields of a component, which is enough for interactive inspection and doctests.
function Base.show(io::IO, ::MIME"text/plain", model::AbstractEpiAwareModel)
    print(io, nameof(typeof(model)))
    fields = fieldnames(typeof(model))
    if isempty(fields)
        print(io, "()")
        return nothing
    end
    println(io, ":")
    for (i, f) in enumerate(fields)
        sep = i == length(fields) ? "└─ " : "├─ "
        println(io, "  ", sep, f, " = ", repr(getfield(model, f)))
    end
    return nothing
end
