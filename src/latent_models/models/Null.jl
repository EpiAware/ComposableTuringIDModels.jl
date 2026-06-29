# Null latent model (generates `nothing`).

@doc raw"
A null latent model that generates `nothing` (no latent variables).

# Examples
```jldoctest Null
using EpiAwarePrototype
null = Null()
mdl = as_turing_model(null, 10)
isnothing(mdl())

# output

true
```
"
struct Null <: AbstractLatentModel end

@model function as_turing_model(model::Null, n)
    return nothing
end
