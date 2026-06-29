# Ascertainment helper: day-of-week reporting effect.

@doc raw"
Build an [`Ascertainment`](@ref) model for a day-of-week reporting effect.

The latent model is wrapped with [`broadcast_dayofweek`](@ref) so a 7-day effect
is broadcast across the expected-observation series, and combined multiplicatively
with the expected observations by default.

# Arguments

  - `model`: the underlying observation model.

# Keyword Arguments

  - `latent_model`: the latent model broadcast over the week (default
    [`HierarchicalNormal`](@ref)`()`).
  - `transform`: the function `(x, y)` combining expected observations with the
    broadcast effect (default `(x, y) -> x .* y`).
  - `latent_prefix`: the prefix applied to the latent model's variables (default
    `\"DayofWeek\"`).

# Examples
```@example ascertainment_dayofweek
using EpiAwarePrototype
obs = ascertainment_dayofweek(PoissonError())
mdl = as_turing_model(obs, missing, fill(10.0, 14))
rand(mdl)
```
"
function ascertainment_dayofweek(model; latent_model = HierarchicalNormal(),
        transform = (x, y) -> x .* y, latent_prefix = "DayofWeek")
    return Ascertainment(
        model, broadcast_dayofweek(latent_model), transform, latent_prefix)
end
