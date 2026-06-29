# Turing helper: fix-and-condition a `DynamicPPL.Model`.

@doc raw"
Condition a `DynamicPPL.Model` by fixing some parameters and conditioning on
others.

```julia
condition_model(model, fix_parameters, condition_parameters)
```

equals `condition(fix(model, fix_parameters), condition_parameters)`. Either
named tuple may be empty.

# Arguments

  - `model`: the `DynamicPPL.Model` to fix and condition.
  - `fix_parameters`: a named tuple of parameters to fix to constant values.
  - `condition_parameters`: a named tuple of parameters to condition on data.

# Examples
```@example condition_model
using EpiAwarePrototype, Distributions
m = as_turing_model(RandomWalk(), 10)
condition_model(m, (rw_init = 0.0,), NamedTuple())
```
"
function condition_model(model::DynamicPPL.Model, fix_parameters::NamedTuple,
        condition_parameters::NamedTuple)
    _model = fix(model, fix_parameters)
    _model = condition(_model, condition_parameters)
    return _model
end
