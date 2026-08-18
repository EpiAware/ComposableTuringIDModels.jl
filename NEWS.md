Changes are documented in Github releases.

## Unreleased

  - Parameter names no longer carry the name of the component that draws them.
    A component names a quantity for what it is, and where two components then
    want the same name the composition that puts them together prefixes them.

    **Breaking**: anything reading these from a chain, or pinning one with
    `fix`, needs the new name. `fix` and `condition` ignore a name the model
    does not have, so a stale call keeps running and pins nothing.

    | component | old | new |
    | --- | --- | --- |
    | `AR` | `ar_init` | `init` |
    | `AR` | `damp_AR` | `damp` |
    | `RandomWalk` | `rw_init` | `init` |
    | `Hierarchy` | `hierarchy_mean` | `mean` |
    | `ARStep` (field) | `damp_AR` | `damp` |

    Unchanged, with reasons: `DiffLatentModel`'s `latent_init` would land on
    the inner process's own `init` (`DiffLatentModel(model = AR())` fails
    `check_model` under that rename); `arima`'s `ar_init` / `diff_init`
    keyword arguments name two different init slots in one constructor;
    `init_incidence`, `cluster_factor`, `import_rates` and `std` name the
    quantity rather than the component.

  - Documented what a component's parameter names namespace to, and how to keep
    two components apart when they want the same name. A component names its
    parameters for what they are, so a Gaussian process and `NormalError` both
    draw a bare `σ`, and composing them puts two parameters on one variable:
    `check_model` fails and `sample` refuses to run. The fix is a prefix where
    the conflict is, through `PrefixLatentModel` or `PrefixObservationModel`,
    and both now have a documented home in the composable-design page along
    with the Gaussian-process docstrings. No parameter is renamed.

  - Infection models can now generate several strata at once. `Stratify` puts
    a stratum axis on a shared process (a partially pooled panel of
    reproduction numbers, say), and `Renewal`'s `mixing` slot couples the
    resulting patches, from a fixed weight matrix to an inferred `Gravity`
    model. `Split` projects any number of infection strata onto any number of
    observation streams through a weight matrix, so `IDModel` builds
    many-to-one and many-to-many infection↔observation mappings from a plain
    two-argument constructor with no separate panel API.

    **Breaking**: `GroupedInfections` and the grouped `IDModel` constructors
    are removed; a `Stratify`-based infection model observed through `Split`
    replaces them.
