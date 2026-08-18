Changes are documented in Github releases.

## Unreleased

  - Documented what a component's parameter names namespace to, and how to keep
    two components apart when they want the same name. A component names its
    parameters for what they are, so a Gaussian process and `NormalError` both
    draw a bare `σ`, and composing them puts two parameters on one variable:
    `check_model` fails and `sample` refuses to run. The fix is a prefix where
    the conflict is, through `PrefixLatentModel` or `PrefixObservationModel`,
    and both now have a documented home in the composable-design page along
    with the Gaussian-process docstrings (issue #268). No parameter is renamed.

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
