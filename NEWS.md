Changes are documented in Github releases.

## Unreleased

  - `HilbertSpaceGP` and `ExactGP` now sample their hyperparameters as `gp_ℓ`
    and `gp_σ` rather than `ℓ` and `σ`. A latent process composes into its host
    prefix-off, so the old bare `σ` shared a namespace with an observation
    error model's own `σ` and one silently overwrote the other: a
    Gaussian-process `R_t` with a `NormalError` failed `check_model` and
    `sample` refused to run (issue #268).

    **Breaking**: code that reads a GP hyperparameter from a chain, or pins one
    with `fix`, must use the new names (`chain[:gp_σ]`,
    `fix(model, (gp_ℓ = 0.5, gp_σ = 0.5))`). The weights (`β` for the
    Hilbert-space basis, `z` for the exact GP) are unchanged.

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
