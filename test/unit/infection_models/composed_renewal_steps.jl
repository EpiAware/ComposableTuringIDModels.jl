# Proof-of-concept for #48 Phase 2: a composable renewal accumulation step.
#
# The acceptance test is behaviour-equivalence: composing the plain
# `ConstantRenewalStep` force-of-infection core with a `SusceptibleDepletion`
# modifier must reproduce the hand-fused `ConstantRenewalWithPopulationStep`
# exactly, under the same `accumulate_scan`. This is the renewal-family
# realisation of the "renewal with susceptible depletion" use case the issue
# calls out, and it is the yardstick for whether the composable contract works.

@testitem "renewal_foi matches the ConstantRenewalStep recurrence" begin
    using ComposableTuringIDModels:
                                    ConstantRenewalStep, renewal_foi
    rev_gen = reverse([0.2, 0.3, 0.5])
    step = ConstantRenewalStep(rev_gen)
    window = [1.0, 2.0, 3.0]
    Rt = 1.4
    # The extracted FOI is exactly the new-incidence term the step commits.
    @test renewal_foi(step, window, Rt) ≈ step(window, Rt)[end]
end

@testitem "ComposedRenewalStep reproduces ConstantRenewalWithPopulationStep" begin
    using ComposableTuringIDModels:
                                    ConstantRenewalStep, ConstantRenewalWithPopulationStep,
                                    ComposedRenewalStep, SusceptibleDepletion,
                                    accumulate_scan,
                                    _renewal_init_state
    rev_gen = reverse([0.2, 0.3, 0.5])
    N = 1000.0
    Rt = [1.6, 1.4, 1.2, 1.0, 0.9, 0.8]

    # Reference: the existing hand-fused with-population step.
    ref_step = ConstantRenewalWithPopulationStep(rev_gen, N)
    ref_init = _renewal_init_state(ref_step, 5.0, 0.1, length(rev_gen))
    ref = accumulate_scan(ref_step, ref_init, Rt)

    # Composed form: plain renewal core + a susceptible-depletion modifier.
    comp_step = ComposedRenewalStep(
        ConstantRenewalStep(rev_gen), (SusceptibleDepletion(N),))
    comp_init = _renewal_init_state(comp_step, 5.0, 0.1, length(rev_gen))
    comp = accumulate_scan(comp_step, comp_init, Rt)

    @test comp ≈ ref
    # Init states must also agree so the two are seeded identically.
    @test comp_init[1] ≈ ref_init[1]
    @test comp_init[2] ≈ ref_init[2]
end

@testitem "ComposedRenewalStep is ForwardDiff-differentiable" begin
    using ComposableTuringIDModels:
                                    ConstantRenewalStep, ComposedRenewalStep,
                                    SusceptibleDepletion,
                                    accumulate_scan, _renewal_init_state
    using ForwardDiff
    rev_gen = reverse([0.2, 0.3, 0.5])
    N = 1000.0
    comp_step = ComposedRenewalStep(
        ConstantRenewalStep(rev_gen), (SusceptibleDepletion(N),))
    # Differentiate the summed infection path with respect to the Rt path; the
    # composed recurrence must stay AD-friendly (no mutation of tracked state).
    f = function (Rt)
        init = _renewal_init_state(comp_step, 5.0, 0.1, length(rev_gen))
        sum(accumulate_scan(comp_step, init, Rt))
    end
    g = ForwardDiff.gradient(f, [1.6, 1.4, 1.2, 1.0])
    @test all(isfinite, g)
    @test length(g) == 4
end
