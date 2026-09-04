# Rebuilding a component from its fields is what `Accessors.@set`, `rewrap` and
# `swap` all go through, so it is checked here for every component at once
# rather than per component.

@testitem "every component rebuilds to itself and takes a set field" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: equal_dimensions
    using OrdinaryDiffEq: AutoVern7, Rodas5P
    using Accessors: Accessors

    # Two properties, which together are the whole contract for a component
    # that transforms or derives a field at construction.
    #
    # 1. Rebuilt from its own stored fields it does not change. `rewrap` and
    #    `swap` replace one field and rebuild, so anything else the component
    #    holds has to survive the trip.
    # 2. With a field set to a user-facing value it equals the component built
    #    with that value from the start. Deriving one model from another relies
    #    on this, and it is where a stale derived field shows up.
    #
    # They pull against each other: the second wants the derivation re-applied,
    # the first wants that to be a no-op. Both hold when the derivation is
    # idempotent and lives in the constructor taking the fields in declaration
    # order, which is the one `ConstructionBase.constructorof` calls by default.
    # A component that satisfies neither needs no machinery, only that shape.

    # `==` on these structs falls back to `===`, which is false for any two
    # separately built arrays, so the comparison has to be structural.
    function structeq(a, b, depth = 0)
        (a isa Type || b isa Type) && return a === b
        typeof(a) === typeof(b) || return false
        a === b && return true
        # An `ODEProblem` holds a solver graph deep enough to exhaust the stack,
        # so below the cap the type match above is the whole comparison.
        depth > 8 && return true
        a isa AbstractArray && return length(a) == length(b) &&
            all(structeq(x, y, depth + 1) for (x, y) in zip(a, b))
        a isa Number && return a == b
        a isa AbstractString && return a == b
        isbits(a) && return a == b
        nf = fieldcount(typeof(a))
        nf == 0 && return true
        return all(
            structeq(getfield(a, i), getfield(b, i), depth + 1) for i in 1:nf
        )
    end

    err, err2 = PoissonError(), NegativeBinomialError()
    d1 = truncated(Normal(0.0, 0.05), 0, 1)
    t1 = truncated(Normal(0.0, 0.05), -1, 1)
    gp, gp2 = truncated(Normal(5.0, 2.0), 0.0, 15.0), truncated(Normal(3.0, 1.0), 0.0, 15.0)
    ud(mu) = UncertainDelay(LogNormal, [Normal(mu, 0.4), Normal(0.4, 0.2)]; D = 10.0)
    di, di2 = DirectInfections(; Z = Normal()), DirectInfections(; Z = RandomWalk())
    sir = SIRParams(;
        tspan = (0.0, 30.0), infectiousness = Normal(0.3, 0.05),
        recovery_rate = Normal(0.1, 0.01), initial_prop_infected = Beta(1, 99)
    )
    seir = SEIRParams(;
        tspan = (0.0, 30.0), infectiousness = Normal(0.3, 0.05),
        incubation_rate = Normal(0.2, 0.01), recovery_rate = Normal(0.1, 0.01),
        initial_prop_infected = Beta(1, 99)
    )

    CASES = [
        ("PoissonError", () -> PoissonError(), (), (), ()),
        ("BinomialError", () -> BinomialError(), (), (), ()),
        ("Null", () -> Null(), (), (), ()),
        (
            "NegativeBinomialError", c -> NegativeBinomialError(; cluster_factor = c),
            (:cluster_factor,), (HalfNormal(0.1),), (HalfNormal(0.5),),
        ),
        ("NormalError", s -> NormalError(s), (:std,), (HalfNormal(0.1),), (HalfNormal(0.5),)),
        ("FixedIntercept", i -> FixedIntercept(i), (:intercept,), (0.1,), (0.5,)),
        ("Intercept", i -> Intercept(i), (:intercept,), (Normal(),), (Normal(1.0, 2.0),)),
        ("IID", e -> IID(e), (:ϵ_t,), (Normal(),), (Normal(1.0, 2.0),)),
        (
            "AR", (d, i, e) -> AR(; damp = d, init = i, ϵ_t = e),
            (:damp, :init, :ϵ_t), ([d1], Normal(), HierarchicalNormal()),
            (fill(d1, 2), Normal(1.0, 2.0), RandomWalk()),
        ),
        (
            "MA", (t, e) -> MA(; θ = t, ϵ_t = e), (:θ, :ϵ_t),
            ([t1], HierarchicalNormal()), (fill(t1, 2), RandomWalk()),
        ),
        (
            "RandomWalk", (i, e) -> RandomWalk(; init = i, ϵ_t = e), (:init, :ϵ_t),
            (Normal(), HierarchicalNormal()), (Normal(1.0, 2.0), IID(Normal())),
        ),
        (
            "HierarchicalNormal", (m, s) -> HierarchicalNormal(; mean = m, std = s),
            (:mean, :std), (0.0, HalfNormal(0.1)), (1.0, HalfNormal(0.5)),
        ),
        (
            "HilbertSpaceGP", (l, s, m, c, k) -> HilbertSpaceGP(;
                length_scale = l, marginal_std = s, m = m, c = c, kernel = k
            ),
            (:length_scale, :marginal_std, :m, :c, :kernel),
            (truncated(Normal(0.0, 0.4), 0.01, Inf), HalfNormal(0.1), 20, 1.5, SqExponentialKernel()),
            (truncated(Normal(0.0, 0.8), 0.01, Inf), HalfNormal(0.5), 30, 2.0, Matern32Kernel()),
        ),
        (
            "ExactGP", (l, s, k, j) -> ExactGP(;
                length_scale = l, marginal_std = s, kernel = k, jitter = j
            ),
            (:length_scale, :marginal_std, :kernel, :jitter),
            (truncated(Normal(0.0, 0.4), 0.01, Inf), HalfNormal(0.1), SqExponentialKernel(), 1.0e-6),
            (truncated(Normal(0.0, 0.8), 0.01, Inf), HalfNormal(0.5), Matern52Kernel(), 1.0e-5),
        ),
        (
            "DiffLatentModel", (m, i) -> DiffLatentModel(; model = m, init = i),
            (:model, :init), (Normal(), [Normal()]), (RandomWalk(), [Normal(), Normal()]),
        ),
        (
            "TransformLatentModel", (m, t) -> TransformLatentModel(m, t),
            (:model, :transform), (Normal(), exp), (RandomWalk(), log1p),
        ),
        (
            "PrefixLatentModel", (m, p) -> PrefixLatentModel(m, p), (:model, :prefix),
            (Normal(), "a"), (RandomWalk(), "b"),
        ),
        (
            "RecordExpectedLatent", m -> RecordExpectedLatent(m), (:model,),
            (Normal(),), (RandomWalk(),),
        ),
        (
            "CombineLatentModels", (m, p) -> CombineLatentModels(m, p), (:models, :prefixes),
            ([Normal(), Normal()], ["a", "b"]), ([RandomWalk(), Normal(1.0, 2.0)], ["c", "d"]),
        ),
        (
            "ConcatLatentModels", (m, a, p) -> ConcatLatentModels(m, a; prefixes = p),
            (:models, :dimension_adaptor, :prefixes),
            ([Normal(), Normal()], equal_dimensions, ["a", "b"]),
            ([RandomWalk(), Normal(1.0, 2.0)], equal_dimensions, ["c", "d"]),
        ),
        (
            "Hierarchy", (m, a) -> Hierarchy(; mean = m, across = a), (:mean, :across),
            (Normal(), IID(Normal())), (Normal(1.0, 2.0), IID(Normal(1.0, 2.0))),
        ),
        (
            "Stratify", (s, a, c) -> Stratify(s, a; combine = c), (:shared, :across, :combine),
            (Normal(), Normal(), +), (RandomWalk(), Hierarchy(), *),
        ),
        ("Replicate", m -> Replicate(m), (:model,), (Normal(),), (RandomWalk(),)),
        (
            "BroadcastLatentModel", (m, p, r) -> BroadcastLatentModel(m; period = p, broadcast_rule = r),
            (:model, :period, :broadcast_rule),
            (Normal(), 2, RepeatEach()), (RandomWalk(), 3, RepeatBlock()),
        ),
        (
            "DirectInfections", (z, t, i) -> DirectInfections(; Z = z, transformation = t, initialisation = i),
            (:Z, :transformation, :initialisation),
            (Normal(), exp, Normal()), (RandomWalk(), identity, Normal(1.0, 2.0)),
        ),
        (
            "ExpGrowthRate", (r, t, i) -> ExpGrowthRate(; rt = r, transformation = t, initialisation = i),
            (:rt, :transformation, :initialisation),
            (Normal(), exp, Normal()), (RandomWalk(), identity, Normal(1.0, 2.0)),
        ),
        ("SeedingPath", m -> SeedingPath(m), (:model,), (Normal(),), (RandomWalk(),)),
        (
            "CombineInfections", (m, n) -> CombineInfections(m, n), (:models, :names),
            ([di, di], ["a", "b"]), ([di2, di], ["c", "d"]),
        ),
        (
            "Gravity", (p, d, a, b, g, w) -> Gravity(p, d; α = a, β = b, γ = g, within = w),
            (:pop, :dist, :α, :β, :γ, :within),
            ([100.0, 200.0], [1.0 2.0; 2.0 1.0], Normal(0, 0.5), Normal(0, 0.5), truncated(Normal(2, 0.5), 0, Inf), 1.0),
            ([50.0, 400.0], [1.0 3.0; 3.0 1.0], Normal(1, 0.5), Normal(1, 0.5), truncated(Normal(3, 0.5), 0, Inf), 2.0),
        ),
        (
            "ImportedCases", (r, t) -> ImportedCases(r; transformation = t),
            (:importation_rate, :transformation),
            (Normal(-1.0, 0.5), exp), (RandomWalk(), identity),
        ),
        (
            "Renewal", (g, tr, r, i, mx, mo) -> Renewal(;
                generation_time = g, transformation = tr, rt = r,
                initialisation = i, mixing = mx, modifiers = mo
            ),
            (:gen_int, :transformation, :rt, :initialisation, :mixing, :modifiers),
            ([0.4, 0.6], exp, RandomWalk(), Normal(), Gravity([100.0, 200.0], [1.0 2.0; 2.0 1.0]), ()),
            (
                [0.2, 0.8], identity, Normal(), Normal(1.0, 2.0),
                Gravity([50.0, 400.0], [1.0 3.0; 3.0 1.0]), (ImportedCases(Normal(-1.0, 0.5)),),
            ),
        ),
        (
            "UncertainDelay", (p, f, d, dd) -> UncertainDelay(f, p; D = d, Δd = dd),
            (:params, :family, :D, :Δd),
            ([Normal(1.5, 0.4), Normal(0.4, 0.2)], LogNormal, 10.0, 1.0),
            ([RandomWalk(), Normal(0.5, 0.2)], Gamma, 20.0, 2.0),
        ),
        (
            "LatentDelay", (m, d) -> LatentDelay(m, d), (:model, :delay),
            (err, [0.5, 0.3, 0.2]), (err2, [0.1, 0.2, 0.7]),
        ),
        (
            "LatentDelay, distribution", (m, d) -> LatentDelay(m, d), (:model, :delay),
            (err, gp), (err2, gp2),
        ),
        (
            "LatentDelay, PMF per time", (m, d) -> LatentDelay(m, d), (:model, :delay),
            (err, [[0.5, 0.3, 0.2] for _ in 1:4]), (err2, [[0.1, 0.2, 0.7] for _ in 1:4]),
        ),
        (
            "LatentDelay, UncertainDelay", (m, d) -> LatentDelay(m, d), (:model, :delay),
            (err, ud(1.5)), (err2, ud(2.5)),
        ),
        (
            "Ascertainment", (m, l, p) -> Ascertainment(m, l; latent_prefix = p),
            (:model, :latent_model, :latent_prefix),
            (err, FixedIntercept(0.1), "Ascertainment"), (err2, Intercept(Normal()), "Reporting"),
        ),
        (
            "Ascertainment, bare prior and no prefix", (m, l, p) -> Ascertainment(m, l; latent_prefix = p),
            (:model, :latent_model, :latent_prefix),
            (err, Normal(0.0, 0.1), ""), (err2, Normal(1.0, 0.2), "Reporting"),
        ),
        (
            "Ascertainment, a prefixed prior re-prefixed", (m, l, p) -> Ascertainment(m, l; latent_prefix = p),
            (:model, :latent_model, :latent_prefix),
            (err, PrefixLatentModel(Normal(0.0, 0.1), "own"), "Ascertainment"),
            (err2, PrefixLatentModel(Normal(1.0, 0.2), "own"), ""),
        ),
        (
            "Aggregate", (m, a) -> Aggregate(m, a), (:model, :aggregation),
            (err, [0, 0, 3]), (err2, [0, 7, 0]),
        ),
        (
            "RightTruncate", (m, c) -> RightTruncate(m, c), (:model, :cdf_model),
            (err, [0.5, 1.0]), (err2, [0.2, 0.6, 1.0]),
        ),
        (
            "ReportTriangle", (m, d) -> ReportTriangle(m, d), (:error_model, :delay_model),
            (err, [0.6, 0.4]), (err2, [0.2, 0.3, 0.5]),
        ),
        ("ReportingCDF", c -> ReportingCDF(c), (:cdf,), ([0.5, 1.0],), ([0.2, 0.6, 1.0],)),
        ("ReportingPMF", p -> ReportingPMF(p), (:pmf,), ([0.5, 0.5],), ([0.2, 0.3, 0.5],)),
        (
            "PrefixObservationModel", (m, p) -> PrefixObservationModel(m, p), (:model, :prefix),
            (err, "a"), (err2, "b"),
        ),
        ("RecordExpectedObs", m -> RecordExpectedObs(m), (:model,), (err,), (err2,)),
        (
            "TransformObservationModel", (m, t) -> TransformObservationModel(m, t),
            (:model, :transform), (err, exp), (err2, log1p),
        ),
        (
            "Split", (s, m) -> Split(s, m), (:streams, :map),
            (err, [1.0 1.0]), (err2, [0.5 0.5]),
        ),
        (
            "IDModel", (i, o) -> IDModel(i, o), (:infection_model, :observation_model),
            (di, err), (di2, err2),
        ),
        (
            "SIRParams", (b, g, i) -> SIRParams(;
                tspan = (0.0, 30.0),
                infectiousness = b, recovery_rate = g, initial_prop_infected = i
            ),
            (:infectiousness, :recovery_rate, :initial_prop_infected),
            (Normal(0.3, 0.05), Normal(0.1, 0.01), Beta(1, 99)),
            (Normal(0.5, 0.05), Normal(0.2, 0.01), Beta(2, 98)),
        ),
        (
            "SEIRParams", (b, a, g, i) -> SEIRParams(;
                tspan = (0.0, 30.0),
                infectiousness = b, incubation_rate = a, recovery_rate = g,
                initial_prop_infected = i
            ),
            (:infectiousness, :incubation_rate, :recovery_rate, :initial_prop_infected),
            (Normal(0.3, 0.05), Normal(0.2, 0.01), Normal(0.1, 0.01), Beta(1, 99)),
            (Normal(0.5, 0.05), Normal(0.3, 0.01), Normal(0.2, 0.01), Beta(2, 98)),
        ),
        (
            "ODEProcess", (p, s, f, o) -> ODEProcess(;
                params = p, solver = s, sol2infs = f, solver_options = o
            ),
            (:params, :solver, :sol2infs, :solver_options),
            (sir, AutoVern7(Rodas5P()), sol -> sol[2, :], Dict(:saveat => 1.0)),
            (seir, Rodas5P(), sol -> sol[3, :], Dict(:saveat => 0.5)),
        ),
    ]
    for (name, build, fields, argsA, argsB) in CASES
        @testset "$name" begin
            a = build(argsA...)
            @test structeq(
                Accessors.modify(identity, a, Accessors.Properties()), a
            )
            for (i, f) in enumerate(fields)
                argsC = ntuple(
                    j -> j == i ? argsB[j] : argsA[j], length(argsA)
                )
                set = Accessors.set(a, Accessors.PropertyLens{f}(), argsB[i])
                @test structeq(set, build(argsC...))
            end
        end
    end
end

@testitem "every concrete component is covered by the rebuild cases" begin
    using ComposableTuringIDModels
    using ComposableTuringIDModels: AbstractComposableModel
    using InteractiveUtils: subtypes

    # The cases above are written out, so without this a component added later
    # could carry a non-idempotent construction-time transform and never be
    # asked about it. Requiring every concrete component to appear makes the
    # omission fail rather than pass quietly.
    # A test item elsewhere in the suite can define its own component, and
    # `subtypes` sees every one loaded in the process, so the walk is restricted
    # to the package's own types.
    function leaves(T, acc = Set{Any}())
        for S in subtypes(T)
            isabstracttype(S) && (leaves(S, acc); continue)
            parentmodule(S) === ComposableTuringIDModels && push!(acc, S)
        end
        return acc
    end

    # `CatalystODEParams` builds only from a `ReactionSystem`, which needs the
    # Catalyst extension; `test/unit/ode/catalyst_ext.jl` checks it there.
    skipped = Set([:CatalystODEParams])
    covered = Set(
        [
            :AR, :Aggregate, :Ascertainment, :BinomialError,
            :BroadcastLatentModel, :CombineInfections, :CombineLatentModels,
            :ConcatLatentModels, :DiffLatentModel, :DirectInfections, :ExactGP,
            :ExpGrowthRate, :FixedIntercept, :Gravity, :HierarchicalNormal,
            :Hierarchy, :HilbertSpaceGP, :IDModel, :IID, :Intercept,
            :LatentDelay, :MA, :NegativeBinomialError, :NormalError, :Null,
            :ODEProcess, :PoissonError, :PrefixLatentModel,
            :PrefixObservationModel, :RandomWalk, :RecordExpectedLatent,
            :RecordExpectedObs, :Renewal, :Replicate, :ReportTriangle,
            :ReportingCDF, :ReportingPMF, :RightTruncate, :SEIRParams,
            :SIRParams, :SeedingPath, :Split, :Stratify, :TransformLatentModel,
            :TransformObservationModel, :UncertainDelay,
        ]
    )
    concrete = Set(nameof.(leaves(AbstractComposableModel)))
    @test setdiff(concrete, covered, skipped) == Set{Symbol}()
    # And nothing named here has since been removed or renamed.
    @test setdiff(covered, concrete) == Set{Symbol}()
end

@testitem "a slot naming variables refuses a rebuild it cannot name" begin
    using ComposableTuringIDModels, Distributions
    using Accessors: Accessors

    # A prefix vector is not derived from the models it names, it is chosen.
    # Re-deriving it on a rebuild would rename the chain's variables to defaults
    # without saying so, so the length mismatch is raised instead.
    combine = CombineLatentModels([Normal(), Normal()], ["a", "b"])
    @test_throws AssertionError Accessors.set(
        combine, Accessors.PropertyLens{:models}(),
        [Normal(), Normal(), Normal()]
    )

    concat = ConcatLatentModels([Normal(), Normal()])
    @test_throws AssertionError Accessors.set(
        concat, Accessors.PropertyLens{:models}(),
        [Normal(), Normal(), Normal()]
    )

    # Setting both together is how the count changes, and the derived
    # `no_models` follows the models rather than the count it was built with.
    grown = ComposableTuringIDModels.ConstructionBase.setproperties(
        concat, (models = [Normal(), Normal(), Normal()], prefixes = ["a", "b", "c"])
    )
    @test grown.no_models == 3
end
