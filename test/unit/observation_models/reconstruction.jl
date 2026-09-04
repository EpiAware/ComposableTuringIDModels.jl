@testitem "an identity rebuild returns the same component" begin
    using ComposableTuringIDModels, Distributions
    using Accessors: Accessors

    # A wrapper whose constructor derives a field is rebuilt by re-deriving it
    # from the stored fields, so a no-op has to return the same model. Nothing
    # throws when it does not, so the check is on the concrete type: a second
    # wrapping shows up there.
    err = PoissonError()
    wrappers = [
        LatentDelay(err, [0.5, 0.3, 0.2]),
        LatentDelay(err, Gamma(2.0, 1.0); D = 5.0),
        Ascertainment(err, FixedIntercept(0.1)),
        Ascertainment(err, Normal(0.0, 0.1)),
        Ascertainment(err, FixedIntercept(0.1); latent_prefix = ""),
        Aggregate(err, [0, 0, 3]),
        RecordExpectedObs(err),
        TransformObservationModel(err),
        PrefixObservationModel(err, "a"),
        RightTruncate(err, [0.5, 1.0]),
        ReportTriangle(err, [0.6, 0.4]),
        Split((cases = err, deaths = NegativeBinomialError())),
        Split(err, [1.0 1.0]),
    ]
    for w in wrappers
        rebuilt = Accessors.modify(identity, w, Accessors.Properties())
        @test typeof(rebuilt) === typeof(w)
        for f in fieldnames(typeof(w))
            @test typeof(getfield(rebuilt, f)) === typeof(getfield(w, f))
        end
    end
end

@testitem "swapping a component leaves the rest of a wrapper alone" begin
    using ComposableTuringIDModels, Distributions
    using Accessors: Accessors, @set

    # `Ascertainment` stores its prior already wrapped in a `PrefixLatentModel`,
    # so a rebuild that re-applied the constructor would wrap it a second time.
    asc = Ascertainment(NormalError(), FixedIntercept(0.1))
    swapped = @set asc.model = PoissonError()
    @test swapped.model isa PoissonError
    @test typeof(getfield(swapped, :latent_model)) ===
        typeof(getfield(asc, :latent_model))
    @test swapped.latent_prefix == asc.latent_prefix

    # `LatentDelay` holds its PMF as it was given, so a rebuild leaves it alone.
    ld = LatentDelay(PoissonError(), [0.5, 0.3, 0.2])
    ld2 = @set ld.model = NegativeBinomialError()
    @test ld2.model isa NegativeBinomialError
    @test ld2.delay == ld.delay

    # `Aggregate` derives its presence mask from its window lengths, so a
    # rebuild must keep the two consistent.
    ag = Aggregate(PoissonError(), [0, 0, 3])
    ag2 = @set ag.model = NormalError()
    @test ag2.aggregation == ag.aggregation
    @test ag2.present == ag.present
end

@testitem "a rebuilt chain generates the same expected series" begin
    using ComposableTuringIDModels, Distributions, Random
    using Accessors: Accessors

    # The consequence a type check stands in for: a corrupted rebuild draws a
    # different prior for the doubly-prefixed effect, so the expected series it
    # generates differs from the original's.
    obs = LatentDelay(
        Ascertainment(PoissonError(), FixedIntercept(log(0.1))), [0.5, 0.3, 0.2]
    )
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal(1.0, 0.1)),
        obs
    )
    rebuilt = IDModel(
        model.infection_model,
        Accessors.modify(identity, obs, Accessors.Properties())
    )
    Random.seed!(42)
    a = as_turing_model(model, missing, 20)()
    Random.seed!(42)
    b = as_turing_model(rebuilt, missing, 20)()
    @test a.expected_y_t == b.expected_y_t
end

@testitem "a derived field is re-derived on both rebuild paths" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: rewrap, wrapped_models
    using Accessors: Accessors

    # A component that derives a field from the others at construction has to
    # answer the same down two paths, both of which go through
    # `ConstructionBase.constructorof`. Rebuilt from its own stored fields it
    # must not change, which is what `rewrap` and `swap` rely on. With a field
    # set to a user-facing value it must match the component built that way from
    # the start, which is what deriving one model from another relies on. A
    # component that re-derives on only one of them is silently wrong on the
    # other, and the two pull against each other: re-deriving is what the second
    # needs and what the first must survive.
    #
    # `==` on these structs falls back to `===`, which is false for any two
    # separately built arrays, so the comparison is structural.
    function structeq(a, b)
        (a isa Type || b isa Type) && return a === b
        typeof(a) === typeof(b) || return false
        a isa AbstractArray && return length(a) == length(b) &&
            all(structeq(x, y) for (x, y) in zip(a, b))
        a isa Number && return a == b
        a isa AbstractString && return a == b
        isbits(a) && return a == b
        nf = fieldcount(typeof(a))
        nf == 0 && return true
        return all(structeq(getfield(a, i), getfield(b, i)) for i in 1:nf)
    end

    err, err2 = PoissonError(), NegativeBinomialError()
    gp = truncated(Normal(5.0, 2.0), 0.0, 15.0)
    gp2 = truncated(Normal(3.0, 1.0), 0.0, 15.0)
    pmfs(p) = [p for _ in 1:4]
    ud(μ) = UncertainDelay(LogNormal, [Normal(μ, 0.4), Normal(0.4, 0.2)]; D = 10.0)

    # Each case is the component's user-facing constructor, the field names in
    # declaration order, and two sets of arguments to build it from.
    cases = [
        (
            "Ascertainment", (m, l, p) -> Ascertainment(m, l; latent_prefix = p),
            (:model, :latent_model, :latent_prefix),
            (err, FixedIntercept(0.1), "Ascertainment"),
            (err2, Intercept(Normal(0.0, 0.2)), "Reporting"),
        ),
        (
            "Ascertainment, bare prior and no prefix",
            (m, l, p) -> Ascertainment(m, l; latent_prefix = p),
            (:model, :latent_model, :latent_prefix),
            (err, Normal(0.0, 0.1), ""),
            (err2, Normal(1.0, 0.2), "Reporting"),
        ),
        (
            "Aggregate", (m, a) -> Aggregate(m, a), (:model, :aggregation),
            (err, [0, 0, 3]), (err2, [0, 7, 0]),
        ),
        (
            "LatentDelay from a PMF", (m, d) -> LatentDelay(m, d),
            (:model, :delay), (err, [0.5, 0.3, 0.2]), (err2, [0.1, 0.2, 0.7]),
        ),
        (
            "LatentDelay from a distribution", (m, d) -> LatentDelay(m, d),
            (:model, :delay), (err, gp), (err2, gp2),
        ),
        (
            "LatentDelay from a PMF per time", (m, d) -> LatentDelay(m, d),
            (:model, :delay),
            (err, pmfs([0.5, 0.3, 0.2])), (err2, pmfs([0.1, 0.2, 0.7])),
        ),
        (
            "LatentDelay from an UncertainDelay", (m, d) -> LatentDelay(m, d),
            (:model, :delay), (err, ud(1.5)), (err2, ud(2.5)),
        ),
        (
            "RightTruncate", (m, c) -> RightTruncate(m, c), (:model, :cdf_model),
            (err, [0.5, 1.0]), (err2, [0.2, 0.6, 1.0]),
        ),
        (
            "ReportTriangle", (m, d) -> ReportTriangle(m, d),
            (:error_model, :delay_model),
            (err, [0.6, 0.4]), (err2, [0.2, 0.3, 0.5]),
        ),
        (
            "UncertainDelay", (p, f) -> UncertainDelay(f, p; D = 10.0),
            (:params, :family),
            ([Normal(1.5, 0.4), Normal(0.4, 0.2)], LogNormal),
            ([RandomWalk(), Normal(0.5, 0.2)], LogNormal),
        ),
        (
            "ImportedCases", (r, t) -> ImportedCases(r; transformation = t),
            (:importation_rate, :transformation),
            (Normal(-1.0, 0.5), exp), (RandomWalk(), exp),
        ),
        (
            "Renewal", (r, g, m) -> Renewal(;
                rt = r, generation_time = g, initialisation = Normal(), mixing = m
            ), (:rt, :gen_int, :mixing),
            (RandomWalk(), [0.4, 0.6], Gravity([100.0, 200.0], [1.0 2.0; 2.0 1.0])),
            (Normal(), [0.2, 0.8], Gravity([50.0, 400.0], [1.0 3.0; 3.0 1.0])),
        ),
        (
            "Renewal with modifiers", (r, g, m) -> Renewal(;
                rt = r, generation_time = g, initialisation = Normal(), modifiers = m
            ), (:rt, :gen_int, :modifiers),
            (RandomWalk(), [0.4, 0.6], ()),
            (Normal(), [0.2, 0.8], (ImportedCases(Normal(-1.0, 0.5)),)),
        ),
        (
            "HierarchicalNormal", (m, s) -> HierarchicalNormal(; mean = m, std = s),
            (:mean, :std),
            (0.0, HalfNormal(0.1)), (1.0, HalfNormal(0.5)),
        ),
        (
            "CombineLatentModels", (m, p) -> CombineLatentModels(m, p),
            (:models, :prefixes),
            ([Normal(), Normal()], ["a", "b"]),
            ([RandomWalk(), Normal(1.0, 2.0)], ["c", "d"]),
        ),
        (
            "ConcatLatentModels", (m) -> ConcatLatentModels(m), (:models,),
            ([Normal(), Normal()],), ([RandomWalk(), Normal(1.0, 2.0)],),
        ),
    ]

    for (name, build, fields, argsA, argsB) in cases
        @testset "$name" begin
            a = build(argsA...)

            # Rebuilt from its own stored fields, unchanged.
            @test structeq(
                Accessors.modify(identity, a, Accessors.Properties()), a
            )

            # One field at a time set to the user-facing value, matching the
            # component built with that value from the start.
            for (i, f) in enumerate(fields)
                argsC = ntuple(
                    j -> j == i ? argsB[j] : argsA[j], length(argsA)
                )
                set = Accessors.set(a, Accessors.PropertyLens{f}(), argsB[i])
                @test structeq(set, build(argsC...))
            end

            # An observation model is also rebuilt by `rewrap`, which is the
            # same path with the wrapped models substituted.
            if a isa AbstractObservationModel
                wrapped = wrapped_models(a)
                @test structeq(rewrap(a, wrapped), a)
                if !isempty(wrapped)
                    swapped = rewrap(a, (NormalError(),))
                    @test only(wrapped_models(swapped)) isa NormalError
                end
            end
        end
    end
end

@testitem "every component pointing constructorof at a rebuild is covered" begin
    using ComposableTuringIDModels
    using ComposableTuringIDModels: ConstructionBase

    # The case list above is written out, so a component that later derives a
    # field could join the package without joining it. A component that derives
    # one declares its own `ConstructionBase.constructorof`, so requiring every
    # declaration to be named here fails until the case is added.
    covered = Set(
        [
            :Ascertainment, :Aggregate, :HierarchicalNormal, :ImportedCases,
            :Renewal, :UncertainDelay,
        ]
    )
    # A `constructorof` method's second argument is `Type{<:T}`, so `T` is the
    # upper bound of the unwrapped type variable.
    subject(m) = nameof(
        Base.unwrap_unionall(
            Base.unwrap_unionall(m.sig.parameters[2]).parameters[1].ub
        )
    )
    declared = Set(
        subject(m) for m in methods(ConstructionBase.constructorof)
            if parentmodule(m) === ComposableTuringIDModels
    )
    @test declared == covered
end

@testitem "a chain derived with @set samples the same variables" begin
    using ComposableTuringIDModels, Distributions
    using Accessors: @set
    using DynamicPPL: VarInfo

    # What losing the prefix wrapper costs: the ascertainment's prior is
    # namespaced by it, so a chain derived from another by setting the prior
    # sampled a differently named variable from the chain it appears to copy.
    # Nothing throws, and the difference shows only in the variable names, so
    # they are what is pinned here.
    chain(prior, delay) = LatentDelay(
        Ascertainment(
            NegativeBinomialError(; cluster_factor = HalfNormal(0.1)), prior
        ), delay
    )
    infection = DirectInfections(;
        Z = RandomWalk(), initialisation = Normal(1.0, 0.1)
    )
    varnames(obs) = keys(VarInfo(as_turing_model(IDModel(infection, obs), missing, 30)))

    cases = chain(FixedIntercept(log(0.6)), LogNormal(1.6, 0.5))
    built = chain(Intercept(Normal(log(0.015), 0.25)), LogNormal(2.8, 0.4))

    derived = @set cases.model.latent_model = Intercept(Normal(log(0.015), 0.25))
    derived = @set derived.delay = LogNormal(2.8, 0.4)

    @test varnames(derived) == varnames(built)
    @test Symbol("Ascertainment.intercept") in Symbol.(varnames(built))
    @test derived.model.latent_model isa PrefixLatentModel
    @test derived.delay == built.delay
end
