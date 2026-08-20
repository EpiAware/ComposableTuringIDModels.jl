@testitem "wrapped_models reports what each component wraps" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: wrapped_models

    # An error model consumes the expected series itself, so it wraps nothing.
    @test wrapped_models(PoissonError()) == ()

    # A modifier wraps exactly one observation model, whatever the slot is
    # called: `ReportTriangle` names its inner model `error_model`.
    inner = PoissonError()
    @test wrapped_models(LatentDelay(inner, [0.5, 0.3, 0.2])) == (inner,)
    @test wrapped_models(Ascertainment(inner, FixedIntercept(0.1))) == (inner,)
    @test wrapped_models(Aggregate(inner, [0, 0, 3])) == (inner,)
    @test wrapped_models(RecordExpectedObs(inner)) == (inner,)
    @test wrapped_models(TransformObservationModel(inner)) == (inner,)
    @test wrapped_models(PrefixObservationModel(inner, "a")) == (inner,)
    @test wrapped_models(RightTruncate(inner, [0.5, 1.0])) == (inner,)
    @test wrapped_models(ReportTriangle(inner, [0.6, 0.4])) == (inner,)

    # A `Split` branches rather than nests: every stream is wrapped at once.
    cases = PoissonError()
    deaths = NegativeBinomialError()
    @test wrapped_models(Split((cases = cases, deaths = deaths))) ==
        (cases, deaths)
    # A strata template is one model replicated across data-named streams.
    @test wrapped_models(Split(cases, [1.0 1.0])) == (cases,)
end

@testitem "observation_components walks a chain and its branches" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: observation_components

    err = PoissonError()
    obs = LatentDelay(Ascertainment(err, FixedIntercept(0.1)), [0.5, 0.3, 0.2])
    parts = observation_components(obs)
    # Root first, then inwards.
    @test parts[1] === obs
    @test parts[2] === obs.model
    @test parts[3] === err
    @test length(parts) == 3

    # Find by type instead of writing a field path.
    delays = filter(x -> x isa LatentDelay, observation_components(obs))
    @test length(delays) == 1
    @test only(delays).delay == [0.2, 0.3, 0.5]

    # The lead-in a downstream package hand-rolls is a sum over the walk.
    lead_in = sum(
        length(d.delay) - 1 for d in delays; init = 0
    )
    @test lead_in == 2

    # Branches are followed, each stream in order.
    split = Split(
        (
            cases = LatentDelay(PoissonError(), [0.5, 0.5]),
            deaths = LatentDelay(NegativeBinomialError(), [0.2, 0.3, 0.5]),
        )
    )
    types = map(typeof, observation_components(split))
    @test count(t -> t <: LatentDelay, types) == 2
    @test count(t -> t <: AbstractObservationErrorModel, types) == 2
    @test first(observation_components(split)) === split

    # A composed model hands over its observation chain.
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()), obs
    )
    @test observation_components(model) == parts
end

@testitem "rewrap rebuilds a component around new wrapped models" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: wrapped_models, rewrap

    replacement = NegativeBinomialError()
    ld = LatentDelay(PoissonError(), [0.5, 0.3, 0.2])
    swapped = rewrap(ld, (replacement,))
    @test swapped isa LatentDelay
    @test only(wrapped_models(swapped)) === replacement
    # The wrapper's own state is carried across untouched: the stored PMF is
    # already reversed, and rebuilding through the public constructor would
    # reverse it a second time.
    @test swapped.delay == ld.delay

    # `Ascertainment` stores its prior already wrapped in a `PrefixLatentModel`,
    # so a rebuild must not prefix it again.
    asc = Ascertainment(PoissonError(), FixedIntercept(0.1))
    asc2 = rewrap(asc, (replacement,))
    @test typeof(getfield(asc2, :latent_model)) ===
        typeof(getfield(asc, :latent_model))
    @test asc2.latent_prefix == asc.latent_prefix

    # An `Aggregate` derives its presence mask, so it too has to be rebuilt
    # through the constructor that recomputes rather than re-supplies it.
    ag = Aggregate(PoissonError(), [0, 0, 3])
    ag2 = rewrap(ag, (replacement,))
    @test ag2.aggregation == ag.aggregation
    @test ag2.present == ag.present

    # A `Split` keeps its stream names and weight map.
    split = Split((cases = PoissonError(), deaths = PoissonError()), [1.0; 1.0;;])
    split2 = rewrap(split, (replacement, replacement))
    @test keys(split2.streams) == (:cases, :deaths)
    @test split2.names == split.names
    @test split2.map == split.map
end

@testitem "every observation component rebuilds to itself" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: wrapped_models, rewrap,
        observation_components

    # An identity rebuild must return an identical component. Any wrapper whose
    # constructor transforms an argument before storing it fails this without
    # anything being thrown, so compare structurally rather than by `isa`.
    function same(a, b)
        typeof(a) === typeof(b) || return false
        return all(
            fieldnames(typeof(a))
        ) do f
            x, y = getfield(a, f), getfield(b, f)
            x isa Function ? x === y :
                (
                    x isa Union{AbstractComposableModel, NamedTuple} ? same(x, y) :
                    x == y
                )
        end
    end
    function same(a::NamedTuple, b::NamedTuple)
        keys(a) == keys(b) || return false
        return all(k -> same(a[k], b[k]), keys(a))
    end

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
        @test same(rewrap(w, wrapped_models(w)), w)
    end

    # Every component in the package's own chains is reachable and rebuildable.
    for w in wrappers
        @test all(
            c -> same(rewrap(c, wrapped_models(c)), c), observation_components(w)
        )
    end
end

@testitem "traversal and rewrap compose into a type-keyed swap" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: wrapped_models, rewrap,
        observation_components

    # The update case, written in terms of the two accessors alone.
    swap(f, m) = f(rewrap(m, map(x -> swap(f, x), wrapped_models(m))))

    obs = LatentDelay(
        Ascertainment(
            Split(
                (
                    cases = PoissonError(),
                    deaths = LatentDelay(PoissonError(), [0.2, 0.8]),
                )
            ),
            FixedIntercept(0.1)
        ),
        [0.5, 0.3, 0.2]
    )
    swapped = swap(x -> x isa PoissonError ? NegativeBinomialError() : x, obs)

    parts = observation_components(swapped)
    @test !any(x -> x isa PoissonError, parts)
    @test count(x -> x isa NegativeBinomialError, parts) == 2
    # Every wrapper survives the rebuild, in the same order and shape.
    @test map(nameof ∘ typeof, parts) ==
        map(nameof ∘ typeof, observation_components(obs))
    # The original is untouched.
    @test count(x -> x isa PoissonError, observation_components(obs)) == 2
end
