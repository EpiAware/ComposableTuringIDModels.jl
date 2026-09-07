@testitem "observation_lead_in sums nested LatentDelay lead-ins" begin
    using ComposableTuringIDModels, Distributions

    # A bare error family consumes nothing.
    @test observation_lead_in(PoissonError()) == 0

    # One delay costs `length(pmf) - 1`.
    d1 = fill(1 / 15, 15)
    @test observation_lead_in(LatentDelay(PoissonError(), d1)) == 14

    # Nested delays sum.
    d2 = fill(1 / 30, 30)
    chain = LatentDelay(LatentDelay(PoissonError(), d1), d2)
    @test observation_lead_in(chain) == 14 + 29
end

@testitem "observation_lead_in walks every observation nesting shape" begin
    using ComposableTuringIDModels, Distributions

    pmf = fill(1 / 5, 5)                       # lead-in 4
    delayed = LatentDelay(PoissonError(), pmf)

    # Modifiers that keep the series length pass the lead-in through.
    @test observation_lead_in(Ascertainment(delayed, FixedIntercept(0.0))) == 4
    @test observation_lead_in(RightTruncate(delayed, fill(1.0, 10))) == 4
    @test observation_lead_in(TransformObservationModel(delayed)) == 4
    @test observation_lead_in(PrefixObservationModel(delayed, "a")) == 4
    @test observation_lead_in(RecordExpectedObs(delayed)) == 4

    # A ReportTriangle's delay is a nowcasting kernel, not a series-shortening
    # convolution, so it contributes nothing of its own.
    @test observation_lead_in(ReportTriangle(PoissonError(), fill(1 / 3, 3))) == 0

    # An aggregation applied to a delayed series is on the time axis, so it
    # passes the lead-in through. A delay applied INSIDE one consumes reporting
    # windows drawn from the series the aggregation already has, so the series
    # does not lengthen and the lead-in in time points is zero.
    weekly = [0, 0, 0, 0, 0, 0, 7]
    @test observation_lead_in(LatentDelay(Aggregate(PoissonError(), weekly), pmf)) == 4
    @test observation_lead_in(Aggregate(PoissonError(), weekly)) == 0
    @test observation_lead_in(Aggregate(delayed, weekly)) == 0

    # Streams of a Split run in parallel: the shared lead-in, not their sum.
    split = Split((cases = delayed, deaths = LatentDelay(PoissonError(), pmf)))
    @test observation_lead_in(split) == 4
    # A Split placed under a delay adds that delay's lead-in.
    @test observation_lead_in(LatentDelay(split, fill(1 / 3, 3))) == 4 + 2
    # A strata-template Split reports its template's lead-in.
    @test observation_lead_in(Split(delayed, [1.0 1.0])) == 4

    # Streams that consume different amounts report one lead-in each, and a
    # surrounding delay adds to every stream.
    mismatched = Split(
        (
            cases = delayed, deaths = LatentDelay(PoissonError(), fill(1 / 9, 9)),
        )
    )
    @test observation_lead_in(mismatched) == (cases = 4, deaths = 8)
    @test observation_lead_in(LatentDelay(mismatched, fill(1 / 3, 3))) ==
        (cases = 6, deaths = 10)
    @test observation_lead_in(Ascertainment(mismatched, FixedIntercept(0.0))) ==
        (cases = 4, deaths = 8)

    # A stacked chain sums each delay it passes through.
    stacked = LatentDelay(
        Ascertainment(
            LatentDelay(RecordExpectedObs(NegativeBinomialError()), fill(1 / 4, 4)),
            FixedIntercept(0.0)
        ),
        pmf
    )
    @test observation_lead_in(stacked) == 4 + 3

    # An IDModel and an IDProblem delegate to their observation model.
    infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    @test observation_lead_in(IDModel(infection, stacked)) == 7
    problem = IDProblem(
        infection = infection, observation_model = stacked, tspan = (1, 30)
    )
    @test observation_lead_in(problem) == 7
end

@testitem "observation_lead_in reads the traversal seam" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: wrapped_models, observation_components

    # The walk is the traversal API's, not a second one: a chain's lead-in is
    # the sum over the `LatentDelay`s the traversal reports.
    obs = LatentDelay(
        Ascertainment(LatentDelay(PoissonError(), fill(1 / 4, 4)), FixedIntercept(0.0)),
        fill(1 / 5, 5)
    )
    delays = filter(x -> x isa LatentDelay, observation_components(obs))
    @test observation_lead_in(obs) ==
        sum(length(d.delay) - 1 for d in delays; init = 0)

    # A modifier is seen through `wrapped_models`, so one defined outside the
    # package needs no lead-in method of its own.
    @test only(wrapped_models(obs)) === obs.model
end

@testitem "observation_lead_in covers uncertain and time-varying delays" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1064)

    # An `UncertainDelay` has a fixed horizon, so its pmf length is known
    # before the parameters are drawn.
    uncertain = UncertainDelay(
        LogNormal, [Normal(1.5, 0.4), truncated(Normal(0.4, 0.2), 0, Inf)]; D = 8.0
    )
    obs_uncertain = LatentDelay(NegativeBinomialError(), uncertain)
    @test observation_lead_in(obs_uncertain) == 7
    sim = as_turing_model(obs_uncertain, missing, fill(100.0, 40))()
    @test length(sim.expected) == 40 - observation_lead_in(obs_uncertain)

    # A deterministic per-time sequence of pmfs: all pmfs share a length.
    n = 20
    pmfs = [fill(1 / 4, 4) for _ in 1:n]
    obs_tv = LatentDelay(PoissonError(), pmfs)
    @test observation_lead_in(obs_tv) == 3
    sim_tv = as_turing_model(obs_tv, missing, fill(100.0, n))()
    @test length(sim_tv.expected) == n - observation_lead_in(obs_tv)

    # A delay whose length cannot be read off the component is an error, not a
    # silent zero.
    opaque = LatentDelay(PoissonError(), RandomWalk())
    @test_throws ArgumentError observation_lead_in(opaque)
end

@testitem "observation_lead_in predicts the expected-series length" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1065)

    n = 40
    pmf = fill(1 / 6, 6)
    chains = [
        LatentDelay(PoissonError(), pmf),
        LatentDelay(Ascertainment(PoissonError(), FixedIntercept(0.0)), pmf),
        LatentDelay(LatentDelay(PoissonError(), pmf), fill(1 / 4, 4)),
        Ascertainment(LatentDelay(PoissonError(), pmf), FixedIntercept(0.0)),
        TransformObservationModel(LatentDelay(PoissonError(), pmf); transform = identity),
        RightTruncate(LatentDelay(PoissonError(), pmf), fill(1.0, n)),
    ]
    for obs in chains
        sim = as_turing_model(obs, missing, fill(100.0, n))()
        @test length(sim.expected) == n - observation_lead_in(obs)
    end
end

@testitem "as_turing_model covers a chain's lead-in itself" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: VarInfo, logjoint
    Random.seed!(1066)

    d1 = fill(1 / 15, 15)
    d2 = fill(1 / 30, 30)
    obs = LatentDelay(LatentDelay(PoissonError(), d1), d2)
    # A tight random walk keeps the prior draw's expected series near 100, so
    # the log-joint stays well scaled and a perturbation of a scored entry is
    # not lost to floating-point rounding.
    infection = Renewal(;
        generation_time = [0.2, 0.3, 0.3, 0.2],
        rt = RandomWalk(; init = Normal(0.0, 0.05), ϵ_t = Normal(0.0, 0.01)),
        initialisation = Normal(log(100.0), 0.1)
    )
    model = IDModel(infection, obs)
    lead_in = observation_lead_in(model)
    @test lead_in == 43

    # `n` is the number of observations. The infection series carries the
    # lead-in on top of it, and every observation is scored.
    n = 60
    sim = as_turing_model(model, missing, n)()
    @test length(sim.generated_y_t) == n
    @test length(sim.expected_y_t) == n
    @test length(sim.I_t) == n + lead_in

    # The regression the automatic lead-in exists to remove: perturbing the
    # FIRST observation moves the log-joint, where it used to leave it
    # bit-identical.
    y = fill(100, n)
    vi = VarInfo(as_turing_model(model, y, n))
    lj(yv) = logjoint(as_turing_model(model, yv, n), vi)
    base = lj(y)
    for i in (1, lead_in ÷ 2, lead_in, lead_in + 1, n)
        perturbed = copy(y)
        perturbed[i] = 10^6
        @test lj(perturbed) != base
    end
    @test data_fits(model, y, n)

end

@testitem "a stratified model grows along its time axis only" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1069)

    # Two infection strata feeding one observation stream through a weight map,
    # with a delay on the stream. The lead-in lengthens the time axis and
    # leaves the strata axis alone.
    model = IDModel(
        DirectInfections(;
            Z = Stratify(RandomWalk(), Hierarchy(; across = IID(Normal(0, 0.5)))),
            initialisation = Normal(log(50), 0.2)
        ),
        Split(LatentDelay(PoissonError(), fill(1 / 5, 5)), [1.0 1.0])
    )
    @test observation_lead_in(model) == 4
    Y = Matrix{Union{Missing, Float64}}(missing, 1, 20)
    @test size(as_turing_model(model, Y)().I_t) == (2, 24)
    @test size(as_turing_model(model, Y, (2, 20))().I_t) == (2, 24)
    @test data_requirements(model, Y, (2, 20)).series_length == 24
end

@testitem "a stratified stream counts a flat series on its time axis" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: Model
    Random.seed!(1070)

    # A bare error model scores its expected series over `eachindex(Y_t)`, so a
    # stratified model hands it one flat `n_strata * n_time` series rather than
    # a matrix. Counting that length against the time axis reads `n_strata`
    # times too many observations and rejects a correct call.
    model = IDModel(
        Renewal(;
            generation_time = [0.2, 0.3, 0.5],
            rt = Replicate(RandomWalk()),
            initialisation = Normal(log(50.0), 0.2)
        ),
        PoissonError()
    )
    y = as_turing_model(model, missing, (3, 8))(MersenneTwister(119)).generated_y_t
    @test y isa AbstractVector
    @test length(y) == 24

    # The call the simulation itself implies must be accepted.
    @test as_turing_model(model, y, (3, 8)) isa Model
    @test data_fits(model, y, (3, 8))

    # A genuinely over-long panel is still rejected, on the time axis.
    @test_throws ArgumentError as_turing_model(model, repeat(y, 2), (3, 8))
    # So is a length that is not a whole number of strata.
    @test_throws ArgumentError as_turing_model(model, y[1:23], (3, 8))
end

@testitem "as_turing_model rejects the legacy series length" begin
    using ComposableTuringIDModels, Distributions

    obs = LatentDelay(PoissonError(), fill(1 / 15, 15))
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()), obs
    )
    y = fill(10, 40)

    # The old idiom passed the infection-series length. Under the new meaning
    # that asks for more observations than were supplied, short by exactly the
    # lead-in, so it is named rather than quietly fitting a longer series.
    err = try
        as_turing_model(model, y, length(y) + observation_lead_in(model))
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("number of observations", err.msg)
    @test occursin("old meaning of `n`", err.msg)

    # More observations than the model can score is an error too: their head
    # would never enter the likelihood, which is the failure being removed.
    over = try
        as_turing_model(model, fill(10, 60), 40)
        nothing
    catch e
        e
    end
    @test over isa ArgumentError
    @test occursin("never be scored", over.msg)

    # Fewer is not an error: the data is right-aligned, so a series that starts
    # later is scored at the end and the earlier expected values are unobserved.
    @test as_turing_model(model, fill(10, 30), 40) isa Any

    # Simulating has no data to disagree with.
    @test as_turing_model(model, missing, 40) isa Any
end

@testitem "data_requirements reports what a model needs" begin
    using ComposableTuringIDModels, Distributions

    obs = LatentDelay(LatentDelay(PoissonError(), fill(1 / 15, 15)), fill(1 / 30, 30))
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()), obs
    )

    req = data_requirements(model, 60)
    @test req.n == 60
    @test req.series_length == 60 + 43
    @test length(req) == 1
    stream = req[:y_t]
    @test stream === req[1] === only(collect(req))
    @test stream.name === :y_t
    @test stream.n_required == 60
    @test stream.n_max == 60
    @test stream.n_scored == 60
    @test stream.lead_in == 43
    @test stream.shape === :series
    @test stream.alignment === :right
    @test stream.fields == ()
    @test stream.n_supplied === nothing
    @test data_fits(req)
    @test_throws KeyError req[:nope]

    # With data, each stream also says how much was supplied.
    @test data_requirements(model, fill(10, 60), 60)[:y_t].n_supplied == 60
    # More than the chain can score does not fit; less does, right-aligned.
    @test !data_fits(model, fill(10, 61), 60)
    @test data_fits(model, fill(10, 59), 60)
    @test data_fits(model, missing, 60)

    # An `IDProblem` reads the observation count from its `tspan`.
    problem = IDProblem(
        infection = model.infection_model, observation_model = obs,
        tspan = (1, 60)
    )
    @test data_requirements(problem).series_length == 103
    @test data_fits(data_requirements(problem, (; y_t = fill(10, 60))))
    @test !data_fits(data_requirements(problem, (; y_t = fill(10, 61))))
end

@testitem "data_requirements reads each data shape's contract" begin
    using ComposableTuringIDModels, Distributions

    pmf = fill(1 / 5, 5)                       # lead-in 4
    n = 30
    plain = LatentDelay(PoissonError(), pmf)

    # A plain vector: the series itself carries the observations.
    @test data_requirements(plain, fill(10, n), n)[:y_t].n_supplied == 30

    # A `MissingObservations` carrier reports the series it stands in for.
    carrier = ComposableTuringIDModels.MissingObservations(
        fill(10.0, n), fill(true, n)
    )
    @test data_requirements(plain, carrier, n)[:y_t].n_supplied == 30

    # A stratified model's data is `strata x time`: the time axis counts.
    @test data_requirements(plain, fill(10, 2, n), (2, n))[:y_t].n_supplied == 30

    # `BinomialError` takes `(; y, N)`. The trials are known covariates, not
    # observations, so the contract names both fields and counts only `y`.
    binom = LatentDelay(BinomialError(), pmf)
    trials = data_requirements(binom, (y = fill(3, n), N = 20), n)[:y_t]
    @test trials.fields == (:y, :N)
    @test trials.n_required == 30
    @test trials.n_supplied == 30
    @test data_requirements(
        binom, (y = fill(3, n), N = fill(20, n)), n
    )[:y_t].n_supplied == 30
    # Simulating supplies nothing to disagree with.
    @test data_fits(binom, (y = missing, N = 20), n)

    # A `ReportTriangle` takes a reference-day × delay matrix, counted down the
    # reference days: the delay columns are not observations of their own.
    tri = ReportTriangle(PoissonError(), fill(1 / 3, 3))
    triangle = data_requirements(tri, fill(2, n, 3), n)[:y_t]
    @test triangle.shape === :triangle
    @test triangle.n_required == 30
    @test triangle.n_supplied == 30
    built = define_y_t(tri, fill(2, n, 3), fill(20.0, n))
    @test data_requirements(tri, built, n)[:y_t].n_supplied == 30
    @test data_fits(tri, missing, n)
    # A reporting triangle does not right-align, so its reference days must
    # match exactly: one per observation, under a delay as without one.
    @test triangle.alignment === :exact
    @test !data_fits(tri, fill(2, n - 1, 3), n)
    delayed = LatentDelay(tri, pmf)
    @test data_fits(delayed, fill(2, n, 3), n)
    @test !data_fits(delayed, fill(2, n - 4, 3), n)
end

@testitem "data_requirements reports a Split stream by stream" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(1067)

    n = 40
    streams = Split(
        (
            cases = LatentDelay(PoissonError(), fill(1 / 5, 5)),
            deaths = LatentDelay(PoissonError(), fill(1 / 20, 20)),
        )
    )
    req = data_requirements(streams, n)
    @test length(req) == 2
    @test [s.name for s in req] == [:cases, :deaths]
    # Every stream is supplied `n` observations; what differs between them is
    # the lead-in each consumes, and the series covers the deepest.
    @test all(s -> s.n_required == n, req)
    @test req[:cases].lead_in == 4
    @test req[:deaths].lead_in == 19
    @test req.series_length == n + 19

    # A stream is checked against its own data, and a shared series reaches
    # every stream.
    supplied = data_requirements(
        streams, (cases = fill(10, n), deaths = fill(1, n - 5)), n
    )
    @test supplied[:cases].n_supplied == n
    @test supplied[:deaths].n_supplied == n - 5
    @test data_fits(supplied)
    @test data_fits(streams, fill(10, n), n)

    # One series serves both streams, so it covers the deeper lead-in and the
    # shallower stream can score earlier data if it has any. More than that
    # does not fit.
    @test req[:cases].n_max == n + 19 - 4
    @test req[:deaths].n_max == n
    @test data_fits(streams, (cases = fill(10, n + 15), deaths = fill(1, n)), n)
    @test !data_fits(streams, (cases = fill(10, n + 16), deaths = fill(1, n)), n)

    # A stream taking `(; y, N)` data is read by its own contract, not unpacked
    # as if `N` were another stream.
    mixed = Split((cases = PoissonError(), positivity = BinomialError()))
    cover = data_requirements(
        mixed, (cases = fill(10, n), positivity = (y = fill(3, n), N = 50)), n
    )
    @test length(cover) == 2
    @test cover[:positivity].fields == (:y, :N)
    @test cover[:positivity].n_supplied == n
    @test data_fits(cover)
end

@testitem "data_requirements counts an Aggregate's scored windows" begin
    using ComposableTuringIDModels, Distributions

    n = 28
    weekly = Aggregate(PoissonError(), [0, 0, 0, 0, 0, 0, 7])

    # The caller supplies a full daily series; only the day closing each week
    # reaches the error model, so four of the 28 are scored.
    req = data_requirements(weekly, fill(10, n), n)[:y_t]
    @test req.n_required == 28
    @test req.n_scored == 4
    @test req.lead_in == 0
    @test req.n_supplied == 28
    @test req.alignment === :exact
    @test data_fits(weekly, fill(10, n), n)
    # An aggregation indexes its data against the expected series rather than
    # right-aligning to it, so a shorter series is a call that fails.
    @test !data_fits(weekly, fill(10, n - 1), n)

    # A delay INSIDE the aggregation costs scored windows rather than series
    # length: it consumes the leading windows, which go unpredicted.
    nested = Aggregate(LatentDelay(PoissonError(), fill(1 / 3, 3)), [0, 0, 0, 0, 0, 0, 7])
    inner = data_requirements(nested, fill(10, n), n)
    @test inner.series_length == n
    @test inner[:y_t].n_required == n
    @test inner[:y_t].n_scored == 4 - 2
    @test data_fits(inner)

    # A delay in front of the aggregation lengthens the series, not the data.
    delayed = LatentDelay(weekly, fill(1 / 8, 8))
    @test observation_lead_in(delayed) == 7
    with_delay = data_requirements(delayed, fill(10, n), n)
    @test with_delay.series_length == n + 7
    @test with_delay[:y_t].n_scored == 4
    @test data_fits(with_delay)
end

@testitem "a requirements report prints what to supply" begin
    using ComposableTuringIDModels, Distributions

    streams = Split(
        (
            cases = LatentDelay(PoissonError(), fill(1 / 5, 5)),
            deaths = BinomialError(),
        )
    )
    out = sprint(show, MIME"text/plain"(), data_requirements(streams, 40))
    @test occursin("40 observations per stream", out)
    @test occursin("infection series built at length 44", out)
    @test occursin(
        "cases: 40 values, 4 of them scored, after a lead-in of 4", out
    ) || occursin("cases: 40 values, after a lead-in of 4", out)
    @test occursin("deaths: 40 values", out)
    @test occursin("up to 44 if earlier data exists", out)

    # Supplying data adds what was supplied, and an aggregation says how much
    # of it is scored.
    weekly = Aggregate(PoissonError(), [0, 0, 0, 0, 0, 0, 7])
    agg = sprint(
        show, MIME"text/plain"(), data_requirements(weekly, fill(10, 28), 28)
    )
    @test occursin("4 of them scored", agg)
    @test occursin("supplied 28", agg)
end

@testitem "forecast needs nothing said about the lead-in" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(1068)

    # A delay chain runs its infection series over `length(y) + lead_in` steps,
    # and the forecast model does the same over the horizon, so a caller says
    # only how far ahead to go.
    obs = LatentDelay(PoissonError(), fill(1 / 6, 6))
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal(1.0, 0.5)),
        obs
    )
    y = fill(5, 20)
    @test data_fits(model, y, length(y))

    chain = sample(
        as_turing_model(model, y, length(y)), Prior(), 20; progress = false
    )
    h = 3
    fc = forecast(model, y, chain, h)
    @test size(fc, 1) == 20

    # An `IDProblem`'s `tspan` is the span of the observations.
    problem = IDProblem(
        infection = model.infection_model, observation_model = obs,
        tspan = (1, length(y))
    )
    @test data_fits(data_requirements(problem, (; y_t = y)))
    pchain = sample(
        as_turing_model(problem, (; y_t = y)), Prior(), 20; progress = false
    )
    pfc = forecast(problem, y, pchain, h)
    @test size(pfc, 1) == 20
end

@testitem "a Split of unequal lead-ins scores both streams in full" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: VarInfo, logjoint
    Random.seed!(1070)

    # Each stream is supplied `n` observations and the series covers the
    # deepest lead-in, so the shallower stream produces more expected values
    # than it has data. Right-aligning in both directions leaves that head
    # unobserved and scores every observation of both streams.
    n = 30
    streams = Split(
        (
            cases = LatentDelay(PoissonError(), fill(1 / 5, 5)),
            deaths = LatentDelay(PoissonError(), fill(1 / 20, 20)),
        )
    )
    model = IDModel(
        Renewal(;
            generation_time = [0.2, 0.3, 0.3, 0.2],
            rt = RandomWalk(; init = Normal(0.0, 0.05), ϵ_t = Normal(0.0, 0.01)),
            initialisation = Normal(log(100.0), 0.1)
        ),
        streams
    )
    y = (cases = fill(100, n), deaths = fill(100, n))
    @test data_fits(model, y, n)
    sim = as_turing_model(model, y, n)()
    @test length(sim.generated_y_t.cases) == n
    @test length(sim.generated_y_t.deaths) == n

    # The first observation of each stream is scored: neither is dropped from
    # the head, whichever lead-in its own chain consumes.
    vi = VarInfo(as_turing_model(model, y, n))
    lj(yv) = logjoint(as_turing_model(model, yv, n), vi)
    base = lj(y)
    for k in (:cases, :deaths), i in (1, n)
        perturbed = merge(y, NamedTuple{(k,)}((replace(y[k], 100 => 100),)))
        moved = copy(perturbed[k])
        moved[i] = 10^6
        @test lj(merge(y, NamedTuple{(k,)}((moved,)))) != base
    end
end

@testitem "a requirements report prints compactly and flags a misfit" begin
    using ComposableTuringIDModels, Distributions

    weekly = Aggregate(PoissonError(), [0, 0, 0, 0, 0, 0, 7])

    # An aggregation is the end of the calendar walk: it consumes reporting
    # windows rather than time points, so it adds no lead-in of its own.
    @test observation_lead_in(weekly) == 0

    r = data_requirements(weekly, fill(10, 28), 28)

    # The compact `show` is what appears inside an array or an error message;
    # the rich one is the report a user reads. Both are exercised, because a
    # struct is only as useful as the two ways it prints.
    @test repr(r) == "DataRequirements(28 observations, 1 stream, series 28)"
    @test repr(r[:y_t]) == "y_t: 28 values, 4 of them scored — supplied 28"
    @test sprint(show, MIME"text/plain"(), r[:y_t]) == repr(r[:y_t])

    # `data_fits` answers for one stream as well as for the whole report.
    @test data_fits(r[:y_t])
    @test data_fits(r)

    # A stream that does not fit says so in capitals rather than quietly, and
    # `data_fits` agrees with the print.
    short = data_requirements(weekly, fill(10, 20), 28)
    out = sprint(show, MIME"text/plain"(), short)
    @test occursin("SUPPLIED 20", out)
    @test !occursin("— supplied 20", out)
    @test !data_fits(short)
    @test !data_fits(short[:y_t])

    # An aggregation is indexed against the expected series rather than
    # right-aligned to it, so a length that differs is rejected outright
    # instead of being scored at the end.
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()), weekly
    )
    msg = try
        as_turing_model(model, fill(10, 20), 28)
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("takes exactly 28 values but 20 were supplied", msg)
    @test occursin("indexed against the expected series", msg)

    # A reporting triangle is the other exact-length component.
    tri = ReportTriangle(PoissonError(), fill(0.25, 4))
    @test ComposableTuringIDModels._alignment(tri) == :exact
    @test ComposableTuringIDModels._alignment(weekly) == :exact
    # Everything else right-aligns, which is the fallback the walk lands on.
    @test ComposableTuringIDModels._alignment(PoissonError()) == :right
    @test ComposableTuringIDModels._needs_input_calendar(weekly)
    @test !ComposableTuringIDModels._needs_input_calendar(PoissonError())
end

@testitem "observation_streams reads a data value's stream axis" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: observation_streams
    # The default reads the value's own axes.
    @test observation_streams(PoissonError(), fill(5.0, 2, 10)) == 2
    @test observation_streams(PoissonError(), fill(5.0, 10)) === nothing
    @test observation_streams(PoissonError(), missing) === nothing
    # An error family's `NamedTuple` is one stream's observations in `y` beside
    # its covariates, so the count comes off `y`'s own axis.
    @test observation_streams(PoissonError(), (y = fill(5.0, 10),)) === nothing
    @test observation_streams(PoissonError(), (y = fill(5.0, 2, 10),)) == 2

    # A modifier passes the question down to whatever consumes the series.
    delayed = LatentDelay(PoissonError(), fill(1 / 4, 4))
    @test observation_streams(delayed, fill(5.0, 3, 10)) == 3

    # A component whose data means something else says so.
    tri = ReportTriangle(PoissonError(), [0.5, 0.3, 0.2])
    @test observation_streams(tri, [10 5 2; 12 6 3]) === nothing
    @test observation_streams(LatentDelay(tri, [0.5, 0.5]), [10 5 2]) === nothing
    binom = BinomialError()
    @test observation_streams(binom, (y = fill(5, 10), N = fill(20, 10))) ===
        nothing
    @test observation_streams(
        binom, (y = fill(5, 2, 10), N = fill(20, 2, 10))
    ) == 2

    # A `Split` slices the value before any stream sees it, so its outer axis is
    # streams whatever a stream then reads.
    @test observation_streams(
        Split(tri), (a = [10 5 2], b = [12 6 3])
    ) == 2

    # An IDModel and an IDProblem delegate to their observation model.
    infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    @test observation_streams(IDModel(infection, tri), [10 5 2; 12 6 3]) ===
        nothing
    problem = IDProblem(
        infection = infection, observation_model = delayed, tspan = (1, 10)
    )
    @test observation_streams(problem, fill(5.0, 3, 10)) == 3
end
