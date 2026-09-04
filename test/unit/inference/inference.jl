@testitem "IDProblem pairs a model with its data and needs no length" begin
    using ComposableTuringIDModels, Distributions, Random, Accessors
    using DynamicPPL: VarInfo
    Random.seed!(71)
    infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    problem = IDProblem(infection, PoissonError(), Vector{Missing}(missing, 20))
    @test problem.model == IDModel(infection, PoissonError())
    sim = as_turing_model(problem)()
    @test length(sim.generated_y_t) == 20
    @test length(sim.Z_t) == 20

    # The problem's route and the model's route build the same model, so
    # holding one is a convenience rather than a fork in the API.
    model = IDModel(infection, PoissonError())
    y = sim.generated_y_t
    @test keys(VarInfo(as_turing_model(IDProblem(model, y)))) ==
        keys(VarInfo(as_turing_model(model, y, length(y))))

    # Refitting to a different series is a new problem, built with `Accessors`,
    # rather than a second argument to `as_turing_model`.
    refit = @set problem.data = fill(5, 12)
    @test isequal(refit.model, problem.model)
    @test isequal(refit.data, fill(5, 12))
    @test isequal(problem.data, Vector{Missing}(missing, 20))
    @test length(as_turing_model(refit)().generated_y_t) == 12
end

@testitem "IDProblem refuses streams that disagree at construction" begin
    using ComposableTuringIDModels, Distributions, Accessors
    infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    # Two streams over one infection stratum, so the observation count is read
    # off `cases` and `deaths` can genuinely contradict it.
    obs = Split(
        (cases = PoissonError(), deaths = PoissonError()),
        reshape([1.0, 1.0], 2, 1)
    )
    ragged = (cases = fill(5, 30), deaths = fill(1, 40))
    @test_throws ArgumentError IDProblem(infection, obs, ragged)
    # The guard is an inner constructor, so respecifying into the same mistake
    # is refused too.
    ok = IDProblem(infection, obs, (cases = fill(5, 30), deaths = fill(1, 30)))
    @test_throws ArgumentError @set ok.data = ragged

    # It never refuses what `as_turing_model` accepts. A shorter stream is
    # right-aligned rather than wrong, and a single series cannot disagree with
    # a count read off itself.
    @test IDProblem(infection, obs, (cases = fill(5, 30), deaths = fill(1, 20))) isa
        IDProblem
    @test IDProblem(infection, PoissonError(), fill(5, 30)) isa IDProblem
    # Both simulation routes stay constructible.
    @test IDProblem(infection, PoissonError(), Vector{Missing}(missing, 20)) isa
        IDProblem
    @test IDProblem(
        infection, obs, (
            cases = Vector{Missing}(missing, 20), deaths = Vector{Missing}(missing, 20),
        )
    ) isa IDProblem
    # And a model reading named fields rather than streams is not mistaken for
    # two disagreeing streams.
    @test IDProblem(
        infection, BinomialError(), (y = fill(5, 30), N = fill(20, 30))
    ) isa IDProblem
end

@testitem "a blank series is a problem awaiting its observations" begin
    using ComposableTuringIDModels, Distributions, Accessors, Random
    using DynamicPPL: VarInfo
    Random.seed!(73)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal(log(20), 0.2)),
        LatentDelay(PoissonError(), fill(1 / 5, 5))
    )
    # There is no data-free problem; a blank series is how a shape is fixed
    # before the observations arrive. It carries a length, so it simulates.
    waiting = IDProblem(model, Vector{Missing}(missing, 20))
    @test length(as_turing_model(waiting)().generated_y_t) == 20
    @test data_requirements(waiting).n == 20

    # A scalar `missing` has no length, so it is not a way to express this.
    @test_throws ArgumentError IDProblem(model, missing)

    # Attaching the data gives what the constructor gives, field for field.
    y = fill(5, 20)
    attached = @set waiting.data = y
    built = IDProblem(model, y)
    @test all(
        f -> isequal(getfield(attached, f), getfield(built, f)),
        fieldnames(IDProblem)
    )
    @test keys(VarInfo(as_turing_model(attached))) ==
        keys(VarInfo(as_turing_model(built)))
    # The blank problem is unchanged, since `@set` builds rather than mutates.
    @test isequal(waiting.data, Vector{Missing}(missing, 20))
end

@testitem "IDProblem refuses a bare length in place of its data" begin
    using ComposableTuringIDModels, Distributions
    problem = IDProblem(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError(), fill(5, 30)
    )
    # The generic method would report as if no data had been supplied, which
    # discards the half of the pairing the type exists for.
    @test_throws ArgumentError data_requirements(problem, 30)
    @test_throws ArgumentError data_requirements(problem, (2, 30))
end

@testitem "IDProblem prints its component tree and a data summary" begin
    using ComposableTuringIDModels, Distributions
    infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    render(x) = sprint(show, MIME"text/plain"(), x)

    plain = render(IDProblem(infection, PoissonError(), fill(5, 30)))
    @test startswith(plain, "IDProblem\n")
    # The model's components hang off the problem rather than under a node of
    # their own, and the data is the last sibling.
    @test occursin("├─ infection: DirectInfections", plain)
    @test occursin("├─ observation: PoissonError", plain)
    @test occursin("└─ data: 30 observations (Int64)", plain)

    # A wholly blank series says it is simulating rather than reporting zero
    # observations.
    blank = render(IDProblem(infection, PoissonError(), Vector{Missing}(missing, 30)))
    @test occursin("└─ data: none, 30 time points (simulating from the prior)", blank)

    gappy = render(
        IDProblem(
            infection, PoissonError(),
            Vector{Union{Missing, Int}}([missing; fill(5, 29)])
        )
    )
    @test occursin("└─ data: 30 observations, 1 missing", gappy)

    # A `NamedTuple` names its entries without claiming they are streams: they
    # are for a `Split`, but a `BinomialError` reads `(y, N)` as one stream's
    # fields.
    streams = render(
        IDProblem(
            infection, Split((cases = PoissonError(), deaths = PoissonError())),
            (cases = fill(5, 30), deaths = fill(1, 30))
        )
    )
    @test occursin("└─ data: 30 observations in each of cases, deaths", streams)

    binom = render(
        IDProblem(
            infection, BinomialError(),
            (y = fill(5, 30), N = fill(20, 30))
        )
    )
    @test occursin("└─ data: 30 observations in each of y, N", binom)

    blank_streams = render(
        IDProblem(
            infection, Split((cases = PoissonError(), deaths = PoissonError())),
            (
                cases = Vector{Missing}(missing, 30),
                deaths = Vector{Missing}(missing, 30),
            )
        )
    )
    @test occursin(
        "└─ data: none, 30 time points (simulating from the prior)", blank_streams
    )

    strata = render(
        IDProblem(
            DirectInfections(;
                Z = Stratify(RandomWalk(), FixedIntercept(0.0)),
                initialisation = IID(Normal())
            ),
            PoissonError(), fill(5, 2, 30)
        )
    )
    @test occursin("└─ data: 2 strata x 30 observations", strata)

    blank_strata = render(
        IDProblem(
            DirectInfections(;
                Z = Stratify(RandomWalk(), FixedIntercept(0.0)),
                initialisation = IID(Normal())
            ),
            PoissonError(), Matrix{Missing}(missing, 2, 30)
        )
    )
    @test occursin(
        "└─ data: none, 30 time points (simulating from the prior)", blank_strata
    )

    # Compact rendering stays one line, as it does for a component.
    @test sprint(show, IDProblem(infection, PoissonError(), fill(5, 30))) ==
        "IDProblem"
end

@testitem "IDProblem reports its own data requirements" begin
    using ComposableTuringIDModels, Distributions
    infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal())
    obs = LatentDelay(PoissonError(), fill(1 / 5, 5))       # lead-in 4

    problem = IDProblem(infection, obs, fill(10, 30))
    @test data_requirements(problem).n == 30
    @test data_requirements(problem).series_length == 34
    # The report is the model's over the same data, printed the same way.
    report(x...) = sprint(show, MIME"text/plain"(), data_requirements(x...))
    @test report(problem) == report(problem.model, fill(10, 30), 30)

    # A stratified problem reports on the stratum axis, because the shape comes
    # from the same helper `as_turing_model` builds with.
    strata = IDProblem(
        DirectInfections(;
            Z = Stratify(RandomWalk(), FixedIntercept(0.0)),
            initialisation = IID(Normal())
        ),
        obs, fill(10, 2, 30)
    )
    @test data_requirements(strata).n == 30
end

@testitem "_obs_data_shape resolves stratified shapes from data or Split" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: _obs_data_shape, _obs_data_shape_missing
    # The two-argument `as_turing_model`, `IDProblem` and `forecast` share this
    # helper, so every branch it can take is exercised directly here rather
    # than only through the callers.
    plain = PoissonError()
    Y = fill(5.0, 2, 10)
    @test _obs_data_shape(plain, Y, 10) == (2, 10)   # AbstractMatrix y_t
    mapped = Split(PoissonError(), [1.0 1.0 1.0])    # 1 stream, 3 inf strata
    @test _obs_data_shape(mapped, fill(5.0, 1, 10), 10) == (3, 10)

    named = Split((a = PoissonError(), b = PoissonError()))
    y_nt = (a = fill(5.0, 10), b = fill(3.0, 10))
    @test _obs_data_shape(named, y_nt, 10) == (2, 10)   # NamedTuple y_t

    # `y_t === missing` falls back to `_obs_data_shape_missing`; its three
    # `Split` branches (map / names / neither) are each checked in turn.
    @test _obs_data_shape(plain, missing, 10) == 10
    # A blank series carries a length but no stream axis, so it reaches the
    # same fallback a scalar `missing` does.
    @test _obs_data_shape(mapped, Vector{Missing}(missing, 10), 10) == (3, 10)
    @test _obs_data_shape_missing(mapped, 10) == (3, 10)   # map branch
    @test _obs_data_shape_missing(named, 10) == (2, 10)    # names branch
    # Data-driven strata mode: no map, no names (both set at data time), so
    # with `y_t = missing` there is nothing to read a stratum count from and
    # the shape falls back to a single series.
    strata_template = Split(PoissonError())
    @test _obs_data_shape_missing(strata_template, 10) == 10   # neither branch
end

@testitem "IDProblem resolves a stratified infection process via Split" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(710)
    T = 12
    # A Split with a weight map fixes the infection-stratum count at the
    # map's column count, regardless of the (1-stream) data.
    infection = DirectInfections(;
        Z = Stratify(RandomWalk(), Hierarchy(; across = IID(Normal(0, 0.5)))),
        initialisation = Normal(log(20), 0.2)
    )
    Y = Matrix{Union{Missing, Float64}}(missing, 1, T)
    problem_map = IDProblem(infection, Split(PoissonError(), [1.0 1.0 1.0]), Y)
    @test size(as_turing_model(problem_map)().I_t) == (3, T)

    # A Split with named streams over a blank series builds one stratum per
    # name, reading the length off the series.
    problem_named = IDProblem(
        infection, Split((a = PoissonError(), b = PoissonError())),
        Vector{Missing}(missing, T)
    )
    @test size(as_turing_model(problem_named)().I_t) == (2, T)
end

@testitem "forecast rejects a non-positive horizon" begin
    using ComposableTuringIDModels, Distributions
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError()
    )
    @test_throws ArgumentError forecast(model, fill(5, 10), :chain, 0)
end

@testitem "forecast extends a RandomWalk model over the horizon" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(101)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal(1.0, 0.5)),
        PoissonError()
    )
    T, h = 15, 6
    y = as_turing_model(model, fill(missing, T), T)().generated_y_t
    chain = sample(as_turing_model(model, y, T), Prior(), 40; progress = false)
    fc = forecast(model, y, chain, h)
    # The horizon points are predicted and integer-valued counts.
    @test size(fc, 1) == 40
    fut = vec(fc[@varname(y_t[T + 1])])
    @test length(fut) == 40
    @test all(x -> x isa Integer && x ≥ 0, fut)
    @test vec(fc[@varname(y_t[T + h])]) |> length == 40
    # The extended latent path continues the fitted trajectory rather than
    # overwriting it: in-sample Zₜ is unchanged and the path reaches T + h.
    gens_fit = vec(returned(as_turing_model(model, y, T), chain))
    gens_fc = vec(
        returned(
            as_turing_model(
                model, vcat(y, fill(missing, h)),
                T + h
            ), fc
        )
    )
    @test length(gens_fc[1].Z_t) == T + h
    @test all(
        d -> isapprox(gens_fit[d].Z_t, gens_fc[d].Z_t[1:T]; atol = 1.0e-8),
        1:10
    )
end

@testitem "forecast refuses a correlated (non-iid) latent stream" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random, LinearAlgebra
    Random.seed!(104)
    # A deliberately correlated latent: its stored stream is a smooth MvNormal,
    # so it does NOT factorise across time. Independent-tail extension would be
    # statistically wrong, and forecast must refuse rather than mis-forecast.
    struct CorrLatent <: AbstractPriorModel end
    @model function ComposableTuringIDModels.as_turing_model(::CorrLatent, n)
        Σ = [exp(-abs(i - j) / 5) for i in 1:n, j in 1:n] + 1.0e-6 * I
        z ~ MvNormal(zeros(n), Σ)
        return z
    end
    model = IDModel(
        DirectInfections(; Z = CorrLatent(), initialisation = Normal()),
        PoissonError()
    )
    T, h = 15, 5
    y = as_turing_model(model, fill(missing, T), T)().generated_y_t
    chain = sample(as_turing_model(model, y, T), Prior(), 30; progress = false)
    @test_throws ErrorException forecast(model, y, chain, h)
end

@testitem "forecast works through an IDProblem and an AR latent" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(102)
    T, h = 18, 5
    infection = DirectInfections(; Z = AR(), initialisation = Normal())
    y = as_turing_model(
        IDProblem(infection, PoissonError(), Vector{Missing}(missing, T))
    )().generated_y_t
    # The fitted problem holds the series, so the forecast needs only the chain
    # and the horizon.
    problem = IDProblem(infection, PoissonError(), y)
    chain = sample(as_turing_model(problem), Prior(), 30; progress = false)
    fc = forecast(problem, chain, h)
    @test size(fc, 1) == 30
    @test length(vec(fc[@varname(y_t[T + h])])) == 30
end

@testitem "forecast extends a stratified model given a matrix y" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(105)
    model = IDModel(
        DirectInfections(;
            Z = Stratify(RandomWalk(), FixedIntercept(0.0)),
            initialisation = IID(Normal(log(20.0), 0.3))
        ),
        PoissonError()
    )
    n_strata, T, h = 2, 12, 4
    y = as_turing_model(
        model, fill(missing, n_strata, T), (n_strata, T)
    )().generated_y_t
    chain = sample(
        as_turing_model(model, y, (n_strata, T)), Prior(), 30;
        progress = false
    )
    fc = forecast(model, y, chain, h)
    @test size(fc, 1) == 30
    # Rebuilding at the extended shape and reading the generated quantities
    # back off the forecast draws confirms the matrix-shaped `y` was resolved
    # to the right stratified shape throughout (build, sample, and forecast).
    gens_fc = vec(
        returned(
            as_turing_model(
                model, hcat(y, fill(missing, n_strata, h)),
                (n_strata, T + h)
            ), fc
        )
    )
    @test all(d -> size(gens_fc[d].I_t) == (n_strata, T + h), 1:10)
end

@testitem "forecast extends a stratified model given a NamedTuple y" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(106)
    model = IDModel(
        DirectInfections(;
            Z = Stratify(RandomWalk(), FixedIntercept(0.0)),
            initialisation = IID(Normal(log(20.0), 0.3))
        ),
        Split((a = PoissonError(), b = PoissonError()))
    )
    T, h = 12, 4
    y = as_turing_model(
        model, (a = missing, b = missing), (2, T)
    )().generated_y_t
    chain = sample(
        as_turing_model(model, y, (2, T)), Prior(), 30; progress = false
    )
    fc = forecast(model, y, chain, h)
    @test size(fc, 1) == 30
    y_ext = map(v -> vcat(v, fill(missing, h)), y)
    gens_fc = vec(returned(as_turing_model(model, y_ext, (2, T + h)), fc))
    @test all(d -> size(gens_fc[d].I_t) == (2, T + h), 1:10)
end

@testitem "forecast extends a renewal with an inferred generation interval" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(103)
    # An INFERRED generation interval (its distribution parameters carry priors)
    # must not break forecasting: the interval's only RV is the fixed-length
    # parameter draw `θ`, not a per-time stream, so the non-centred horizon
    # extension leaves it untouched.
    gen = UncertainDelay(
        LogNormal,
        [Normal(1.9, 0.2), truncated(Normal(0.5, 0.2), 0, Inf)]; D = 10.0
    )
    model = IDModel(
        Renewal(;
            generation_time = gen, rt = RandomWalk(),
            initialisation = Normal()
        ),
        PoissonError()
    )
    T, h = 18, 5
    y = as_turing_model(model, fill(missing, T), T)().generated_y_t
    chain = sample(as_turing_model(model, y, T), Prior(), 30; progress = false)
    fc = forecast(model, y, chain, h)
    @test size(fc, 1) == 30
    @test length(vec(fc[@varname(y_t[T + h])])) == 30
    @test all(x -> x isa Integer && x ≥ 0, vec(fc[@varname(y_t[T + 1])]))
end

@testitem "forecast's factorisation probe fails clearly on an explosive prior" tags = [
    :sample,
] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(2024)
    # An un-damped `DiffLatentModel(RandomWalk())` on `log R_t` is
    # prior-explosive: incidence can reach `Inf` on some prior draws, and
    # `forecast`'s `_assert_factorised` probe samples the whole model,
    # observation included. A non-finite rate reaching `SafePoisson` must
    # raise a `DomainError` naming the value, not an `InexactError` from
    # deep inside the sampler.
    diffrt = DiffLatentModel(model = RandomWalk(), init = [Normal(0.0, 0.2)])
    ren = Renewal(
        generation_time = Gamma(6.5, 0.62); rt = diffrt,
        initialisation = Normal(log(1.0), 0.1)
    )
    model = IDModel(ren, PoissonError())
    y = fill(5, 20)
    chain = sample(
        as_turing_model(model, y, length(y)), Prior(), 12; progress = false
    )
    @test_throws DomainError forecast(model, y, chain, 6)
end

@testitem "forecast generates its horizon from the posterior" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(104)
    model = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal(1.0, 0.5)),
        PoissonError()
    )
    T, h = 15, 5
    y = Vector{Union{Missing, Int}}(
        as_turing_model(model, fill(missing, T), T)().generated_y_t
    )
    y[[4, 9]] .= missing
    chain = sample(as_turing_model(model, y, T), Prior(), 30; progress = false)
    # Nothing about `y_t` is a parameter of the fit, gaps included, so the
    # horizon cannot be read out of the posterior.
    @test !any(k -> occursin("y_t", string(k)), collect(keys(chain)))

    fc = forecast(model, y, chain, h)
    # It is generated instead: the horizon comes back as counts, per draw.
    @test all(x -> x isa Integer && x >= 0, vec(fc[@varname(y_t[T + 1])]))
    @test length(vec(fc[@varname(y_t[T + h])])) == 30
    # Generating rather than reading gives the in-sample points too, so the
    # same chain carries the in-sample posterior predictive at the gaps.
    @test length(vec(fc[@varname(y_t[4])])) == 30
    @test all(x -> x isa Integer && x >= 0, vec(fc[@varname(y_t[9])]))
end
