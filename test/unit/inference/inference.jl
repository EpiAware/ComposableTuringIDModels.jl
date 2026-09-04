@testitem "IDProblem assembles and simulates a composed model" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(71)
    problem = IDProblem(
        infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        observation_model = PoissonError(),
        tspan = (1, 20)
    )
    m = as_turing_model(problem, (; y_t = missing))
    sim = m()
    @test length(sim.generated_y_t) == 20
    @test length(sim.Z_t) == 20
end

@testitem "_obs_data_shape resolves stratified shapes from data or Split" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: _obs_data_shape
    # `IDProblem` and `forecast` share this helper, so every branch it can
    # take is exercised directly here rather than only through the callers.
    plain = PoissonError()
    Y = fill(5.0, 2, 10)
    @test _obs_data_shape(plain, Y, 10) == (2, 10)   # AbstractMatrix y_t
    mapped = Split(PoissonError(), [1.0 1.0 1.0])    # 1 stream, 3 inf strata
    @test _obs_data_shape(mapped, fill(5.0, 1, 10), 10) == (3, 10)

    named = Split((a = PoissonError(), b = PoissonError()))
    y_nt = (a = fill(5.0, 10), b = fill(3.0, 10))
    @test _obs_data_shape(named, y_nt, 10) == (2, 10)   # NamedTuple y_t

    # A `Split` that fixes its own stratum count resolves the same shape with
    # or without data.
    @test _obs_data_shape(plain, missing, 10) == 10
    @test _obs_data_shape(mapped, missing, 10) == (3, 10)
    @test _obs_data_shape(named, missing, 10) == (2, 10)

    # A strata template takes its stream names from the data, so with
    # `y_t = missing` there is nothing to read a stratum count from and the
    # shape falls back to a single series.
    strata_template = Split(PoissonError())
    @test _obs_data_shape(strata_template, missing, 10) == 10
    @test _obs_data_shape(strata_template, y_nt, 10) == (2, 10)
end

@testitem "a data axis stratifies only where the component says it is streams" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: _obs_data_shape
    # A reporting triangle's rows are reference days and its columns reporting
    # delays, so neither axis is a stream axis.
    # The dense matrix and the built carrier are two spellings of one dataset
    # and must give the same shape.
    tri = ReportTriangle(PoissonError(), [0.5, 0.3, 0.2])
    N = [10 5 2; 12 6 3; 14 7 4]
    carrier = define_y_t(tri, N, fill(20.0, 3))
    @test _obs_data_shape(tri, N, 3) == 3
    @test _obs_data_shape(tri, carrier, 3) == 3
    @test _obs_data_shape(tri, missing, 3) == 3
    # A modifier in front of the triangle does not change what its data means.
    @test _obs_data_shape(LatentDelay(tri, [0.5, 0.5]), N, 3) == 3

    # `BinomialError` takes its trials beside its counts, so the two fields are
    # one stream rather than two, and a panel keeps its stream axis.
    binom = BinomialError()
    @test _obs_data_shape(binom, (y = fill(5, 10), N = fill(20, 10)), 10) == 10
    @test _obs_data_shape(
        binom, (y = fill(5, 2, 10), N = fill(20, 2, 10)), 10
    ) == (2, 10)

    # Every error family takes its observations in `y`, so wrapping a series in
    # a `NamedTuple` must not change the model it builds.
    plain = PoissonError()
    @test _obs_data_shape(plain, (y = fill(5, 10),), 10) ==
        _obs_data_shape(plain, fill(5, 10), 10)
    @test _obs_data_shape(plain, (y = fill(5, 2, 10),), 10) ==
        _obs_data_shape(plain, fill(5, 2, 10), 10)
end

@testitem "a nested Split's weight map fixes the stratum count" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: _obs_data_shape
    # The map states the model's stratum count, so a modifier wrapped round the
    # `Split` does not hide it and the count holds with no data to read.
    nested = LatentDelay(Split(PoissonError(), [1.0 1.0 1.0]), [0.5, 0.3, 0.2])
    @test _obs_data_shape(nested, fill(5.0, 1, 10), 10) == (3, 10)
    @test _obs_data_shape(nested, missing, 10) == (3, 10)
end

@testitem "the resolved shape is inferred rather than a union" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: _obs_data_shape
    # The shape is built in a `@model` body, so a union return type would make
    # every downstream variable in the model unpredictable to the compiler.
    tri = ReportTriangle(PoissonError(), [0.5, 0.3, 0.2])
    cases = [
        (PoissonError(), fill(5.0, 12)),
        (PoissonError(), fill(5.0, 2, 12)),
        (BinomialError(), (y = fill(5, 12), N = fill(20, 12))),
        (tri, [10 5 2; 12 6 3; 14 7 4]),
        (tri, missing),
        (
            LatentDelay(Split(PoissonError(), [1.0 1.0 1.0]), [0.5, 0.5]),
            fill(5.0, 1, 12),
        ),
        (
            Split((a = PoissonError(), b = PoissonError())),
            (a = fill(5.0, 12), b = fill(5.0, 12)),
        ),
    ]
    for (obs, y_t) in cases
        types = Base.return_types(
            _obs_data_shape, Tuple{typeof(obs), typeof(y_t), Int}
        )
        @test length(types) == 1
        @test isconcretetype(only(types))
    end
end

@testitem "IDProblem builds a single-series model for a BinomialError" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: VarInfo
    Random.seed!(246)
    # `(y, N)` is one stream's counts and its trials, so the infection process
    # is a single unstratified series.
    problem = IDProblem(
        infection = DirectInfections(;
            Z = RandomWalk(), initialisation = Normal(log(50), 0.2)
        ),
        observation_model = BinomialError(),
        tspan = (1, 10)
    )
    m = as_turing_model(problem, (; y_t = (y = fill(5, 10), N = fill(20, 10))))
    @test Set(string.(keys(VarInfo(m)))) ==
        Set(["init", "std", "ϵ_t", "init_incidence"])
    @test length(m().Z_t) == 10
end

@testitem "IDProblem builds one series for a reporting triangle either way" begin
    using ComposableTuringIDModels, Distributions, Random
    using DynamicPPL: VarInfo
    Random.seed!(2460)
    # A dense count matrix and the carrier built from it are the same dataset,
    # so they must build the same model.
    obs = ReportTriangle(PoissonError(), [0.5, 0.3, 0.2])
    N = [10 5 2; 12 6 3; 14 7 4]
    problem = IDProblem(
        infection = DirectInfections(;
            Z = RandomWalk(), initialisation = Normal(log(50), 0.2)
        ),
        observation_model = obs,
        tspan = (1, 3)
    )
    from_matrix = as_turing_model(problem, (; y_t = N))
    from_carrier = as_turing_model(
        problem, (; y_t = define_y_t(obs, N, fill(20.0, 3)))
    )
    @test keys(VarInfo(from_matrix)) == keys(VarInfo(from_carrier))
    @test Set(string.(keys(VarInfo(from_matrix)))) ==
        Set(["init", "std", "ϵ_t", "init_incidence"])
end

@testitem "IDProblem resolves a stratified infection process via Split" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(710)
    T = 12
    # A Split with a weight map fixes the infection-stratum count at the
    # map's column count, regardless of the (1-stream) data.
    problem_map = IDProblem(
        infection = DirectInfections(;
            Z = Stratify(
                RandomWalk(), Hierarchy(; across = IID(Normal(0, 0.5)))
            ),
            initialisation = Normal(log(20), 0.2)
        ),
        observation_model = Split(PoissonError(), [1.0 1.0 1.0]),
        tspan = (1, T)
    )
    Y = Matrix{Union{Missing, Float64}}(missing, 1, T)
    sim_map = as_turing_model(problem_map, (; y_t = Y))()
    @test size(sim_map.I_t) == (3, T)

    # A Split with named streams and `y_t = missing` builds one stratum per
    # name.
    problem_named = IDProblem(
        infection = DirectInfections(;
            Z = Stratify(
                RandomWalk(), Hierarchy(; across = IID(Normal(0, 0.5)))
            ),
            initialisation = Normal(log(20), 0.2)
        ),
        observation_model = Split((a = PoissonError(), b = PoissonError())),
        tspan = (1, T)
    )
    sim_named = as_turing_model(problem_named, (; y_t = missing))()
    @test size(sim_named.I_t) == (2, T)
end

@testitem "apply_method runs a NUTSampler over an IDProblem" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(72)
    problem = IDProblem(
        infection = DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        observation_model = PoissonError(),
        tspan = (1, 20)
    )
    ydata = as_turing_model(problem, (; y_t = missing))().generated_y_t
    res = apply_method(problem, NUTSampler(; ndraws = 40, nchains = 1), (; y_t = ydata))
    @test res isa IDObservables
    @test res.samples !== nothing
    # Generated quantities are now recovered from the chain (previously always
    # `missing`); `returned` yields the model's `(; generated_y_t, I_t, Z_t)`
    # per draw.
    @test res.generated !== missing
end

@testitem "spread_draws produces tidy draw/chain/iteration columns" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, MCMCChains, Random
    Random.seed!(74)
    @model f() = (x ~ Normal())
    chn = MCMCChains.Chains(sample(f(), NUTS(), 30; progress = false))
    df = spread_draws(chn)
    @test all(c -> c in names(df), ["draw", "chain", "iteration"])
    @test size(df, 1) == 30
end

@testitem "generated_observables wraps model, data, and solution" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(75)
    m = as_turing_model(
        IDModel(
            DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
            PoissonError()
        ), missing, 10
    )
    obs = generated_observables(m, (; y_t = missing), rand(m))
    @test obs isa IDObservables
    @test obs.model === m
    # `returned` recovers the model's generated quantities for any draw it can
    # consume (a chain or, as here, a single draw); only solutions it cannot
    # consume (e.g. an optimiser result) leave `generated` as `missing`.
    @test obs.generated !== missing
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
    problem = IDProblem(
        infection = DirectInfections(; Z = AR(), initialisation = Normal()),
        observation_model = PoissonError(),
        tspan = (1, 18)
    )
    T, h = 18, 5
    y = as_turing_model(problem, (; y_t = missing))().generated_y_t
    chain = sample(
        as_turing_model(problem, (; y_t = y)), Prior(), 30;
        progress = false
    )
    fc = forecast(problem, y, chain, h)
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

@testitem "generated_observables leaves non-chain solutions missing" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(76)
    m = as_turing_model(
        IDModel(
            DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
            PoissonError()
        ), missing, 10
    )
    # A solution `returned` cannot consume (here a bare marker) has no generated
    # quantities, so the field stays `missing`.
    obs = generated_observables(m, (; y_t = missing), :no_solution)
    @test obs.generated === missing
    # And the untyped fallback: a non-model `model` also yields `missing`.
    @test generated_observables(:not_a_model, (; y_t = missing), :no_solution).generated ===
        missing
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
