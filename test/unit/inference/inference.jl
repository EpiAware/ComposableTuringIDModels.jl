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
    using ComposableTuringIDModels: _obs_data_shape, _obs_data_shape_missing
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

    # `y_t === missing` falls back to `_obs_data_shape_missing`; its three
    # `Split` branches (map / names / neither) are each checked in turn.
    @test _obs_data_shape(plain, missing, 10) == 10
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
    # Regression for #140. An un-damped `DiffLatentModel(RandomWalk())` on
    # `log R_t` is prior-explosive: incidence can reach `Inf` on some prior
    # draws, and `forecast`'s `_assert_factorised` probe samples the whole
    # model (including the observation) 256 times to check the correctness
    # guard. That used to surface as an opaque `InexactError: BigInt(Inf)`
    # from deep inside `SafePoisson`'s sampler (`_safe_int_floor`); it must
    # now fail with a clear `DomainError` naming the non-finite rate instead.
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
    threw_inexact = try
        forecast(model, y, chain, 6)
        false
    catch e
        e isa InexactError
    end
    @test !threw_inexact
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
