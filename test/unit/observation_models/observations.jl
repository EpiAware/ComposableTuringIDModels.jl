@testitem "observation error models simulate and condition" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(11)
    Y_t = fill(10.0, 15)
    for obs in (PoissonError(), NegativeBinomialError())
        # Simulate from the prior (missing observations). Every observation model
        # returns the uniform `(; y_t, expected)` contract.
        res = as_turing_model(obs, missing, Y_t)()
        @test keys(res) == (:y_t, :expected)
        @test res.expected == Y_t
        sim = res.y_t
        @test length(sim) == length(Y_t)
        @test all(>=(0), sim)
        # Condition on simulated data: the model still builds and evaluates.
        cond = as_turing_model(obs, sim, Y_t)
        @test cond().y_t == sim
    end
end

@testitem "NormalError is a continuous observation-error model" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(111)

    ne = NormalError()
    # Role + interface conformance (the post-role-hierarchy contract).
    @test ne isa AbstractObservationErrorModel
    @test ne isa AbstractObservationModel
    @test implements_observation_interface(ne)

    # Default std prior is positive (a HalfNormal).
    @test minimum(ne.std) >= 0

    Y_t = fill(10.0, 20)
    # Simulate from the prior: continuous (real) observations. Simulating from a
    # `missing` series returns a `Union{Missing,Float64}` vector, so check the
    # sampled values are real rather than the container eltype.
    sim = as_turing_model(ne, missing, Y_t)().y_t
    @test length(sim) == length(Y_t)
    @test all(x -> x isa Real, sim)
    # Condition on the simulated data: the model still builds and evaluates,
    # returning the same observations.
    cond = as_turing_model(ne, sim, Y_t)
    @test cond().y_t == sim
    # The standard deviation is an inferred parameter.
    draw = rand(as_turing_model(ne, sim, Y_t))
    @test any(k -> occursin("σ", string(k)), keys(draw))

    # A custom std prior is honoured.
    ne2 = NormalError(; std = truncated(Normal(0, 2), 0, Inf))
    @test length(as_turing_model(ne2, missing, fill(5.0, 6))().y_t) == 6

    # Expected-mean alignment: conditioning on data with a matching expected
    # series evaluates to that data (the Gaussian likelihood is centred on Y_t).
    μ = 50.0
    obs = μ .+ 0.5 .* randn(200)
    m = as_turing_model(
        NormalError(; std = truncated(Normal(0, 1), 0, Inf)),
        obs, fill(μ, 200)
    )
    @test m().y_t == obs
end

@testitem "BinomialError reads trials from NamedTuple data" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(222)

    be = BinomialError()
    # The struct carries no data.
    @test isempty(fieldnames(BinomialError))
    # Role + interface conformance. The expected series is a probability, and the
    # number of trials is supplied via the NamedTuple data `y_t`.
    @test be isa AbstractObservationErrorModel
    @test be isa AbstractObservationModel
    @test implements_observation_interface(
        be; y_t = (y = missing, N = 20),
        Y_t = fill(0.3, 10)
    )

    # Scalar N (in the data) is broadcast across the series; successes lie in 0..N.
    p = fill(0.3, 12)
    sim = as_turing_model(be, (y = missing, N = 20), p)().y_t
    @test length(sim) == length(p)
    @test all(x -> 0 <= x <= 20, sim)
    # Conditioning on the data returns it.
    @test as_turing_model(be, (y = sim, N = 20), p)().y_t == sim

    # A per-time-point trials vector in the data is honoured.
    Nvec = collect(5:14)            # length 10
    simv = as_turing_model(be, (y = missing, N = Nvec), fill(0.8, 10))().y_t
    @test all(i -> 0 <= simv[i] <= Nvec[i], eachindex(simv))
    @test as_turing_model(be, (y = simv, N = Nvec), fill(0.8, 10))().y_t == simv

    # `y_t` must be a NamedTuple carrying `N`: a plain vector or a NamedTuple
    # without `N` is rejected.
    @test_throws Exception as_turing_model(be, fill(3, 10), fill(0.3, 10))()
    @test_throws Exception as_turing_model(be, (y = missing,), fill(0.3, 10))()

    # A trials vector whose length does not match the series is rejected.
    @test_throws Exception as_turing_model(
        be, (y = missing, N = [5, 5, 5]),
        fill(0.2, 10)
    )()

    # The success probability is clamped away from 0/1 (no degenerate likelihood).
    edge = as_turing_model(be, (y = missing, N = 8), [0.0, 1.0, 0.5, 0.5])().y_t
    @test all(x -> 0 <= x <= 8, edge)
end

@testitem "define_y_t unpacks counts for vector or NamedTuple data" begin
    using ComposableTuringIDModels
    Y_t = fill(10.0, 5)
    # Plain vector passes through.
    @test define_y_t(PoissonError(), [1, 2, 3, 4, 5], Y_t) == [1, 2, 3, 4, 5]
    # NamedTuple: the `y` field is unpacked.
    @test define_y_t(PoissonError(), (y = [1, 2, 3, 4, 5],), Y_t) == [1, 2, 3, 4, 5]
    # `missing` (plain or in the `y` field) becomes a length-Y_t missing vector.
    @test all(ismissing, define_y_t(PoissonError(), missing, Y_t))
    @test length(define_y_t(PoissonError(), missing, Y_t)) == 5
    @test all(ismissing, define_y_t(PoissonError(), (y = missing,), Y_t))
    # BinomialError shares the default unpacking for its `y` field.
    @test define_y_t(BinomialError(), (y = [3, 4], N = 10), fill(0.5, 2)) == [3, 4]

    # A `MissingObservations` carrier (what `concrete_observations` returns for
    # a partially-missing vector) is rebuilt into the ragged vector, at the top
    # level and unpacked from a `y` field alike.
    using ComposableTuringIDModels: MissingObservations
    carrier = MissingObservations([1, 0, 3, 4, 5], Bool[1, 0, 1, 1, 1])
    rebuilt = define_y_t(PoissonError(), carrier, Y_t)
    @test isequal(rebuilt, [1, missing, 3, 4, 5])
    @test isequal(
        define_y_t(PoissonError(), (y = carrier,), Y_t),
        [1, missing, 3, 4, 5]
    )
end

@testitem "LatentDelay shortens expectations and wraps an error model" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(12)
    obs = LatentDelay(PoissonError(), truncated(Normal(5.0, 2.0), 0.0, Inf))
    Y_t = fill(10.0, 30)
    sim = as_turing_model(obs, missing, Y_t)().y_t
    @test length(sim) == length(Y_t)
    # LatentDelay deliberately leaves the head of the series unobserved
    # (partially observed data), so only the non-missing tail is filled.
    observed = filter(!ismissing, sim)
    @test !isempty(observed)
    @test all(>=(0), observed)
end

@testitem "safe count distributions tolerate very large means" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(13)
    # exp(48) overflows Int; the safe samplers must not throw.
    bigλ = exp(48.0)
    @test rand(SafePoisson(bigλ)) >= 0
    σ² = bigλ + 0.05 * bigλ^2
    p = bigλ / σ²
    r = bigλ * p / (1 - p)
    @test rand(SafeNegativeBinomial(r, p)) >= 0
end

@testitem "SafeNegativeBinomial rejects an invalid shape or success probability" begin
    using ComposableTuringIDModels, Distributions, Random

    # Construction itself does not validate: `r`/`p` may transiently be
    # out-of-domain Dual numbers on the AD path (e.g. built from a sampled
    # cluster factor while evaluating `logpdf`), and validating eagerly there
    # would throw on every such draw. `Base.rand` is the point where an
    # invalid value must be caught, before it reaches the opaque `sqrt`/`Gamma`
    # error deep inside sampling.
    @test SafeNegativeBinomial(-1.0, 0.5) isa SafeNegativeBinomial

    rng = Random.default_rng()

    # A negative or zero `r` is rejected at `rand`, with a clear `DomainError`
    # naming `r`.
    err = try
        rand(rng, SafeNegativeBinomial(-1.0, 0.5))
        nothing
    catch e
        e
    end
    @test err isa DomainError
    @test err.val == -1.0
    @test occursin("r (the shape", err.msg)
    @test occursin("SafeNegativeBinomial", err.msg)

    err0 = try
        rand(rng, SafeNegativeBinomial(0.0, 0.5))
        nothing
    catch e
        e
    end
    @test err0 isa DomainError
    @test err0.val == 0.0

    # `p` must lie in `(0, 1]`.
    errp = try
        rand(rng, SafeNegativeBinomial(1.0, 0.0))
        nothing
    catch e
        e
    end
    @test errp isa DomainError
    @test errp.val == 0.0
    @test occursin("p (the success", errp.msg)

    errp2 = try
        rand(rng, SafeNegativeBinomial(1.0, 1.5))
        nothing
    catch e
        e
    end
    @test errp2 isa DomainError
    @test errp2.val == 1.5

    # `p == 1` is the valid degenerate boundary and must still sample.
    @test rand(rng, SafeNegativeBinomial(1.0, 1.0)) == 0
end

@testitem "NegativeBinomialMeanClust propagates an invalid cluster factor" begin
    using ComposableTuringIDModels, Distributions, Random

    # `α <= 0` drives a non-positive `r = 1/α` through to `SafeNegativeBinomial`,
    # which must still raise a named `DomainError` when sampled, rather than
    # failing later with an opaque error.
    err = try
        rand(Random.default_rng(), NegativeBinomialMeanClust(10.0, -0.1))
        nothing
    catch e
        e
    end
    @test err isa DomainError
end

@testitem "SafePoisson rejects a negative mean" begin
    using ComposableTuringIDModels, Distributions, Random

    # Construction does not validate, for the same AD-safety reason as
    # `SafeNegativeBinomial` above.
    @test SafePoisson(-1.0) isa SafePoisson

    rng = Random.default_rng()
    err = try
        rand(rng, SafePoisson(-1.0))
        nothing
    catch e
        e
    end
    @test err isa DomainError
    @test err.val == -1.0
    @test occursin("SafePoisson", err.msg)

    # `λ == 0` is the valid degenerate boundary and must still sample.
    @test rand(rng, SafePoisson(0.0)) == 0
end

@testitem "a gap in fitted data is marginalised out, not sampled" begin
    using ComposableTuringIDModels, Distributions, Random
    using ComposableTuringIDModels: MissingObservations
    using DynamicPPL: DynamicPPL, VarInfo, logjoint
    Random.seed!(45)
    model = IDModel(
        DirectInfections(;
            Z = RandomWalk(), initialisation = Normal(log(50), 0.2)
        ),
        NegativeBinomialError()
    )
    n = 20
    y_full = fill(50, n)
    y_gap = Vector{Union{Missing, Int}}(y_full)
    y_gap[[3, 8, 14]] .= missing

    vi_full = VarInfo(as_turing_model(model, y_full, n))
    vi_gap = VarInfo(as_turing_model(model, y_gap, n))
    # Punching holes in the data adds no parameters: an absent entry leaves the
    # likelihood rather than becoming a latent to sample.
    @test collect(keys(vi_gap)) == collect(keys(vi_full))
    @test length(vi_gap[:]) == length(vi_full[:])
    @test !any(vn -> occursin("y_t", string(vn)), keys(vi_gap))

    # The returned series marks the gaps rather than filling them, and passes
    # the observed entries through as given.
    out = as_turing_model(model, y_gap, n)()
    @test all(ismissing, out.generated_y_t[[3, 8, 14]])
    @test out.generated_y_t[[1, 2, 20]] == y_full[[1, 2, 20]]

    # Nothing at a gap reaches the likelihood: two carriers differing only in
    # the placeholder they hold at the absent entries score identically, while
    # the fully observed series scores differently (the observed entries do
    # still contribute).
    present = trues(n)
    present[[3, 8, 14]] .= false
    carrier = MissingObservations(copy(y_full), present)
    shifted = MissingObservations(
        map((v, p) -> p ? v : v + 1000, y_full, present), present
    )
    m_full = as_turing_model(model, y_full, n)
    vi = VarInfo(m_full)
    @test logjoint(as_turing_model(model, carrier, n), vi) ==
        logjoint(as_turing_model(model, shifted, n), vi)
    @test logjoint(m_full, vi) != logjoint(as_turing_model(model, carrier, n), vi)
end

@testitem "NUTS samples a discrete stream with gaps" tags = [:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(46)
    # A count observation has no bijector to link, so a sampled gap would make
    # this model unsamplable by NUTS. Marginalising the gaps out leaves only
    # the continuous latents.
    model = IDModel(
        DirectInfections(;
            Z = RandomWalk(), initialisation = Normal(log(50), 0.2)
        ),
        NegativeBinomialError()
    )
    n = 20
    y = Vector{Union{Missing, Int}}(fill(50, n))
    y[[4, 9, 15]] .= missing
    chain = sample(
        as_turing_model(model, y, n), NUTS(), 20; progress = false
    )
    @test size(chain, 1) == 20
    @test !any(k -> occursin("y_t", string(k)), keys(chain))
end

@testitem "a gap is filled predictively after fitting, not during" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    using DynamicPPL: VarInfo
    Random.seed!(47)
    model = IDModel(
        DirectInfections(;
            Z = RandomWalk(), initialisation = Normal(log(50), 0.2)
        ),
        PoissonError()
    )
    n = 12
    y = Vector{Union{Missing, Int}}(fill(50, n))
    y[[3, 7]] .= missing

    # Fitting carries no `y_t` at all, gap or otherwise.
    fit = as_turing_model(model, y, n)
    @test !any(vn -> occursin("y_t", string(vn)), keys(VarInfo(fit)))
    chain = sample(fit, Prior(), 20; progress = false)

    # Replaying the posterior through a wholly unobserved series generates
    # every point, the gaps among them, from the error model at the expected
    # values that replay produces.
    pred = predict(as_turing_model(model, fill(missing, n), n), chain)
    drawn = reduce(hcat, (vec(pred[@varname(y_t[i])]) for i in 1:n))
    @test size(drawn) == (20, n)
    @test all(x -> x isa Integer && x >= 0, drawn)

    # The same replay recovers the expected series the draws are taken at.
    gens = vec(returned(as_turing_model(model, y, n), chain))
    @test length(gens[1].expected_y_t) == n
    @test all(>=(0), gens[1].expected_y_t)
end
