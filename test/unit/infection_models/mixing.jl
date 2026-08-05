# Tests for the coupling seam: `renewal_pressure`'s methods (each checked
# against a hand-written sum, not a re-implementation of the same broadcast),
# a fixed off-diagonal mixing matrix moving incidence between strata, the
# drawn `Gravity` operator (and the plain `gravity` function underneath it),
# and `pairwise_gen_int`.

@testitem "renewal_pressure: single series, UniformScaling" begin
    using ComposableTuringIDModels: renewal_pressure
    using LinearAlgebra: I
    g = [0.2, 0.3, 0.5]
    window = [10.0, 20.0, 30.0]
    hand = sum(window[i] * g[i] for i in eachindex(g))
    @test renewal_pressure(I, g, window) ≈ hand
end

@testitem "renewal_pressure: strata window, shared generation interval" begin
    using ComposableTuringIDModels: renewal_pressure
    using LinearAlgebra: I
    g = [0.2, 0.3, 0.5]
    window = [10.0 20.0 30.0; 1.0 2.0 3.0]   # 2 strata x 3 lags
    hand = [sum(window[s, i] * g[i] for i in eachindex(g)) for s in 1:2]
    @test renewal_pressure(I, g, window) ≈ hand
end

@testitem "renewal_pressure: strata window, one interval per stratum" begin
    using ComposableTuringIDModels: renewal_pressure
    using LinearAlgebra: I
    g = [0.2 0.3 0.5; 0.5 0.3 0.2]            # 2 strata x 3 lags
    window = [10.0 20.0 30.0; 1.0 2.0 3.0]
    hand = [sum(window[s, i] * g[s, i] for i in 1:3) for s in 1:2]
    @test renewal_pressure(I, g, window) ≈ hand
end

@testitem "renewal_pressure: a mixing matrix redistributes histories" begin
    using ComposableTuringIDModels: renewal_pressure
    g = [0.2, 0.3, 0.5]
    window = [10.0 20.0 30.0; 1.0 2.0 3.0]
    K = [0.9 0.1; 0.2 0.8]
    conv = [sum(window[s, i] * g[i] for i in eachindex(g)) for s in 1:2]
    hand = [sum(K[gg, h] * conv[h] for h in 1:2) for gg in 1:2]
    @test renewal_pressure(K, g, window) ≈ hand
end

@testitem "renewal_pressure: a per-pair strata x strata x lags interval" begin
    using ComposableTuringIDModels: renewal_pressure
    n_strata, n_lags = 2, 3
    window = [10.0 20.0 30.0; 1.0 2.0 3.0]   # oldest -> newest columns
    K = zeros(n_strata, n_strata, n_lags)
    K[1, 1, :] = [0.5, 0.3, 0.2]
    K[1, 2, :] = [0.10, 0.05, 0.05]
    K[2, 1, :] = [0.00, 0.02, 0.03]
    K[2, 2, :] = [0.6, 0.2, 0.1]
    n_time = size(window, 2)
    hand = [sum(K[gg, h, i] * window[h, n_time - i + 1]
            for h in 1:n_strata, i in 1:n_lags)
            for gg in 1:n_strata]
    # The reversed generation interval `g` is unused: the array carries the
    # per-pair intervals itself.
    @test renewal_pressure(K, [0.0], window) ≈ hand
end

@testitem "renewal_pressure: a mixing matrix AND a per-stratum interval" begin
    using ComposableTuringIDModels: renewal_pressure
    # `K` is checked with a shared `g`, and a per-stratum `g` is checked with
    # `I`, but never both non-identity together, so the general path (a
    # mixing matrix on top of per-stratum convolved histories) is unchecked.
    g = [0.2 0.3 0.5; 0.5 0.3 0.2]            # 2 strata x 3 lags
    window = [10.0 20.0 30.0; 1.0 2.0 3.0]
    K = [0.9 0.1; 0.2 0.8]
    conv = [sum(window[s, i] * g[s, i] for i in 1:3) for s in 1:2]
    hand = [sum(K[gg, h] * conv[h] for h in 1:2) for gg in 1:2]
    @test renewal_pressure(K, g, window) ≈ hand
end

@testitem "an off-diagonal mixing matrix moves incidence between strata" begin
    using ComposableTuringIDModels, Distributions
    using LinearAlgebra: I
    gen_int = [0.2, 0.3, 0.5]
    n_strata, n_time = 2, 25
    K = [0.9 0.3; 0.1 0.7]   # off-diagonal: cross-stratum pressure
    rt = Stratify(FixedIntercept(log(1.2)), FixedIntercept(0.0))
    seeds = [Dirac(log(5.0)), Dirac(log(0.05))]   # stratum 2 barely seeded

    uncoupled = Renewal(; generation_time = gen_int, rt = rt,
        initialisation = seeds)
    coupled = Renewal(; generation_time = gen_int, rt = rt,
        initialisation = seeds, mixing = K)
    plain = Renewal(; generation_time = gen_int, rt = rt,
        initialisation = seeds, mixing = I)

    I_u = as_turing_model(uncoupled, (n_strata, n_time))().I_t
    I_c = as_turing_model(coupled, (n_strata, n_time))().I_t
    I_p = as_turing_model(plain, (n_strata, n_time))().I_t

    # `mixing = I` explicitly matches the uncoupled default exactly.
    @test I_p ≈ I_u
    # Only cross-stratum coupling can lift the near-zero-seeded stratum.
    @test all(I_c[2, :] .> I_u[2, :])
end

@testitem "gravity builds an asymmetric matrix from pop and distance" begin
    using ComposableTuringIDModels
    pop = [1e6, 2e5]
    dist = [0.0 50.0; 50.0 0.0]
    K = gravity(pop, dist; α = 1.0, β = 1.0, γ = 2.0, within = 1.0)
    @test K[1, 1] == 1.0 == K[2, 2]
    @test K[1, 2] ≈ pop[1]^1.0 * pop[2]^1.0 / dist[1, 2]^2.0
    @test K[2, 1] ≈ pop[2]^1.0 * pop[1]^1.0 / dist[2, 1]^2.0
    within = gravity(pop, dist; within = 0.5)
    @test within[1, 1] == 0.5 == within[2, 2]
end

@testitem "gravity rejects a dist matrix that doesn't match pop" begin
    using ComposableTuringIDModels
    @test_throws AssertionError gravity(
        [1e6, 2e5, 3e5], [0.0 1.0; 1.0 0.0])
end

@testitem "a drawn Gravity mixing samples its exponents" begin
    using ComposableTuringIDModels, Distributions, Random
    using LinearAlgebra: diag
    pop = [1e6, 2e5, 5e5]
    dist = [0.0 50.0 80.0; 50.0 0.0 30.0; 80.0 30.0 0.0]
    g = Gravity(pop, dist)
    Random.seed!(704)
    K1 = as_turing_model(g, (3, 10))()
    Random.seed!(705)
    K2 = as_turing_model(g, (3, 10))()
    @test size(K1) == (3, 3)
    @test K1 != K2   # different seeds draw different exponents
    @test all(==(1.0), diag(K1))       # default `within`
end

@testitem "fixing Gravity's exponents matches the plain gravity matrix" begin
    using ComposableTuringIDModels, Distributions
    pop = [1e6, 2e5]
    dist = [0.0 40.0; 40.0 0.0]
    fixed_K = gravity(pop, dist; α = 0.5, β = 0.8, γ = 1.5)
    g = Gravity(pop, dist; α = FixedIntercept(0.5), β = FixedIntercept(0.8),
        γ = FixedIntercept(1.5))
    drawn_K = as_turing_model(g, (2, 10))()
    @test drawn_K ≈ fixed_K
end

@testitem "a Renewal with a Gravity mixing resolves through the seam" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(706)
    pop = [1e6, 2e5]
    dist = [0.0 40.0; 40.0 0.0]
    gen_int = [0.2, 0.3, 0.5]
    grav = Gravity(pop, dist; α = FixedIntercept(0.5),
        β = FixedIntercept(0.5), γ = FixedIntercept(2.0))
    model = Renewal(; generation_time = gen_int,
        rt = Stratify(RandomWalk(), Hierarchy()),
        initialisation = IID(Normal(log(20.0), 0.5)), mixing = grav)
    out = as_turing_model(model, (2, 20))()
    @test size(out.I_t) == (2, 20)
    @test all(isfinite, out.I_t)
    @test all(>=(0), out.I_t)
end

@testitem "a Gravity mixing's inferred exponent is named under mixing" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(707)
    pop = [1e6, 2e5]
    dist = [0.0 40.0; 40.0 0.0]
    gen_int = [0.2, 0.3, 0.5]
    grav = Gravity(pop, dist; α = Normal(0, 0.5), β = FixedIntercept(0.5),
        γ = FixedIntercept(2.0))
    model = Renewal(; generation_time = gen_int,
        rt = Stratify(RandomWalk(), FixedIntercept(0.0)),
        initialisation = Normal(), mixing = grav)
    mdl = as_turing_model(model, (2, 20))
    names = string.(collect(keys(rand(mdl))))
    # `scan_step ~ ...` composes flat (no `scan_step.` prefix), so the drawn
    # exponent's own prefix (`mixing.α`) reaches the top level unchanged.
    @test any(startswith(n, "mixing.") for n in names)
    @test !any(startswith(n, "scan_step.") for n in names)
end

@testitem "an unresolved MixingStep says so instead of a MethodError" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: MixingStep, _renewal_init_state
    step = MixingStep(reverse([0.2, 0.3, 0.5]),
        Gravity([1e6, 2e5], [0.0 1.0; 1.0 0.0]))
    @test_throws "has no scan interface" step([1.0 2.0; 1.0 2.0], [1.0, 1.0])
    @test_throws "has no scan interface" _renewal_init_state(
        step, [1.0, 1.0], [0.1, 0.1], 3)
end

@testitem "pairwise_gen_int builds a strata x strata x lags mixing array" begin
    using ComposableTuringIDModels
    K = [0.9 0.1; 0.05 0.95]
    g = [0.2, 0.3, 0.5]
    P = pairwise_gen_int(K, g)
    @test size(P) == (2, 2, 3)
    @test P[1, 1, :] ≈ K[1, 1] .* g
    @test P[2, 1, :] ≈ K[2, 1] .* g

    # Per-pair intervals: a between-stratum transmission can carry a longer
    # effective interval than a within-stratum one.
    G = zeros(2, 2, 4)
    G[1, 1, :] = [0.5, 0.3, 0.2, 0.0]
    G[1, 2, :] = [0.1, 0.2, 0.3, 0.4]
    G[2, 1, :] = [0.1, 0.2, 0.3, 0.4]
    G[2, 2, :] = [0.5, 0.3, 0.2, 0.0]
    P2 = pairwise_gen_int(K, G)
    @test P2[1, 2, :] ≈ K[1, 2] .* G[1, 2, :]

    @test_throws "non-negative" pairwise_gen_int(K, [-0.1, 0.5, 0.6])
    @test_throws "sum to 1" pairwise_gen_int(K, [0.2, 0.2, 0.2])
end

@testitem "pairwise_gen_int plugs into a Renewal's mixing slot" begin
    using ComposableTuringIDModels, Distributions, Random
    using LinearAlgebra: I
    # `pairwise_gen_int` is only unit-tested as a bare function elsewhere; a
    # `strata x strata x lags` array built from it is never actually handed
    # to `Renewal(...; mixing = ...)`. With `K = I` the per-pair intervals
    # degenerate to the shared interval on the diagonal alone, so this must
    # reproduce the uncoupled (`mixing = I`) model exactly — the identity
    # that pins the array's lag convention, worth a regression test rather
    # than only something checked by hand.
    gen_int = [0.2, 0.3, 0.5]
    n_strata, n_time = 2, 20
    rt = Stratify(RandomWalk(), Hierarchy())
    seeds = IID(Normal(log(20.0), 0.5))

    identity_pairwise = pairwise_gen_int([1.0 0.0; 0.0 1.0], gen_int)
    coupled_identity = Renewal(; generation_time = gen_int, rt = rt,
        initialisation = seeds, mixing = identity_pairwise)
    uncoupled = Renewal(; generation_time = gen_int, rt = rt,
        initialisation = seeds, mixing = I)

    Random.seed!(900)
    I_pairwise = as_turing_model(coupled_identity, (n_strata, n_time))().I_t
    Random.seed!(900)
    I_uncoupled = as_turing_model(uncoupled, (n_strata, n_time))().I_t
    @test I_pairwise ≈ I_uncoupled

    # A genuinely off-diagonal pairwise interval still runs and differs.
    K = [0.9 0.1; 0.2 0.8]
    coupled = Renewal(; generation_time = gen_int, rt = rt,
        initialisation = seeds, mixing = pairwise_gen_int(K, gen_int))
    Random.seed!(900)
    I_coupled = as_turing_model(coupled, (n_strata, n_time))().I_t
    @test all(isfinite, I_coupled)
    @test I_coupled != I_uncoupled
end
