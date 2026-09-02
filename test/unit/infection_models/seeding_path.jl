# The renewal seeding window: the level slot it is decayed from, the
# `SeedingPath` that estimates it instead, and its exposure as `I_seed`.

@testitem "the default seeding window decays a level and is returned as I_seed" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: renewal_init_window
    using DynamicPPL: fix
    gen_int = [0.2, 0.3, 0.5]
    renewal = Renewal(;
        generation_time = gen_int, rt = FixedIntercept(log(1.5)),
        initialisation = Normal()
    )
    out = fix(as_turing_model(renewal, 10), (init_incidence = log(10.0),))()
    r = R_to_r(1.5, gen_int)
    @test out.I_seed ≈ renewal_init_window(renewal.recurrent_step, 10.0, r, 3)
    # The window ends at the level it was seeded from, and is not part of `I_t`.
    @test last(out.I_seed) ≈ 10.0
    @test length(out.I_t) == 10
end

@testitem "a SeedingPath initialisation is used as the seeding window itself" begin
    using ComposableTuringIDModels, Distributions
    using LinearAlgebra: dot
    gen_int = [0.2, 0.3, 0.5]
    # Every component is deterministic, so the window and the first scanned
    # step are exact: a flat seeding path of 10 infections and `R_t` = 1.5.
    seeded = Renewal(;
        generation_time = gen_int, rt = FixedIntercept(log(1.5)),
        initialisation = SeedingPath(FixedIntercept(log(10.0)))
    )
    out = as_turing_model(seeded, 10)()
    @test out.I_seed ≈ fill(10.0, 3)
    # No geometric decay is applied: the drawn path is the window the scan
    # convolves, so the first infection is `R_1` times that convolution.
    @test out.I_t[1] ≈ 1.5 * dot(fill(10.0, 3), reverse(gen_int))
    @test length(out.I_t) == 10
end

@testitem "a SeedingPath is drawn at the window length, not at length 1" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(263)
    gen_int = [0.1, 0.2, 0.3, 0.4]
    seeded = Renewal(;
        generation_time = gen_int, rt = RandomWalk(),
        initialisation = SeedingPath(RandomWalk(; init = Normal(log(50), 0.5)))
    )
    out = as_turing_model(seeded, 15)()
    # The whole run-up is estimated: one value per generation-interval lag,
    # varying rather than collapsed to a single level.
    @test length(out.I_seed) == length(gen_int)
    @test all(>(0), out.I_seed)
    @test length(unique(out.I_seed)) == length(gen_int)
end

@testitem "a stratified SeedingPath gives one seeding window per stratum" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(264)
    gen_int = [0.2, 0.3, 0.5]
    seeded = Renewal(;
        generation_time = gen_int,
        rt = Stratify(RandomWalk(), Hierarchy()),
        initialisation = SeedingPath(Stratify(RandomWalk(), Hierarchy()))
    )
    out = as_turing_model(seeded, (3, 12))()
    @test size(out.I_seed) == (3, 3)
    @test size(out.I_t) == (3, 12)
    @test all(isfinite, out.I_t)
end

@testitem "a SeedingPath composes with a renewal modifier" begin
    using ComposableTuringIDModels, Distributions
    gen_int = [0.2, 0.3, 0.5]
    # The modifier's substate is built from the drawn window, so susceptible
    # depletion starts from a pool matching the seeded incidence.
    seeded = Renewal(
        gen_int, SusceptibleDepletion(100.0);
        rt = FixedIntercept(log(1.5)),
        initialisation = SeedingPath(FixedIntercept(log(10.0)))
    )
    out = as_turing_model(seeded, 10)()
    @test out.I_seed ≈ fill(10.0, 3)
    # The pool is full at the first step, so depletion shows later in the run.
    plain = Renewal(;
        generation_time = gen_int, rt = FixedIntercept(log(1.5)),
        initialisation = SeedingPath(FixedIntercept(log(10.0)))
    )
    @test last(out.I_t) < last(as_turing_model(plain, 10)().I_t)
end

@testitem "IDModel exposes the renewal seeding window" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(265)
    model = IDModel(
        Renewal(;
            generation_time = [0.2, 0.3, 0.5], rt = RandomWalk(),
            initialisation = SeedingPath(
                RandomWalk(; init = Normal(log(50), 0.5))
            )
        ),
        PoissonError()
    )
    sim = as_turing_model(model, missing, 12)()
    @test haskey(sim, :I_seed)
    @test length(sim.I_seed) == 3
    # An infection model with no seeding window returns the tuple it always did.
    direct = IDModel(
        DirectInfections(; Z = RandomWalk(), initialisation = Normal()),
        PoissonError()
    )
    @test keys(as_turing_model(direct, missing, 12)()) ==
        (:generated_y_t, :expected_y_t, :I_t, :Z_t)
end

@testitem "SeedingPath is rejected outside a renewal initialisation slot" begin
    using ComposableTuringIDModels, Distributions
    sp = SeedingPath(RandomWalk())
    @test_throws "only" as_turing_model(sp, 5)
    @test_throws "only" as_turing_model(sp, (2, 5))
    # A bare distribution widens to a constant path, as every PATH slot does.
    @test SeedingPath(Normal()).model isa Intercept
end

@testitem "a level initialisation rejects a draw that is not a single value" begin
    using ComposableTuringIDModels
    using ComposableTuringIDModels: _seed
    @test _seed([2.0], 5) == 2.0
    @test _seed(2.0, 5) == 2.0
    # The slot is a level, so a longer draw is an error naming the alternative
    # rather than a silently collapsed path.
    @test_throws "SeedingPath" _seed([1.0, 2.0], 5)
end
