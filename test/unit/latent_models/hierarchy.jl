# Tests for the numeric `Hierarchy` partial-pooling construct: it returns a
# length-`n_groups` numeric per-group level path (shared mean + group deviations),
# with the cross-group relationship supplied through the prior interface
# (iid-Normal ⇒ classic exchangeable pooling, RandomWalk ⇒ correlated neighbours).

@testitem "Hierarchy constructs and is an AbstractLatentModel" begin
    using ComposableTuringIDModels, Distributions
    h = Hierarchy(Normal(), IID(Normal(0.0, 1.0)))
    @test h isa Hierarchy
    @test h isa AbstractLatentModel
    # Both slots are coerced to the prior interface.
    @test h.mean isa AbstractPriorModel
    @test h.across isa IID
    # Keyword constructor; across defaults to an iid Normal.
    hk = Hierarchy()
    @test hk isa Hierarchy
    @test hk.across isa IID
    # A bare Distribution in `mean` is coerced via as_prior.
    @test Hierarchy(; mean = Normal(0, 2)).mean isa BroadcastPrior
end

@testitem "Hierarchy returns a numeric length-n_groups path" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(381)
    h = Hierarchy(Normal(), IID(Normal(0.0, 1.0)))
    # n_groups is NOT stored on the struct; it is read at build time.
    @test !(:n_groups in fieldnames(typeof(h)))
    for G in (2, 4, 7)
        vals = as_turing_model(h, G)()
        @test length(vals) == G
        @test all(x -> x isa Real, vals)      # numeric, not model variants
    end
end

@testitem "iid-Normal across recovers classic (exchangeable) partial pooling" begin
    using ComposableTuringIDModels, Distributions, Random, Statistics
    Random.seed!(383)
    # Zero shared mean so the group values are the (exchangeable) deviations.
    h = Hierarchy(Dirac(0.0), IID(Normal(0.0, 1.0)))
    G = 6
    reps = 3000
    M = reduce(hcat, [as_turing_model(h, G)() for _ in 1:reps])   # G × reps
    # Exchangeable draws ⇒ near-zero correlation between neighbouring groups.
    neighbour_cors = [cor(M[g, :], M[g + 1, :]) for g in 1:(G - 1)]
    @test mean(abs.(neighbour_cors)) < 0.1
    @test abs(mean(M)) < 0.1
    @test isapprox(std(vec(M)), 1.0; atol = 0.1)
end

@testitem "RandomWalk across gives correlated neighbouring group effects" begin
    using ComposableTuringIDModels, Distributions, Random, Statistics
    Random.seed!(384)
    h_rw = Hierarchy(Dirac(0.0), RandomWalk())
    h_iid = Hierarchy(Dirac(0.0), IID(Normal(0.0, 1.0)))
    G = 6
    reps = 3000
    rw = reduce(hcat, [as_turing_model(h_rw, G)() for _ in 1:reps])
    iid = reduce(hcat, [as_turing_model(h_iid, G)() for _ in 1:reps])
    rw_neighbour = mean([cor(rw[g, :], rw[g + 1, :]) for g in 1:(G - 1)])
    iid_neighbour = mean([cor(iid[g, :], iid[g + 1, :]) for g in 1:(G - 1)])
    # The random walk relates neighbours: strong positive lag-1 correlation,
    # well above the (near-zero) iid case.
    @test rw_neighbour > 0.5
    @test rw_neighbour > iid_neighbour + 0.4
end

@testitem "Hierarchy conforms to the latent interface and composes numerically" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(385)
    h = Hierarchy(Normal(), IID(Normal(0.0, 1.0)))
    @test implements_latent_interface(h; n = 4)
    @test !implements_observation_interface(h)
    @test !implements_infection_interface(h)
    # Because it returns numeric values, it threads straight into an infection
    # model's latent slot (no group-axis contract change).
    inf = DirectInfections(; Z = h, initialisation = Normal())
    out = as_turing_model(inf, 5)()
    @test length(out.I_t) == 5
    @test all(x -> x isa Real, out.Z_t)
end

@testitem "Hierarchy drives a per-group quantity in a stacked composed model" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random, Statistics
    using Turing: to_submodel, returned
    Random.seed!(388)
    # A FULL composed model (IDModel: infection + observation) whose per-group
    # reporting level is partially pooled by a Hierarchy. The group dimension is
    # read from the data matrix (its columns), NOT stored on any component.
    idmodel = IDModel(
        DirectInfections(;
            Z = RandomWalk(), initialisation = Normal(log(50.0), 0.2)),
        PoissonError())
    hierarchy = Hierarchy(;
        mean = Normal(0.0, 0.5), across = IID(Normal(0.0, 0.5)))

    @model function grouped_epidemic(idmodel, hierarchy, Y)
        n_time, n_groups = size(Y)
        # Namespace the group prior so its IID innovations do not collide with the
        # infection process's own `ϵ_t` under the prefix-off convention: this is
        # the group-requirement threading made explicit.
        group_levels ~ to_submodel(
            as_turing_model(PrefixLatentModel(hierarchy, "groups"), n_groups),
            false)
        infections ~ to_submodel(
            as_turing_model(idmodel.infection_model, n_time), false)
        I_t = infections.I_t
        ys = Vector{Any}(undef, n_groups)
        for g in 1:n_groups
            expected_g = exp(group_levels[g]) .* I_t
            og = PrefixObservationModel(idmodel.observation_model, "group$g")
            y_g ~ to_submodel(as_turing_model(og, Y[:, g], expected_g), false)
            ys[g] = y_g
        end
        return (; I_t, group_levels, y = ys)
    end

    n_time, n_groups = 20, 4
    # n_groups threads from the number of columns of the (missing) data matrix.
    Ymiss = Matrix{Union{Missing, Float64}}(missing, n_time, n_groups)
    sim = grouped_epidemic(idmodel, hierarchy, Ymiss)()
    @test length(sim.I_t) == n_time
    @test length(sim.group_levels) == n_groups
    Ydata = reduce(hcat, [Int.(sim.y[g].y_t) for g in 1:n_groups])
    @test size(Ydata) == (n_time, n_groups)

    posterior = grouped_epidemic(idmodel, hierarchy, Float64.(Ydata))
    chain = sample(posterior, NUTS(0.8; adtype = Turing.AutoForwardDiff()), 80;
        progress = false)
    @test size(chain, 1) == 80
    # The per-group levels are recovered as generated quantities.
    draws = reduce(hcat, [g.group_levels for g in vec(returned(posterior, chain))])
    @test size(draws) == (n_groups, 80)
    post_mean = vec(mean(draws; dims = 2))
    @test cor(sim.group_levels, post_mean) > 0.5
end

@testitem "Hierarchy plugs into a component prior slot (Ascertainment)" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(390)
    # Because `Hierarchy` is an `AbstractPriorModel`, it drops straight into any
    # component prior slot through the priors weave — no bespoke grouped model.
    # A prior slot passes its own length `n`, so when that length is the group
    # dimension the pooling is across groups. Here the `Ascertainment` slot passes
    # `n = length(Y_t)`: with one observation per group, each group gets its own
    # partially-pooled ascertainment intercept.
    h = Hierarchy(; mean = Normal(0.0, 0.3), across = IID(Normal(0.0, 0.4)))
    asc = Ascertainment(PoissonError(), h)
    @test asc isa AbstractObservationModel
    @test implements_observation_interface(asc; Y_t = fill(300.0, 6))
    G = 6
    draw = as_turing_model(asc, missing, fill(300.0, G))()
    # One observation per group; the per-group expected carries the pooled effect.
    @test length(draw.y_t) == G
    @test length(draw.expected) == G
end

@testitem "Hierarchy as an Ascertainment prior recovers per-group intercepts" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random, Statistics
    using Turing: returned
    Random.seed!(391)
    # A per-group ascertainment intercept partially pooled by a `Hierarchy` passed
    # straight into the `Ascertainment` prior slot. Each group contributes one
    # Poisson count whose expected value is a shared baseline scaled by the pooled
    # per-group effect `exp(ℓ_g)`. Sampling recovers the per-group intercepts.
    h = Hierarchy(; mean = Normal(0.0, 0.3), across = IID(Normal(0.0, 0.4)))
    asc = Ascertainment(PoissonError(), h)
    G = 8
    base = fill(300.0, G)
    sim = as_turing_model(asc, missing, base)()
    ydata = Int.(sim.y_t)
    true_levels = log.(sim.expected ./ base)
    post = as_turing_model(asc, Float64.(ydata), base)
    chain = sample(post, NUTS(0.8; adtype = Turing.AutoForwardDiff()), 200;
        progress = false)
    @test size(chain, 1) == 200
    # Per-group intercepts recovered from the returned expected series.
    rets = vec(returned(post, chain))
    expected_draws = reduce(hcat, [r.expected for r in rets])
    post_mean = vec(mean(log.(expected_draws ./ base); dims = 2))
    @test cor(true_levels, post_mean) > 0.8
end

@testitem "a partially-pooled model samples under NUTS (ForwardDiff)" tags=[:sample] begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    using DynamicPPL: to_submodel
    Random.seed!(386)
    # Each group has its own partially-pooled level mapped to Poisson counts.
    # n_groups is read from the data (the length of the observed group vector).
    h = Hierarchy(Normal(), IID(Normal(0.0, 1.0)))

    @model function pooled_counts(h, y, n_groups)
        levels ~ to_submodel(as_turing_model(h, n_groups), false)
        for g in 1:n_groups
            y[g] ~ Poisson(exp(levels[g]))
        end
        return y
    end

    G = 5
    ydata = Int.(pooled_counts(h, fill(missing, G), G)())
    chn = sample(pooled_counts(h, ydata, G),
        NUTS(0.8; adtype = Turing.AutoForwardDiff()), 60; progress = false)
    @test chn !== nothing
    @test size(chn, 1) == 60
end
