@testitem "AbstractPriorModel role: latent-as-prior, not bare distributions" begin
    using ComposableTuringIDModels, Distributions
    @test AbstractLatentModel <: AbstractPriorModel
    # Latent models satisfy the prior contract, so they are prior models.
    @test RandomWalk() isa AbstractPriorModel
    @test AR() isa AbstractPriorModel
    # A bare Distribution is not a prior model (it composes through the seam).
    @test !(Normal() isa AbstractPriorModel)
    @test implements_prior_interface(RandomWalk())
    @test !implements_prior_interface(Normal())
end

@testitem "as_turing_submodel composes a component or a raw prior" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(100)
    # The seam threads a component or a raw prior through a submodel.
    @model function use_prior(prior, n)
        θ ~ as_turing_submodel(prior, n)
        return θ
    end
    # A bare Distribution flows through the seam as a single scalar RV.
    @test use_prior(Normal(), 4)() isa Real
    # A process flows through as a length-n path; IID() gives n i.i.d. draws.
    @test length(use_prior(RandomWalk(), 6)()) == 6
    @test length(use_prior(IID(Normal()), 6)()) == 6
    @test length(use_prior([Normal(), Normal()], 2)()) == 2
end

@testitem "as_turing_submodel prefix kwarg namespaces inner variables" begin
    using ComposableTuringIDModels, Turing
    using DynamicPPL: VarInfo
    # prefix = true names the inner variables under the slot's LHS name, so a
    # RandomWalk's `rw_init`/`ϵ_t` become `z.rw_init...`/`z.ϵ_t...`.
    @model function prefixed(m, n)
        z ~ as_turing_submodel(m, n; prefix = true)
        return z
    end
    prefixed_keys = string.(keys(VarInfo(prefixed(RandomWalk(), 6))))
    @test all(startswith("z."), prefixed_keys)
    @test any(occursin("rw_init", k) for k in prefixed_keys)
    # prefix = false (the default) keeps the inner variables flat, unprefixed.
    @model function flat(m, n)
        z ~ as_turing_submodel(m, n; prefix = false)
        return z
    end
    flat_keys = string.(keys(VarInfo(flat(RandomWalk(), 6))))
    @test !any(startswith("z."), flat_keys)
    @test any(occursin("rw_init", k) for k in flat_keys)
end

@testitem "as_turing_model over a Distribution draws a single scalar" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(101)
    # The single seam draws a bare Distribution as ONE scalar RV (a constant, no
    # length-n allocation); `n` is ignored.
    @test as_turing_model(Normal(), 5)() isa Real
    @test as_turing_model(Normal(), 1)() isa Real
    # For n i.i.d. draws use the explicit IID() component instead.
    v = as_turing_model(IID(Normal()), 5)()
    @test length(v) == 5
    @test !all(==(first(v)), v)
end

@testitem "as_turing_model over a vector gives one draw per element" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(102)
    # Homogeneous vector: filldist path — two INDEPENDENT draws (not one shared
    # value repeated). This is the AR/MA per-lag coefficient semantics.
    vh = as_turing_model([Normal(), Normal()], 2)()
    @test length(vh) == 2
    @test vh[1] != vh[2]
    # Heterogeneous vector: product_distribution, element i from distribution i.
    vt = as_turing_model([Normal(0, 1), Normal(5, 0.1)], 2)()
    @test length(vt) == 2
    @test vt[2] > vt[1]                  # second is tight around 5, first around 0
    # Length must match the vector length.
    @test_throws Exception as_turing_model([Normal(), Normal()], 3)()
end

@testitem "vector damp/θ prior ⇒ independent per-lag AR/MA coefficients" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(204)
    # Two identical per-lag priors draw two INDEPENDENT coefficients — a vector
    # prior is not a single shared coefficient repeated across lags.
    d = rand(as_turing_model(
        AR(; damp = [Normal(0.0, 1.0), Normal(0.0, 1.0)], init = Normal()), 8))
    damp = reduce(vcat, [d[k] for k in keys(d) if occursin("damp_AR", string(k))])
    @test length(damp) == 2
    @test damp[1] != damp[2]
    Random.seed!(205)
    dm = rand(as_turing_model(
        MA(; θ = [Normal(0.0, 1.0), Normal(0.0, 1.0)]), 8))
    θ = reduce(vcat, [dm[k] for k in keys(dm) if occursin("θ", string(k))])
    @test length(θ) == 2
    @test θ[1] != θ[2]
end

@testitem "bare latent-model-as-prior threads under a linked log-density (#80)" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    using DynamicPPL: LogDensityFunction, VarInfo, link, getlogjoint
    import LogDensityProblems as LDP
    Random.seed!(180)
    # The exact #80 repro: a bare `AR(damp = RandomWalk())` (now a time-varying
    # coefficient path) used to sample via `rand` but ERROR as a linked
    # log-density, because the damping RandomWalk's inner `std`/`ϵ_t`/`rw_init`
    # collided with the AR innovation's. The prior slot namespaces the whole
    # submodel (prefix-on `as_turing_submodel`), so they cannot collide.
    m = as_turing_model(AR(; damp = RandomWalk()), 8)
    @test rand(m) !== nothing                        # sampled fine before too
    vi = link(VarInfo(m), m)
    ldf = LogDensityFunction(m, getlogjoint, vi)
    val = LDP.logdensity(ldf, zeros(LDP.dimension(ldf)))   # previously threw
    @test isfinite(val)
    # And it samples under NUTS end-to-end.
    chn = sample(m, NUTS(0.8; adtype = Turing.AutoForwardDiff()), 40;
        progress = false)
    @test size(chn, 1) == 40
end

@testitem "order inference: single distribution ⇒ 1, vector ⇒ length" begin
    using ComposableTuringIDModels, Distributions
    using ComposableTuringIDModels: _prior_order
    # A single distribution (or richer prior) implies order 1, a vector fixes it
    # to the vector length.
    @test _prior_order(Normal()) == 1
    @test _prior_order([Normal(), Normal()]) == 2
    @test _prior_order(RandomWalk()) == 1
    # A single (non-vector) damping distribution ⇒ order-1 AR.
    @test AR(; damp = Normal()).p == 1
    @test AR(; damp = [truncated(Normal(0, 0.05), 0, 1)]).p == 1
    # A length-k vector damp prior ⇒ order p = k, with `init` sized to match even
    # when left at its default (regression: previously threw on the init length).
    @test AR(; damp = [Normal(), Normal()]).p == 2
    @test length(as_turing_model(AR(; damp = [Normal(), Normal()]), 8)()) == 8
    @test MA(; θ = [Normal(), Normal()]).q == 2
    @test MA(; θ = Normal()).q == 1
    # A vector-of-distributions init sets the differencing order.
    @test DiffLatentModel(; model = RandomWalk(), init = [Normal(), Normal()]).d == 2
end

@testitem "raw priors are stored unwrapped in component fields" begin
    using ComposableTuringIDModels, Distributions
    # Scalar prior slots hold the bare distribution (drawn with a native tilde).
    @test NormalError(truncated(Normal(0, 1), 0, Inf)).std isa Distribution
    @test NegativeBinomialError(HalfNormal(0.1)).cluster_factor isa Distribution
    @test Intercept(; intercept = Normal()).intercept isa Distribution
    @test HierarchicalNormal(truncated(Normal(0, 1), 0, Inf)).std isa Distribution
    # Vector prior slots hold the raw vector.
    @test AR(; damp = [Normal(), Normal()],
        init = [Normal(), Normal()]).damp isa AbstractVector
    # A process prior slot holds the latent model unchanged.
    @test AR(; damp = RandomWalk()).damp isa RandomWalk
    # Positional bare-`Distribution` constructor forms build valid latents.
    for m in (AR(Normal(), Normal()),
        HierarchicalNormal(truncated(Normal(0, 1), 0, Inf)),
        HierarchicalNormal(0.5, truncated(Normal(0, 1), 0, Inf)),
        Intercept(; intercept = Normal()),
        MA(truncated(Normal(0.0, 0.05), -1, 1)),
        DiffLatentModel(RandomWalk(), Normal(); d = 2))
        @test m isa AbstractLatentModel
    end
end

@testitem "bare Distribution in a PATH slot auto-wraps as a constant Intercept" begin
    using ComposableTuringIDModels, Distributions, Random
    Random.seed!(210)
    # A length-`n` PATH slot (an innovation `ϵ_t`, or a latent `Z` / `rt`) given a
    # bare `Distribution` is wrapped in an `Intercept` at construction, so the
    # stored field is a constant path — never a scalar that would silently produce
    # a length-1 result.
    @test RandomWalk(; ϵ_t = Normal()).ϵ_t isa Intercept
    @test AR(; ϵ_t = Normal()).ϵ_t isa Intercept
    @test MA(; ϵ_t = Normal()).ϵ_t isa Intercept
    @test DirectInfections(; Z = Normal()).Z isa Intercept
    @test ExpGrowthRate(; rt = Normal()).rt isa Intercept
    @test Renewal(; generation_time = [0.2, 0.3, 0.5], rt = Normal()).rt isa
          Intercept
    # And the resulting model yields a length-`n` path, not a scalar.
    @test length(as_turing_model(RandomWalk(; ϵ_t = Normal()), 6)()) == 6
    @test length(as_turing_model(AR(; ϵ_t = Normal()), 6)()) == 6
    @test length(as_turing_model(MA(; ϵ_t = Normal()), 6)()) == 6
    @test length(as_turing_model(DirectInfections(; Z = Normal()), 6)().I_t) == 6
    @test length(as_turing_model(
        Renewal(; generation_time = [0.2, 0.3, 0.5], rt = Normal()), 8)().I_t) ==
          8
    # A process, an explicit `IID`, or an explicit `Intercept` passes through
    # unchanged (no double-wrapping).
    @test RandomWalk(; ϵ_t = RandomWalk()).ϵ_t isa RandomWalk
    @test RandomWalk(; ϵ_t = IID(Normal())).ϵ_t isa IID
    @test RandomWalk(; ϵ_t = Intercept(Normal())).ϵ_t isa Intercept
    # A per-step PARAMETER slot is UNAFFECTED: a bare `Distribution` stays a scalar
    # constant (drawn with a native tilde), not an `Intercept`.
    @test AR(; damp = Normal()).damp isa Distribution
    @test !(AR(; damp = Normal()).damp isa Intercept)
    @test NormalError(truncated(Normal(0, 1), 0, Inf)).std isa Distribution
    # Manipulators that thread an inner LATENT path wrap a bare `Distribution` too.
    @test DiffLatentModel(; model = Normal(), init = [Normal()]).model isa Intercept
    @test TransformLatentModel(Normal(), x -> x).model isa Intercept
    @test BroadcastLatentModel(
        Normal(); period = 7, broadcast_rule = RepeatEach()).model isa Intercept
    @test RecordExpectedLatent(Normal()).model isa Intercept
    @test Hierarchy(; across = Normal()).across isa Intercept
    @test Hierarchy(; across = IID(Normal())).across isa IID   # process passes through
    # Combine/Concat members are namespaced (prefix-wrapped) around the wrapped
    # constant, and still generate a length-`n` series with a bare member.
    @test length(as_turing_model(ConcatLatentModels([Normal(), AR()]), 10)()) == 10
    @test length(as_turing_model(CombineLatentModels([Normal(), AR()]), 10)()) == 10
end

@testitem "priors compose as submodels (distribution and latent)" begin
    using ComposableTuringIDModels, Distributions, Turing, Random
    Random.seed!(103)
    # A RandomWalk used directly as a (time-varying) prior returns a length-n path.
    @test length(as_turing_model(RandomWalk(), 8)()) == 8
    # A scalar parameter is a plain native-tilde draw (no submodel): NormalError's
    # σ has no `.θ` namespace path.
    draw = rand(as_turing_model(NormalError(), missing, fill(10.0, 5)))
    @test any(k -> occursin("σ", string(k)), keys(draw))
end
