# Tests for the optional Catalyst extension (`EpiAwarePrototypeCatalystExt`).
# `Catalyst` + `ModelingToolkit` are test-only deps (see test/Project.toml);
# loading them triggers the extension, which supplies the `ReactionSystem`
# constructor / sampling for the public, exported `CatalystODEParams` type.
# `CatalystODEParams` is model-agnostic, so the tests exercise it on more than one
# reaction network (SIR and SEIR), and sampling / solution indexing are symbolic
# (no positional-index bookkeeping).

@testitem "Catalyst extension loads and CatalystODEParams is public" begin
    using EpiAwarePrototype, Catalyst, ModelingToolkit
    ext = Base.get_extension(EpiAwarePrototype, :EpiAwarePrototypeCatalystExt)
    @test ext !== nothing
    # The type is a first-class exported public component (defined in `src/`),
    # not something dug out of the extension module.
    @test CatalystODEParams isa Type
    @test isdefined(EpiAwarePrototype, :CatalystODEParams)
end

@testitem "CatalystODEParams errors helpfully before Catalyst is loaded" begin
    # The fallback constructor lives in `src/`; without a `ReactionSystem` it
    # raises an informative error rather than a bare MethodError.
    using EpiAwarePrototype
    @test_throws ArgumentError CatalystODEParams(:not_a_reaction_system;
        tspan = (0.0, 1.0), u0_priors = [], p_priors = [])
end

@testitem "CatalystODEParams samples (u0, p) for an arbitrary (SIR) network" begin
    using EpiAwarePrototype, Catalyst, ModelingToolkit, OrdinaryDiffEq,
          Distributions, Random
    Random.seed!(101)

    sir = @reaction_network begin
        β, S + I --> 2I
        γ, I --> R
    end
    params = CatalystODEParams(sir;
        tspan = (0.0, 30.0),
        u0_priors = [sir.S => Beta(99, 1), sir.I => Beta(1, 99), sir.R => 0.0],
        p_priors = [sir.β => LogNormal(log(0.3), 0.05),
            sir.γ => LogNormal(log(0.1), 0.05)])

    # Sampling returns symbolic `symbol => value` maps, one per species / rate.
    u0, p = as_turing_model(params, nothing)()
    @test length(u0) == 3
    @test length(p) == 2
    @test all(x -> x isa Pair, u0)
    @test all(x -> x isa Pair, p)
    # The fixed (Real) R initial condition rides through as a constant value
    # (compare on the values, not the symbolic keys, to avoid symbolic `==`).
    @test 0.0 in [pr.second for pr in u0]

    # Distribution-valued specs are sampled with flat, symbol-named keys; the
    # fixed (Real) R initial condition is NOT sampled, so it is absent.
    draw = rand(as_turing_model(params, nothing))
    nms = string.(collect(keys(draw)))
    @test all(n -> n in nms, ["β", "γ", "S", "I"])
    @test !("R" in nms)
end

@testitem "CatalystODEParams enforces a prior for every species and parameter" begin
    using EpiAwarePrototype, Catalyst, ModelingToolkit, OrdinaryDiffEq, Distributions
    sir = @reaction_network begin
        β, S + I --> 2I
        γ, I --> R
    end
    # Missing the R species spec.
    @test_throws ArgumentError CatalystODEParams(sir;
        tspan = (0.0, 30.0),
        u0_priors = [sir.S => Beta(99, 1), sir.I => Beta(1, 99)],
        p_priors = [sir.β => LogNormal(log(0.3), 0.05),
            sir.γ => LogNormal(log(0.1), 0.05)])
end

@testitem "Catalyst SEIR trajectory matches the hand-coded SEIR" begin
    using EpiAwarePrototype, Catalyst, ModelingToolkit, OrdinaryDiffEq, Distributions

    tspan = (0.0, 60.0)
    seir = @reaction_network begin
        β, S + I --> E + I
        α, E --> I
        γ, I --> R
    end
    # Dirac priors so both models are evaluated at identical states/rates: the
    # only difference is hand-coded vs Catalyst-generated dynamics.
    fixed = (β = 0.31, α = 0.095, γ = 0.105, initial_infs = 0.02)
    I0 = fixed.initial_infs * fixed.α / (fixed.α + fixed.γ)
    E0 = fixed.initial_infs * fixed.γ / (fixed.α + fixed.γ)

    # Fixed (Real) specs so the Catalyst trajectory is deterministic and matches
    # the hand-coded SEIR evaluated at the same states/rates.
    catalyst = CatalystODEParams(seir;
        tspan,
        u0_priors = [seir.S => 1.0 - fixed.initial_infs, seir.E => E0,
            seir.I => I0, seir.R => 0.0],
        p_priors = [seir.β => fixed.β, seir.α => fixed.α, seir.γ => fixed.γ])
    handcoded = SEIRParams(; tspan,
        infectiousness = Dirac(fixed.β), incubation_rate = Dirac(fixed.α),
        recovery_rate = Dirac(fixed.γ), initial_prop_infected = Dirac(fixed.initial_infs))

    # Symbolic solution indexing: pull the infectious compartment by its handle,
    # no stored-index lookup.
    cat_proc = ODEProcess(params = catalyst, sol2infs = sol -> sol[seir.I, :],
        solver_options = Dict(:saveat => 1.0))
    hand_proc = ODEProcess(params = handcoded, sol2infs = sol -> sol[3, :],
        solver_options = Dict(:saveat => 1.0))

    cat_I = as_turing_model(cat_proc, nothing)().I_t
    hand_I = as_turing_model(hand_proc, nothing)().I_t

    @test length(cat_I) == length(hand_I)
    # #46 reported agreement to ≈3.5e-9 at matched priors/seed.
    @test maximum(abs.(cat_I .- hand_I)) < 1e-6
    @test isapprox(cat_I, hand_I; atol = 1e-7)
end

@testitem "CatalystODEParams composes into an ODEProcess and exposes no latent" begin
    using EpiAwarePrototype, Catalyst, ModelingToolkit, OrdinaryDiffEq,
          Distributions, LogExpFunctions, Random
    Random.seed!(102)

    seir = @reaction_network begin
        β, S + I --> E + I
        α, E --> I
        γ, I --> R
    end
    params = CatalystODEParams(seir;
        tspan = (0.0, 50.0),
        u0_priors = [seir.S => 0.99, seir.E => 0.005,
            seir.I => 0.005, seir.R => 0.0],
        p_priors = [seir.β => LogNormal(log(0.3), 0.05),
            seir.α => LogNormal(log(0.1), 0.05),
            seir.γ => LogNormal(log(0.1), 0.05)])
    N = 1000.0
    proc = ODEProcess(params = params,
        sol2infs = sol -> softplus.(N .* sol[seir.I, :]),
        solver_options = Dict(:saveat => 1.0))

    out = as_turing_model(proc, nothing)()
    @test length(out.I_t) == 51
    @test all(>=(0), out.I_t)
    @test isnothing(out.Z_t)

    draw = rand(as_turing_model(proc, nothing))
    nms = string.(collect(keys(draw)))
    @test all(n -> n in nms, ["β", "α", "γ"])
end

@testitem "Catalyst SEIR + observation samples under ForwardDiff NUTS" tags=[:forwarddiff] begin
    using EpiAwarePrototype, Catalyst, ModelingToolkit, OrdinaryDiffEq,
          Distributions, LogExpFunctions, Turing, ADTypes, Random
    using DynamicPPL: @varname
    Random.seed!(103)

    N = 763
    n_days = 14
    seir = @reaction_network begin
        β, S + I --> E + I
        α, E --> I
        γ, I --> R
    end
    params = CatalystODEParams(seir;
        tspan = (0.0, Float64(n_days)),
        u0_priors = [seir.S => 0.99, seir.E => Beta(2, 200),
            seir.I => Beta(2, 200), seir.R => 0.0],
        p_priors = [seir.β => LogNormal(-0.5, 0.4),
            seir.α => Gamma(8, 0.05), seir.γ => Gamma(8, 0.03125)])
    obs = TransformObservationModel(PoissonError(), x -> softplus.(N .* x))
    process = ODEProcess(params = params,
        sol2infs = sol -> sol[seir.I, :],
        solver_options = Dict(:saveat => 1.0))
    model = EpiAwareModel(process, obs)

    sim = as_turing_model(model, fill(missing, n_days + 1), n_days + 1)()
    y_obs = sim.generated_y_t
    @test length(y_obs) == n_days + 1

    # ForwardDiff is the supported AD path for ODE infection models (#46);
    # Mooncake-driven NUTS through the solver is a separate, pre-existing gap.
    # `sample` returns a FlexiChains chain, indexed by variable name directly.
    # This exercises symbolic-map remake carrying ForwardDiff `Dual`s.
    chain = sample(as_turing_model(model, y_obs, n_days + 1),
        NUTS(; adtype = AutoForwardDiff()), 20; progress = false)
    βs = vec(chain[@varname(β)])
    @test length(βs) == 20
    @test all(isfinite, βs)
end
