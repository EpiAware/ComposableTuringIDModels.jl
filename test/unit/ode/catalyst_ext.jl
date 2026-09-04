# Tests for the optional Catalyst extension (`ComposableTuringIDModelsCatalystExt`).
# `Catalyst` + `ModelingToolkit` are test-only deps (see test/Project.toml);
# loading them triggers the extension, which supplies the `ReactionSystem`
# constructor / sampling for the public, exported `CatalystODEParams` type.
# `CatalystODEParams` is model-agnostic, so the tests exercise it on more than one
# reaction network (SIR and SEIR), and sampling / solution indexing are symbolic
# (no positional-index bookkeeping).

@testitem "Catalyst extension loads and CatalystODEParams is public" begin
    using ComposableTuringIDModels, Catalyst, ModelingToolkit
    ext = Base.get_extension(ComposableTuringIDModels, :ComposableTuringIDModelsCatalystExt)
    @test ext !== nothing
    # The type is a first-class exported public component (defined in `src/`),
    # not something dug out of the extension module.
    @test CatalystODEParams isa Type
    @test isdefined(ComposableTuringIDModels, :CatalystODEParams)
end

@testitem "CatalystODEParams errors helpfully before Catalyst is loaded" begin
    # The fallback constructor lives in `src/`; without a `ReactionSystem` it
    # raises an informative error rather than a bare MethodError.
    using ComposableTuringIDModels
    @test_throws ArgumentError CatalystODEParams(
        :not_a_reaction_system;
        tspan = (0.0, 1.0), u0_priors = [], p_priors = []
    )
end

@testitem "CatalystODEParams samples (u0, p) for an arbitrary (SIR) network" begin
    using ComposableTuringIDModels, Catalyst, ModelingToolkit, OrdinaryDiffEq,
        Distributions, Random
    Random.seed!(101)

    sir = @reaction_network begin
        β, S + I --> 2I
        γ, I --> R
    end
    params = CatalystODEParams(
        sir;
        tspan = (0.0, 30.0),
        u0_priors = [sir.S => Beta(99, 1), sir.I => Beta(1, 99), sir.R => 0.0],
        p_priors = [
            sir.β => LogNormal(log(0.3), 0.05),
            sir.γ => LogNormal(log(0.1), 0.05),
        ]
    )

    # Sampling returns plain vectors, one entry per species / rate, ordered the
    # way the compiled problem stores them rather than the way the priors were
    # written.
    u0, p = as_turing_model(params, nothing)()
    @test length(u0) == 3
    @test length(p) == 2
    @test u0 isa Vector{<:Real}
    @test p isa Vector{<:Real}
    # The fixed (Real) R initial condition rides through as a constant value.
    @test 0.0 in u0
    # Catalyst sorts species internally, so the stored layout is what says where
    # each one landed; the sampled vector has to agree with it.
    for (name, slot) in pairs(params.layout.slots)
        spec = only(sp for sp in params.u0_specs if sp.name == name)
        spec.spec isa Real && @test u0[slot] == spec.spec
    end

    # Distribution-valued specs are sampled with flat, symbol-named keys; the
    # fixed (Real) R initial condition is NOT sampled, so it is absent.
    draw = rand(as_turing_model(params, nothing))
    nms = string.(collect(keys(draw)))
    @test all(n -> n in nms, ["β", "γ", "S", "I"])
    @test !("R" in nms)
end

@testitem "CatalystODEParams builds a FullSpecialize problem" begin
    using ComposableTuringIDModels, Catalyst, ModelingToolkit, OrdinaryDiffEq,
        Distributions, SciMLBase
    # An MTK-built problem defaults to `AutoDespecialize`, which wraps the
    # generated Jacobian in a `Float64`-typed `FunctionWrapper`. Reverse mode
    # runs the solve in `Float64` and differentiates it with nested
    # `ForwardDiff`, so the moment a stiff proposal switches
    # `AutoVern7(Rodas5P())` onto `Rodas5P` the Jacobian is called with a `Dual`
    # matrix no wrapper covers. A single gradient evaluation stays on the
    # non-stiff half and never reaches it, so only pinning the specialisation
    # guards this.
    sir = @reaction_network begin
        β, S + I --> 2I
        γ, I --> R
    end
    params = CatalystODEParams(
        sir;
        tspan = (0.0, 30.0),
        u0_priors = [sir.S => Beta(99, 1), sir.I => Beta(1, 99), sir.R => 0.0],
        p_priors = [
            sir.β => LogNormal(log(0.3), 0.05),
            sir.γ => LogNormal(log(0.1), 0.05),
        ]
    )
    @test SciMLBase.specialization(params.prob.f) === SciMLBase.FullSpecialize
end

@testitem "CatalystODEParams enforces a prior for every species and parameter" begin
    using ComposableTuringIDModels, Catalyst, ModelingToolkit, OrdinaryDiffEq, Distributions
    sir = @reaction_network begin
        β, S + I --> 2I
        γ, I --> R
    end
    # Missing the R species spec.
    @test_throws ArgumentError CatalystODEParams(
        sir;
        tspan = (0.0, 30.0),
        u0_priors = [sir.S => Beta(99, 1), sir.I => Beta(1, 99)],
        p_priors = [
            sir.β => LogNormal(log(0.3), 0.05),
            sir.γ => LogNormal(log(0.1), 0.05),
        ]
    )
end

@testitem "Catalyst SEIR trajectory matches a hand-written SEIR" begin
    using ComposableTuringIDModels, Catalyst, ModelingToolkit, OrdinaryDiffEq, Distributions

    # A plain SEIR vector field, written out here so the Catalyst-generated
    # dynamics are checked against an independent implementation rather than
    # against nothing. It carries no Jacobian and the package exports no such
    # model: this is a reference solve, and it is the only hand-written
    # compartmental drift left anywhere in the repository outside SIR.
    function seir_vf(du, u, p, t)
        S, E, Iv, R = u
        β, α, γ = p
        du[1] = -β * S * Iv
        du[2] = β * S * Iv - α * E
        du[3] = α * E - γ * Iv
        du[4] = γ * Iv
        return nothing
    end

    tspan = (0.0, 60.0)
    # Fixed (Real) specs so the Catalyst trajectory is deterministic and the two
    # paths are evaluated at identical states and rates.
    fixed = (β = 0.31, α = 0.095, γ = 0.105, initial_infs = 0.02)
    E0 = fixed.initial_infs * fixed.γ / (fixed.α + fixed.γ)
    I0 = fixed.initial_infs * fixed.α / (fixed.α + fixed.γ)
    S0 = 1.0 - fixed.initial_infs

    seir = @reaction_network begin
        β, S + I --> E + I
        α, E --> I
        γ, I --> R
    end
    catalyst = CatalystODEParams(
        seir;
        tspan,
        u0_priors = [seir.S => S0, seir.E => E0, seir.I => I0, seir.R => 0.0],
        p_priors = [seir.β => fixed.β, seir.α => fixed.α, seir.γ => fixed.γ]
    )
    # Symbolic solution indexing: pull the infectious compartment by its handle,
    # no positional bookkeeping in user code.
    cat_proc = ODEProcess(
        params = catalyst, sol2infs = sol -> sol[seir.I, :],
        solver_options = Dict(:saveat => 1.0)
    )
    cat_I = as_turing_model(cat_proc, nothing)().I_t

    reference = solve(
        ODEProblem(seir_vf, [S0, E0, I0, 0.0], tspan, [fixed.β, fixed.α, fixed.γ]),
        AutoVern7(Rodas5P()); saveat = 1.0
    )
    ref_I = reference[3, :]

    @test length(cat_I) == length(ref_I)
    @test maximum(abs.(cat_I .- ref_I)) < 1.0e-6
    @test isapprox(cat_I, ref_I; atol = 1.0e-7)
end

@testitem "CatalystODEParams composes into an ODEProcess and exposes no latent" begin
    using ComposableTuringIDModels, Catalyst, ModelingToolkit, OrdinaryDiffEq,
        Distributions, LogExpFunctions, Random
    Random.seed!(102)

    seir = @reaction_network begin
        β, S + I --> E + I
        α, E --> I
        γ, I --> R
    end
    params = CatalystODEParams(
        seir;
        tspan = (0.0, 50.0),
        u0_priors = [
            seir.S => 0.99, seir.E => 0.005,
            seir.I => 0.005, seir.R => 0.0,
        ],
        p_priors = [
            seir.β => LogNormal(log(0.3), 0.05),
            seir.α => LogNormal(log(0.1), 0.05),
            seir.γ => LogNormal(log(0.1), 0.05),
        ]
    )
    N = 1000.0
    proc = ODEProcess(
        params = params,
        sol2infs = sol -> softplus.(N .* sol[seir.I, :]),
        solver_options = Dict(:saveat => 1.0)
    )

    out = as_turing_model(proc, nothing)()
    @test length(out.I_t) == 51
    @test all(>=(0), out.I_t)
    @test isnothing(out.Z_t)

    draw = rand(as_turing_model(proc, nothing))
    nms = string.(collect(keys(draw)))
    @test all(n -> n in nms, ["β", "α", "γ"])
end

@testitem "Catalyst SEIR + observation samples under ForwardDiff NUTS" tags = [:forwarddiff] begin
    using ComposableTuringIDModels, Catalyst, ModelingToolkit, OrdinaryDiffEq,
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
    params = CatalystODEParams(
        seir;
        tspan = (0.0, Float64(n_days)),
        u0_priors = [
            seir.S => 0.99, seir.E => Beta(2, 200),
            seir.I => Beta(2, 200), seir.R => 0.0,
        ],
        p_priors = [
            seir.β => LogNormal(-0.5, 0.4),
            seir.α => Gamma(8, 0.05), seir.γ => Gamma(8, 0.03125),
        ]
    )
    obs = TransformObservationModel(PoissonError(), x -> softplus.(N .* x))
    process = ODEProcess(
        params = params,
        sol2infs = sol -> sol[seir.I, :],
        solver_options = Dict(:saveat => 1.0)
    )
    model = IDModel(process, obs)

    sim = as_turing_model(model, fill(missing, n_days + 1), n_days + 1)()
    y_obs = sim.generated_y_t
    @test length(y_obs) == n_days + 1

    # ForwardDiff is the supported AD path for ODE infection models;
    # Mooncake-driven NUTS through the solver is a separate, pre-existing gap.
    # `sample` returns a FlexiChains chain, indexed by variable name directly.
    # This exercises symbolic-map remake carrying ForwardDiff `Dual`s.
    chain = sample(
        as_turing_model(model, y_obs, n_days + 1),
        NUTS(; adtype = AutoForwardDiff()), 20; progress = false
    )
    βs = vec(chain[@varname(β)])
    @test length(βs) == 20
    @test all(isfinite, βs)
end
