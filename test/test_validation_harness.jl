# CI-safe tests for the validation-harness matrix schema (test/validation/matrix.jl).
# The runner (scripts/validation/validate.jl) is exercised LOCALLY against machine-
# local domain data; here we validate the dependency-free declarative layer that
# CI can check without any external inputs.

@testset "validation harness — matrix schema" begin
    include(joinpath(@__DIR__, "validation", "matrix.jl"))

    M = validation_matrix()

    @testset "matrix well-formed" begin
        @test !isempty(M)
        ids = [c.id for c in M]
        @test length(ids) == length(unique(ids))          # unique slugs
        for c in M
            @test c.mode in (:sp, :cn, :fates_sp, :fates_bgc)
            @test c.backend in (:cpu, :metal, :cuda, :amdgpu)
            @test c.init in (:cold, :warm)
            @test haskey(DEPTH_DEFAULT_STEPS, c.depth)
            @test !isempty(c.oracles)
            @test all(o -> o in ORACLE_KINDS, c.oracles)
            @test c.flags isa NamedTuple
            @test c.tier in TIERS
        end
        @test any(c -> c.tier === :pr, M)        # a fast per-PR lane exists
        @test any(c -> :parity in c.oracles, M)  # T1 anchor present
    end

    @testset "pairwise covering array is complete + terminates" begin
        for flags in ([:a, :b, :c, :d, :e], [:f, :g, :h, :i, :j, :k])
            rows = pairwise_binary(flags)
            n = length(flags)
            bool = [[f in keys(r) for f in flags] for r in rows]
            push!(bool, fill(false, n))           # the all-off baseline row exists too
            covered = Set()
            for br in bool, i in 1:n-1, j in i+1:n
                push!(covered, (i, j, br[i], br[j]))
            end
            @test length(covered) == binomial(n, 2) * 4   # every 2-way value-pair present
            @test length(rows) <= 4 * binomial(n, 2)      # bounded (no runaway)
        end
        @test isempty(pairwise_binary([:only_one]))       # <2 flags ⇒ no pairs
    end

    @testset "tier filtering is monotone (pr ⊆ nightly ⊆ weekly)" begin
        M = validation_matrix()
        rank = Dict(:pr => 1, :nightly => 2, :weekly => 3)
        sel(t) = filter(c -> rank[c.tier] <= rank[t], M)
        @test Set(c.id for c in sel(:pr)) ⊆ Set(c.id for c in sel(:nightly))
        @test Set(c.id for c in sel(:nightly)) ⊆ Set(c.id for c in sel(:weekly))
        @test length(sel(:weekly)) == length(M)
    end

    @testset "vcfg validates its fields" begin
        @test_throws ErrorException vcfg("bad-mode"; mode=:nope)
        @test_throws ErrorException vcfg("bad-oracle"; oracles=[:not_a_real_oracle])
        @test_throws ErrorException vcfg("bad-depth"; depth=:eon)
        @test_throws ErrorException vcfg("bad-backend"; backend=:tpu)
        @test_throws ErrorException vcfg("bad-init"; init=:lukewarm)
        # a well-formed entry constructs and carries defaults
        c = vcfg("ok"; mode=:cn, flags=(use_luna=true,), depth=:day)
        @test c.id == "ok"
        @test c.mode === :cn
        @test c.flags.use_luna === true
        @test c.domain === :bow
        @test c.oracles == [:conservation]
    end

    @testset "depth step counts monotone" begin
        @test DEPTH_DEFAULT_STEPS[:smoke] < DEPTH_DEFAULT_STEPS[:day] <
              DEPTH_DEFAULT_STEPS[:season] < DEPTH_DEFAULT_STEPS[:year] <
              DEPTH_DEFAULT_STEPS[:multiyear]
    end
end
