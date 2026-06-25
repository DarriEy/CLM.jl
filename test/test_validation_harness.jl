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
        end
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
