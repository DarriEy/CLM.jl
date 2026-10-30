using Test
using CLM

# ============================================================================
# @Const contract audit (src/infrastructure/kconst_check.jl).
#
# KernelAbstractions' @Const(x) promises x is read-only and aliases nothing the
# kernel writes; on the CPU it was lowered through Base.Experimental.@aliasscope,
# which Julia >= 1.11 miscompiles (JuliaLang/julia#63129, KA #652; it hit CNFUN's
# _fun_p2_kernel!, PR #365). KA >= 0.9.45 no longer emits the aliasscope on
# Julia >= 1.11 (Project.toml compat), but the contract still matters on 1.10 and
# for GPU backends. This test:
#   1. statically: no @kernel in src/ writes a @Const argument (or a view of one);
#   2. checks the static and the opt-in dynamic checker catch known hazards;
#      (incl. that disjoint @view slices of one packed array -- CNFUN's layout -- pass).
# ============================================================================

const _KCT_FILE = joinpath(@__DIR__, "kconst_toy_kernels.jl")
include(_KCT_FILE)

@testset "@Const contract (kconst_check)" begin
    @testset "static: no @Const argument written in src/" begin
        reg = CLM.kconst_kernel_signatures()
        @test length(reg) > 500                            # the parser sees the kernels
        @test count(v -> any(v[2]), values(reg)) > 500
        @test reg[:_fun_p2_kernel!][1][1] === :mask_soilp
        w = CLM.kconst_static_writes()
        isempty(w) || foreach(x -> @info("static @Const write", x...), w)
        @test isempty(w)
    end

    @testset "static checker catches the toy hazards" begin
        w = CLM.kconst_static_writes(_KCT_FILE)
        ks = Set((x.kernel, x.const_arg, x.via) for x in w)
        @test (:_kct_write_const!, :a, :a) in ks
        @test (:_kct_write_view!, :A, :v) in ks
        @test !any(x -> x.kernel === :_kct_copy!, w)
    end

    @testset "dynamic checker (opt-in _launch! mode)" begin
        @test CLM.kconst_check_mode() === :off
        CLM._KCONST_REGISTRY[] = merge(CLM.kconst_kernel_signatures(),
                                       CLM.kconst_kernel_signatures(_KCT_FILE))
        try
            CLM.set_kconst_check!(:warn)
            CLM.clear_kconst_findings!()
            f() = CLM.KCONST_FINDINGS
            kinds() = Set((x.kernel, x.kind, x.const_arg, x.other_arg) for x in f())

            # clean launch: nothing reported
            a = rand(8); out = zeros(8)
            Base.invokelatest(CLM._launch!, _kct_copy!, out, a)
            @test isempty(f()) && out == 2 .* a

            # same array as output and @Const input
            x = rand(8)
            @test_logs (:warn,) match_mode = :any Base.invokelatest(CLM._launch!, _kct_copy!, x, x)
            @test (:_kct_copy!, :alias_written, :a, :out) in kinds()
            @test (:_kct_copy!, :const_written, :a, :_) in kinds()

            # overlapping views of one parent (shifted by one element)
            CLM.clear_kconst_findings!()
            P = rand(9)
            Base.invokelatest(CLM._launch!, _kct_copy!, view(P, 2:9), view(P, 1:8))
            @test (:_kct_copy!, :alias_written, :a, :out) in kinds()

            # DISJOINT column views of one packed array (the CNFUN pattern): clean
            CLM.clear_kconst_findings!()
            Q = rand(8, 3)
            Base.invokelatest(CLM._launch!, _kct_copy!, view(Q, :, 1), view(Q, :, 2))
            @test isempty(f())

            # A kernel writing its @Const argument DIRECTLY cannot run on the CPU
            # (Base.Experimental.Const has no setindex!, nor do views of it), so that
            # case is covered statically above; dynamically it shows up only through
            # an alias, as :const_written.
            @test_throws Exception Base.invokelatest(CLM._launch!, _kct_write_const!, zeros(8), rand(8))

            # :error mode throws
            CLM.set_kconst_check!(:error)
            @test_throws ErrorException Base.invokelatest(CLM._launch!, _kct_copy!, x, x)
        finally
            CLM.set_kconst_check!(:off)
            CLM._KCONST_REGISTRY[] = nothing
            CLM.clear_kconst_findings!()
        end
        @test Base.invokelatest(CLM.kconst_check_mode) === :off
    end
end
