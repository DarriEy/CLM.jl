#!/usr/bin/env julia
# ==========================================================================
# Conservation checks that must run with bounds checking at its DEFAULT.
#
# WHY THIS FILE EXISTS, AND WHY IT IS NOT PART OF runtests.jl
# ----------------------------------------------------------
# `Pkg.test` (and therefore julia-actions/julia-runtest, and therefore CI)
# forces `--check-bounds=yes`. That is not a neutral choice here: disabling
# `@inbounds` also blocks the LLVM vectorization that miscompiles a
# loop-carried read-after-write inside a `@kernel` body. So the entire test
# suite is STRUCTURALLY BLIND to that bug class — a kernel can be silently
# non-conservative in production while CI stays green.
#
# Three confirmed instances, all invisible to `Pkg.test`:
#   * the Richards flux loop (fixed in 518648f) — soil absorbed ~72% of
#     infiltration;
#   * `_lvt_column_kernel!` (litter vertical transport) — soil C/N transport
#     off by ~6 mass units, up to 400% relative at nlevdecomp = 5;
#   * `_ch4tran_column_kernel!` — CH4 −7.2%, O2 −5.1%.
#
# The mechanism (see the comment in `_tridiag_multi_kernel!`): KA's CPU
# backend emits `@simd ivdep` over workitems, and the resulting
# parallel-loop metadata reaches NESTED loops in the kernel body too. LLVM
# then vectorizes an inner level loop and ignores its loop-carried
# dependency. It is NOT a thread race — it reproduces with a single column.
#
# So: this file is run by its own CI job WITHOUT `--check-bounds=yes`. Adding
# these tests to runtests.jl would defeat their purpose.
#
# Run locally with:
#     julia --project=. test/conservation/run_default_bounds.jl
# and confirm the same numbers appear under `--check-bounds=yes` — agreement
# between the two modes is the actual invariant being tested.
# ==========================================================================

using Test
using CLM

include(joinpath(@__DIR__, "helpers.jl"))

@assert Base.JLOptions().check_bounds != 1 """
run_default_bounds.jl was started with --check-bounds=yes, which defeats its
purpose: that flag disables @inbounds and thereby suppresses the very
miscompilation these tests exist to catch. Run it with bounds checking left
at its default.
"""

@testset "conservation under default bounds checking" begin

    # ---------------------------------------------------------------------
    # Litter vertical transport: the Patankar discretisation with zero-flux
    # boundary rows conserves mass exactly, so
    #     Σⱼ (conc_after − conc_before)·dzⱼ == Σⱼ sourceⱼ·dzⱼ
    # Deep columns matter — the miscompile needs enough levels to vectorize.
    # ---------------------------------------------------------------------
    @testset "litter_vert_transp! conserves tracer mass" begin
        for (nc, nlevdecomp) in ((12, 20), (12, 10), (12, 5), (1, 20), (12, 30))
            r = _lvt_conservation_case(nc, nlevdecomp)
            @test r.max_abs_err < 1e-9
            @test r.max_rel_err < 1e-9
        end
    end

    # ---------------------------------------------------------------------
    # SNICAR: the host `snicar_rt!` is a plain per-column loop, not a kernel,
    # so it is an independent reference for the device kernel. They agreed
    # BIT-EXACTLY under --check-bounds=yes even while the device kernel was
    # broken, which is what proves the transliteration faithful and the
    # optimizer culpable. Here they must agree with bounds checking at its
    # default too.
    # ---------------------------------------------------------------------
    @testset "snicar device kernel matches the host reference" begin
        for flg_slr_in in (1, 2)
            d = _snicar_host_vs_device(flg_slr_in)
            @test d.max_albedo_diff < 1e-12
            @test d.max_flx_abs_diff < 1e-12
        end
    end

    # ---------------------------------------------------------------------
    # CH4/O2 transport. Its own conservation diagnostic must close, and — the
    # sharper test — the answer must not depend on the bounds-checking mode.
    # A reference produced under `--check-bounds=yes` is committed alongside
    # this file; agreement with it IS the invariant, because a miscompiled
    # reduction shows up as a difference between the two modes and nothing
    # else. `ch4_ebul_total` was accumulated straight into an array element
    # inside a level loop and then fed into the CH4 source at the water table,
    # so the corruption landed in the tracer solve: the error peaked exactly
    # at jwt, decayed upward, vanished below, and left O2 untouched.
    # ---------------------------------------------------------------------
    @testset "ch4_tran! conserves and is bounds-mode invariant" begin
        r = _ch4_tran_case()
        @test r.max_errch4 < 1e-9
        @test all(isfinite, r.conc_ch4) && all(isfinite, r.conc_o2)
        # values recorded under --check-bounds=yes on the fixed code
        @test sum(r.conc_ch4) ≈ 0.021826987644967245 rtol = 1e-13
        @test sum(r.conc_o2)  ≈ 1.7981397559679242   rtol = 1e-13
    end
end
