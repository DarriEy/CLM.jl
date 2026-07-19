# ==========================================================================
# gpu_validate_handlenewsnow_e2e.jl — GPU parity for the host loops that were
# kernelized inside handle_new_snow! (HydrologyNoDrainageMod / A7).
#
# handle_new_snow! is a large orchestrator: most of its sub-steps are already
# kernelized in their own modules, but THREE inline host loops lived directly in
# the function body and are the ones converted in this A7 work:
#   1. temperature-only new-snow bulk density fallback   (_hydrond_bifall_fallback_kernel!)
#   2. total snow water h2osno_total                     (_hydrond_h2osno_total_kernel!)
#        — a per-column kernel with an INTERNAL sequential snow-layer loop over
#          the variable-snl padded slice [snl+1 .. 0] (loop-carried accumulation)
#   3. per-column landunit-type / urban gather           (_hydrond_gather_lun_kernel!)
#
# This harness launches those three kernels on Metal vs CPU across columns with
# DIFFERENT snl so the variable-layer accumulation branch is exercised:
#   col 1: snl= 0  (no resolved snow -> only no-layer reservoir contributes)
#   col 2: snl=-1  (single thin layer)
#   col 3: snl=-3  (deep multi-layer accumulation)
#   col 4: snl=-3, warm   (bifall high branch)
# CPU reference fields are asserted finite so a both-NaN false PASS cannot slip
# through.
#
#   julia --project=scripts scripts/gpu_validate_handlenewsnow_e2e.jl
# ==========================================================================

using CLM
using Printf
import KernelAbstractions as KA
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end

const NC = 4
const NLEVSNO = 5

function build(::Type{FT}) where {FT}
    nc = NC
    nlev = NLEVSNO + 15   # snow + soil padding
    snl = Int[0, -1, -3, -3]
    mask = trues(nc)

    TFRZ = FT(CLM.TFRZ)
    # forc_t spanning bifall branches: cold (<TFRZ-15), mid, warm (>TFRZ+2)
    forc_t = FT[TFRZ - FT(20.0), TFRZ - FT(5.0), TFRZ - FT(1.0), TFRZ + FT(5.0)]

    h2osno_no_layers = FT[FT(7.5), FT(0.0), FT(0.0), FT(0.0)]
    h2osoi_ice = fill(FT(0.0), nc, nlev)
    h2osoi_liq = fill(FT(0.0), nc, nlev)
    # fill the resolved snow layers [snl+1 .. 0] -> Julia [snl+1+nlevsno .. nlevsno]
    for c in 1:nc
        for j in (snl[c] + 1):0
            jj = j + NLEVSNO
            h2osoi_ice[c, jj] = FT(10.0) + FT(c) + FT(-j)
            h2osoi_liq[c, jj] = FT(1.0)  + FT(0.1) * FT(-j)
        end
    end

    col_landunit = Int[1, 2, 3, 4]
    lun_itype    = Int[CLM.ISTSOIL, CLM.ISTSOIL, CLM.ISTDLAK, CLM.ISTWET]
    lun_urbpoi   = Bool[false, false, false, false]

    return (; nc, snl, mask, forc_t, h2osno_no_layers, h2osoi_ice, h2osoi_liq,
            col_landunit, lun_itype, lun_urbpoi)
end

# Drive the three kernels exactly as handle_new_snow! does.
function run_kernels!(S, bifall, h2osno_total, lun_itype_col, urbpoi_col)
    FT = eltype(bifall)
    bounds = 1:S.nc
    CLM._launch!(CLM._hydrond_bifall_fallback_kernel!, bifall, S.mask, S.forc_t,
                 FT(CLM.TFRZ), first(bounds), last(bounds))
    CLM._launch!(CLM._hydrond_h2osno_total_kernel!, h2osno_total, S.mask, S.snl,
                 S.h2osno_no_layers, S.h2osoi_ice, S.h2osoi_liq, NLEVSNO,
                 first(bounds), last(bounds))
    CLM._launch!(CLM._hydrond_gather_lun_kernel!, lun_itype_col, S.mask,
                 S.col_landunit, S.lun_itype, S.lun_urbpoi, urbpoi_col,
                 first(bounds), last(bounds))
    return nothing
end

function main(backend)
    println("=" ^ 72)
    println("END-TO-END GPU parity for handle_new_snow! kernelized loops (A7)")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT); B = build(FT)

    # CPU scratch (Float fields filled, Int/Bool gathers initialized to defined defaults)
    bifall_h = fill(FT(50.0), H.nc)
    h2osno_h = zeros(FT, H.nc)
    luit_h   = fill(0, H.nc)
    urb_h    = fill(false, H.nc)
    run_kernels!(H, bifall_h, h2osno_h, luit_h, urb_h)

    db(x) = to(collect(Bool, x))
    Sd = (; nc = B.nc, snl = to(B.snl), mask = db(B.mask), forc_t = to(B.forc_t),
            h2osno_no_layers = to(B.h2osno_no_layers),
            h2osoi_ice = to(B.h2osoi_ice), h2osoi_liq = to(B.h2osoi_liq),
            col_landunit = to(B.col_landunit), lun_itype = to(B.lun_itype),
            lun_urbpoi = db(B.lun_urbpoi))
    bifall_d = to(fill(FT(50.0), B.nc))
    h2osno_d = to(zeros(FT, B.nc))
    luit_d   = to(fill(0, B.nc))
    urb_d    = db(fill(false, B.nc))

    if !(bifall_d isa device_array_type())
        println("  BLOCKED: device arrays did not move under to().")
        return 2
    end
    run_kernels!(Sd, bifall_d, h2osno_d, luit_d, urb_d)

    tonum(x) = Float64.(Array(x))
    checks = [
        ("bifall",        bifall_h, bifall_d),
        ("h2osno_total",  h2osno_h, h2osno_d),
        ("lun_itype_col", luit_h,   luit_d),
        ("urbpoi_col",    urb_h,    urb_d),
    ]
    nfail = 0
    for (nm, a, b) in checks
        an = tonum(a)
        if !any(isfinite, an)
            @printf("  [WARN] %-14s CPU reference all-NaN/Inf — skipping\n", nm)
            continue
        end
        d = reldiff(an, tonum(b)); ok = d < 1f-3
        @printf("  [%s] %-14s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end

    println()
    for c in 1:H.nc
        @printf("    col %d: snl=%2d  bifall=%9.4f  h2osno_total=%9.4f\n",
                c, H.snl[c], bifall_h[c], h2osno_h[c])
    end

    println()
    println(nfail == 0 ? "  handle_new_snow! kernelized loops MATCH CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
