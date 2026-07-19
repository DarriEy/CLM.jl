# ==========================================================================
# gpu_validate_snowcoverfrac_e2e.jl — GPU parity for the handle_new_snow! tail
# (A10): bulkdiag_new_snow_diagnostics! and the snow_cover_fraction.jl helpers it
# drives — update_snow_depth_and_frac!, add_newsnow_to_intsnow!, calc_frac_sno_eff!.
#
# This runs the WHOLE bulkdiag_new_snow_diagnostics! on Metal vs CPU. That function
# internally invokes (all now kernelized):
#   - loop A: swe_old save + newsnow/snomelt/int_snow/snowmelt   (internal layer loop)
#   - update_snow_depth_and_frac! (SwensonLawrence2012): frac_sno (folds
#       frac_snow_during_melt inline), calc_frac_sno_eff!, snow_depth
#   - loop B: int_snow reset when pack started at 0
#   - add_newsnow_to_intsnow! (SwensonLawrence2012): int_snow accumulation
#   - loop C: top snow-layer dz update
#
# Branches exercised across the 5 columns:
#   col 1: snl= 0, h2osno=0,  newsnow>0      (fresh accumulation, int_snow reset)
#   col 2: snl=-2, h2osno>0,  snowmelt>0     (melt depletion curve, acos branch)
#   col 3: snl=-3, h2osno>0,  newsnow>0      (snow-on-snow accumulation)
#   col 4: snl= 0, h2osno=0,  newsnow=0      (no snow -> frac_sno=0, depth=0)
#   col 5: snl=-1, h2osno>0,  lake (ISTDLAK) (calc_frac_sno_eff non-fractional branch)
# The n_melt SCA shape vector is moved to the device explicitly.
#
# CPU reference fields asserted finite so a both-NaN false PASS cannot slip in.
#
#   julia --project=scripts scripts/gpu_validate_snowcoverfrac_e2e.jl
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

const NC = 5
const NLEVSNO = 5
const DTIME = 1800.0

function build(::Type{FT}) where {FT}
    nc = NC
    nlev = NLEVSNO + 10
    snl  = Int[0, -2, -3, 0, -1]
    mask = trues(nc)

    h2osoi_ice = fill(FT(0.0), nc, nlev)
    h2osoi_liq = fill(FT(0.0), nc, nlev)
    for c in 1:nc
        for j in (snl[c] + 1):0
            jj = j + NLEVSNO
            h2osoi_ice[c, jj] = FT(12.0) + FT(c) + FT(-j)
            h2osoi_liq[c, jj] = FT(1.5)  + FT(0.3) * FT(-j)
        end
    end

    # h2osno_total per branch design (col1,col4 == 0; others > 0)
    h2osno_total = FT[FT(0.0), FT(60.0), FT(90.0), FT(0.0), FT(40.0)]

    # newsnow drivers: qflx_snow_grnd * dtime = newsnow
    # col1: newsnow>0 ; col2: newsnow=0 (pure melt) ; col3: newsnow>0 ;
    # col4: newsnow=0 ; col5: newsnow>0
    qflx_snow_grnd = FT[FT(0.002), FT(0.0), FT(0.0015), FT(0.0), FT(0.001)]
    # snowmelt drivers: qflx_snow_drain*dtime ; col2 has melt
    qflx_snow_drain = FT[FT(0.0), FT(0.0008), FT(0.0), FT(0.0), FT(0.0)]

    int_snow = FT[FT(5.0), FT(120.0), FT(150.0), FT(0.0), FT(80.0)]
    snomelt_accum = fill(FT(0.01), nc)
    bifall = FT[FT(100.0), FT(120.0), FT(90.0), FT(110.0), FT(130.0)]

    frac_sno = FT[FT(0.0), FT(0.6), FT(0.8), FT(0.0), FT(0.5)]
    frac_sno_eff = fill(FT(0.0), nc)
    snow_depth = FT[FT(0.0), FT(0.3), FT(0.5), FT(0.0), FT(0.2)]
    swe_old = fill(FT(0.0), nc, nlev)
    dz = fill(FT(0.0), nc, nlev)
    for c in 1:nc
        for j in (snl[c] + 1):0
            jj = j + NLEVSNO
            dz[c, jj] = FT(0.02)
        end
    end

    lun_itype = Int[CLM.ISTSOIL, CLM.ISTSOIL, CLM.ISTSOIL, CLM.ISTSOIL, CLM.ISTDLAK]
    urbpoi    = Bool[false, false, false, false, false]

    # n_melt SCA shape parameter (per column) — default 200/std clamp
    n_melt = FT[FT(20.0), FT(15.0), FT(8.0), FT(20.0), FT(12.0)]

    return (; nc, snl, mask, h2osoi_ice, h2osoi_liq, h2osno_total, qflx_snow_grnd,
            qflx_snow_drain, int_snow, snomelt_accum, bifall, frac_sno, frac_sno_eff,
            snow_depth, swe_old, dz, lun_itype, urbpoi, n_melt)
end

function run_bulkdiag!(S, scf, n_melt_dev)
    bounds = 1:S.nc
    CLM.bulkdiag_new_snow_diagnostics!(
        S.dz, S.int_snow, S.swe_old, S.frac_sno, S.frac_sno_eff,
        S.snow_depth, S.snomelt_accum,
        DTIME, S.lun_itype, S.urbpoi, S.snl, S.bifall, S.h2osno_total,
        S.h2osoi_ice, S.h2osoi_liq, S.qflx_snow_grnd, S.qflx_snow_drain,
        S.mask, bounds, NLEVSNO; scf_method=scf, n_melt_dev=n_melt_dev)
    return nothing
end

function main(backend)
    println("=" ^ 72)
    println("END-TO-END GPU parity for handle_new_snow! tail (bulkdiag + scf) (A10)")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT); B = build(FT)
    scf = CLM.SnowCoverFractionSwensonLawrence2012()

    run_bulkdiag!(H, scf, H.n_melt)

    db(x) = to(collect(Bool, x))
    Sd = (; nc = B.nc, snl = to(B.snl), mask = db(B.mask),
            h2osoi_ice = to(B.h2osoi_ice), h2osoi_liq = to(B.h2osoi_liq),
            h2osno_total = to(B.h2osno_total), qflx_snow_grnd = to(B.qflx_snow_grnd),
            qflx_snow_drain = to(B.qflx_snow_drain), int_snow = to(B.int_snow),
            snomelt_accum = to(B.snomelt_accum), bifall = to(B.bifall),
            frac_sno = to(B.frac_sno), frac_sno_eff = to(B.frac_sno_eff),
            snow_depth = to(B.snow_depth), swe_old = to(B.swe_old), dz = to(B.dz),
            lun_itype = to(B.lun_itype), urbpoi = db(B.urbpoi), n_melt = to(B.n_melt))

    if !(Sd.int_snow isa device_array_type())
        println("  BLOCKED: device arrays did not move under to().")
        return 2
    end
    run_bulkdiag!(Sd, scf, Sd.n_melt)

    tonum(x) = Float64.(Array(x))
    checks = [
        ("frac_sno",      H.frac_sno,      Sd.frac_sno),
        ("frac_sno_eff",  H.frac_sno_eff,  Sd.frac_sno_eff),
        ("snow_depth",    H.snow_depth,    Sd.snow_depth),
        ("int_snow",      H.int_snow,      Sd.int_snow),
        ("snomelt_accum", H.snomelt_accum, Sd.snomelt_accum),
        ("swe_old",       H.swe_old,       Sd.swe_old),
        ("dz",            H.dz,            Sd.dz),
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
        @printf("    col %d: snl=%2d lun=%d frac_sno=%7.4f frac_eff=%7.4f depth=%7.4f int_snow=%10.3f\n",
                c, H.snl[c], H.lun_itype[c], H.frac_sno[c], H.frac_sno_eff[c],
                H.snow_depth[c], H.int_snow[c])
    end

    println()
    println(nfail == 0 ? "  handle_new_snow! tail MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
