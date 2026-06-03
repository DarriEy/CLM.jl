# ==========================================================================
# gpu_validate_ch4diag_e2e.jl — end-to-end Metal parity for the CH4 "diagnostic"
# function-group of methane.jl:
#   • get_jwt!           — water-table layer search (per-column, in-thread j-loops)
#   • ch4_annualupdate!  — annual mean update (3 kernels: counter / col / patch /
#                          reset, all own-index)
#   • ch4_totcolch4!     — total column CH4 (per-column reduction over levels)
#
# These mutate a handful of CH4Data fields. CH4Data in this worktree is the
# concrete {FT} struct (not yet @adapt_structure'd), so this harness moves the
# touched fields (and the loose array args) to the device individually with the
# backend converter, exercises the SAME public calls on CPU vs device, and
# compares the mutated outputs field-by-field. The reference fields are asserted
# finite so a both-NaN false PASS cannot slip through.
#
# Two column geometries are exercised so every branch runs:
#   col 1: fully unsaturated         -> jwt == nlevsoi
#   col 2: bottom layers saturated   -> jwt at the water-table boundary
#   col 3: top layer saturated       -> jwt == 0
# annualupdate is run with the counter both below and at the year threshold; and
# totcolch4 is run with allowlakeprod true (lake mask active) and false.
#
#   julia --project=. scripts/gpu_validate_ch4diag_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

# reldiff: NaN-aware (both-NaN agrees; one-sided NaN flags divergence).
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

# Build a small CH4 instance at precision FT with deliberately diverse columns.
function build(::Type{FT}) where {FT}
    nc = 3; np = 4; ng = 2; nlevsoi = 5

    ch4 = CLM.CH4Data{FT}(
        conc_ch4_sat_col   = fill(FT(1.0e-4), nc, nlevsoi),
        conc_ch4_unsat_col = fill(FT(1.0e-5), nc, nlevsoi),
        finundated_col     = fill(FT(0.1), nc),

        totcolch4_col      = zeros(FT, nc),
        totcolch4_bef_col  = zeros(FT, nc),

        annsum_counter_col = zeros(FT, nc),
        tempavg_somhr_col  = zeros(FT, nc),
        annavg_somhr_col   = fill(FT(1.0e-6), nc),
        tempavg_finrw_col  = zeros(FT, nc),
        annavg_finrw_col   = fill(FT(0.1), nc),

        annavg_agnpp_patch  = fill(FT(1.0e-5), np),
        annavg_bgnpp_patch  = fill(FT(1.0e-5), np),
        tempavg_agnpp_patch = fill(FT(1.0e-6), np),
        tempavg_bgnpp_patch = fill(FT(1.0e-6), np),
    )

    params = CLM.CH4Params()

    mask_soil   = trues(nc)
    mask_soilp  = trues(np)
    mask_nolake = Bool[true, true, false]   # col 3 is a lake column
    mask_lake   = Bool[false, false, true]

    watsat = fill(FT(0.45), nc, nlevsoi)
    t_soisno = fill(FT(CLM.TFRZ + 15.0), nc, nlevsoi)
    dz = fill(FT(0.1), nc, nlevsoi)

    # h2osoi_vol: per-column branch coverage for get_jwt!
    h2osoi_vol = fill(FT(0.3), nc, nlevsoi)        # col 1: all unsaturated -> jwt==nlevsoi
    h2osoi_vol[2, 4] = FT(0.44); h2osoi_vol[2, 5] = FT(0.44)  # col 2: bottom sat -> boundary
    h2osoi_vol[3, 1] = FT(0.44)                    # col 3: top sat -> jwt==0

    patch_column = [1, 1, 2, 3]
    is_fates = falses(nc)
    somhr = fill(FT(1.0e-6), nc)
    agnpp = fill(FT(1.0e-5), np)
    bgnpp = fill(FT(1.0e-5), np)
    dt = FT(1800.0)
    secsperyear = FT(365.0 * 86400.0)

    return (; nc, np, ng, nlevsoi, ch4, params, mask_soil, mask_soilp,
            mask_nolake, mask_lake, watsat, t_soisno, dz, h2osoi_vol,
            patch_column, is_fates, somhr, agnpp, bgnpp, dt, secsperyear)
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for CH4 diagnostics (get_jwt!/annualupdate!/totcolch4!)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT)   # CPU reference
    B = build(FT)   # device source (independent copy)

    # Device-array helpers: Bool masks must be Vector{Bool} on the device.
    dvec(x) = to(collect(x))
    dmask(m) = to(collect(Bool, m))

    nfail = 0
    function check(nm, a, b)
        if !cpu_has_finite(a)
            @printf("  [WARN] %-22s CPU reference all-NaN/Inf — skipping (no signal)\n", nm)
            return
        end
        d = reldiff(a, b)
        ok = d < (FT === Float32 ? 1f-3 : 1e-10)
        @printf("  [%s] %-22s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end

    # ---------------- get_jwt! ----------------
    jwt_h = zeros(Int, H.nc)
    CLM.get_jwt!(jwt_h, H.mask_soil, H.watsat, H.h2osoi_vol, H.t_soisno, H.nlevsoi, H.params)

    jwt_d = to(zeros(Int, B.nc))
    CLM.get_jwt!(jwt_d, dmask(B.mask_soil), dvec(B.watsat), dvec(B.h2osoi_vol),
                 dvec(B.t_soisno), B.nlevsoi, B.params)

    if !(jwt_d isa typeof(to(zeros(Int, 1))))
        println("  BLOCKED: jwt did not move to the device under the converter.")
        return 2
    end
    check("get_jwt", jwt_h, jwt_d)
    @printf("    jwt_cpu = %s   jwt_dev = %s\n", string(Array(jwt_h)), string(Array(jwt_d)))

    # ---------------- ch4_totcolch4! (allowlakeprod true) ----------------
    tot_h = zeros(FT, H.nc)
    CLM.ch4_totcolch4!(tot_h, H.ch4, H.mask_nolake, H.mask_lake, H.dz, H.nlevsoi, true)

    tot_d = to(zeros(FT, B.nc))
    ch4_d_tot = CLM.CH4Data{FT}(
        finundated_col     = dvec(B.ch4.finundated_col),
        conc_ch4_sat_col   = dvec(B.ch4.conc_ch4_sat_col),
        conc_ch4_unsat_col = dvec(B.ch4.conc_ch4_unsat_col),
    )
    CLM.ch4_totcolch4!(tot_d, ch4_d_tot, dmask(B.mask_nolake), dmask(B.mask_lake),
                       dvec(B.dz), B.nlevsoi, true)
    check("totcolch4 (lakeprod)", tot_h, tot_d)

    # allowlakeprod false: lake column should be zeroed
    tot_h2 = zeros(FT, H.nc)
    CLM.ch4_totcolch4!(tot_h2, H.ch4, H.mask_nolake, H.mask_lake, H.dz, H.nlevsoi, false)
    tot_d2 = to(zeros(FT, B.nc))
    CLM.ch4_totcolch4!(tot_d2, ch4_d_tot, dmask(B.mask_nolake), dmask(B.mask_lake),
                       dvec(B.dz), B.nlevsoi, false)
    check("totcolch4 (no lake)", tot_h2, tot_d2)

    # ---------------- ch4_annualupdate! (below threshold) ----------------
    CLM.ch4_annualupdate!(H.ch4, H.mask_soil, H.mask_soilp, H.patch_column,
                          H.is_fates, H.somhr, H.agnpp, H.bgnpp, H.dt, H.secsperyear)

    ch4_d = CLM.CH4Data{FT}(
        annsum_counter_col  = dvec(B.ch4.annsum_counter_col),
        annavg_somhr_col    = dvec(B.ch4.annavg_somhr_col),
        tempavg_somhr_col   = dvec(B.ch4.tempavg_somhr_col),
        annavg_finrw_col    = dvec(B.ch4.annavg_finrw_col),
        tempavg_finrw_col   = dvec(B.ch4.tempavg_finrw_col),
        finundated_col      = dvec(B.ch4.finundated_col),
        annavg_agnpp_patch  = dvec(B.ch4.annavg_agnpp_patch),
        tempavg_agnpp_patch = dvec(B.ch4.tempavg_agnpp_patch),
        annavg_bgnpp_patch  = dvec(B.ch4.annavg_bgnpp_patch),
        tempavg_bgnpp_patch = dvec(B.ch4.tempavg_bgnpp_patch),
    )
    CLM.ch4_annualupdate!(ch4_d, dmask(B.mask_soil), dmask(B.mask_soilp),
                          dvec(B.patch_column), dmask(B.is_fates), dvec(B.somhr),
                          dvec(B.agnpp), dvec(B.bgnpp), B.dt, B.secsperyear)

    check("annsum_counter",  H.ch4.annsum_counter_col,  ch4_d.annsum_counter_col)
    check("tempavg_somhr",   H.ch4.tempavg_somhr_col,   ch4_d.tempavg_somhr_col)
    check("tempavg_finrw",   H.ch4.tempavg_finrw_col,   ch4_d.tempavg_finrw_col)
    check("tempavg_agnpp",   H.ch4.tempavg_agnpp_patch, ch4_d.tempavg_agnpp_patch)
    check("tempavg_bgnpp",   H.ch4.tempavg_bgnpp_patch, ch4_d.tempavg_bgnpp_patch)

    # ---------------- ch4_annualupdate! (at/over threshold -> the >=year branch) ----------------
    H2 = build(FT); B2 = build(FT)
    H2.ch4.annsum_counter_col .= FT(2.0) .* H2.secsperyear     # force the year-elapsed branch
    B2.ch4.annsum_counter_col .= FT(2.0) .* B2.secsperyear
    H2.ch4.tempavg_somhr_col  .= FT(3.0e-6)
    B2.ch4.tempavg_somhr_col  .= FT(3.0e-6)
    H2.ch4.tempavg_finrw_col  .= FT(4.0e-7)
    B2.ch4.tempavg_finrw_col  .= FT(4.0e-7)

    CLM.ch4_annualupdate!(H2.ch4, H2.mask_soil, H2.mask_soilp, H2.patch_column,
                          H2.is_fates, H2.somhr, H2.agnpp, H2.bgnpp, H2.dt, H2.secsperyear)

    ch4_d2 = CLM.CH4Data{FT}(
        annsum_counter_col  = dvec(B2.ch4.annsum_counter_col),
        annavg_somhr_col    = dvec(B2.ch4.annavg_somhr_col),
        tempavg_somhr_col   = dvec(B2.ch4.tempavg_somhr_col),
        annavg_finrw_col    = dvec(B2.ch4.annavg_finrw_col),
        tempavg_finrw_col   = dvec(B2.ch4.tempavg_finrw_col),
        finundated_col      = dvec(B2.ch4.finundated_col),
        annavg_agnpp_patch  = dvec(B2.ch4.annavg_agnpp_patch),
        tempavg_agnpp_patch = dvec(B2.ch4.tempavg_agnpp_patch),
        annavg_bgnpp_patch  = dvec(B2.ch4.annavg_bgnpp_patch),
        tempavg_bgnpp_patch = dvec(B2.ch4.tempavg_bgnpp_patch),
    )
    CLM.ch4_annualupdate!(ch4_d2, dmask(B2.mask_soil), dmask(B2.mask_soilp),
                          dvec(B2.patch_column), dmask(B2.is_fates), dvec(B2.somhr),
                          dvec(B2.agnpp), dvec(B2.bgnpp), B2.dt, B2.secsperyear)

    check("annavg_somhr (year)",  H2.ch4.annavg_somhr_col,  ch4_d2.annavg_somhr_col)
    check("annavg_finrw (year)",  H2.ch4.annavg_finrw_col,  ch4_d2.annavg_finrw_col)
    check("annsum_reset (year)",  H2.ch4.annsum_counter_col, ch4_d2.annsum_counter_col)
    check("annavg_agnpp (year)",  H2.ch4.annavg_agnpp_patch, ch4_d2.annavg_agnpp_patch)
    check("annavg_bgnpp (year)",  H2.ch4.annavg_bgnpp_patch, ch4_d2.annavg_bgnpp_patch)

    println()
    println(nfail == 0 ? "  CH4 DIAGNOSTICS MATCH CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
