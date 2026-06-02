# ==========================================================================
# gpu_validate_snowcompaction_e2e.jl — end-to-end Metal parity for the WHOLE
# snow_compaction! driver (the per-column sequential snow-layer compaction).
#
# Builds a small set of snow columns with DIFFERENT snow-layer counts
# (snl = -1, -2, -3, plus a no-snow column), runs snow_compaction! on the CPU,
# then runs the SAME call with every array adapted to the Metal device, and
# compares the mutated `dz` (the only output) field.
#
# Each snow column is set up so the destructive-metamorphism, overburden, melt,
# AND wind-drift compaction branches are all exercised (imelt=1 on a non-lake /
# non-urban landunit, wind_dependent_snow_density = true, forc_wind > 0).
#
#   julia --project=. scripts/gpu_validate_snowcompaction_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

const NLEVSNO = 5

# reldiff that treats both-NaN as agreement (a one-sided NaN flags divergence),
# and is relative where the magnitude is meaningful.
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        x = Float64(A[i]); y = Float64(B[i])
        (isnan(x) && isnan(y)) && continue
        denom = max(abs(x), abs(y), 1e-12)
        m = max(m, abs(x - y) / denom)
    end
    return m
end

function assert_finite(name, a)
    A = Array(a)
    for i in eachindex(A)
        @assert isfinite(A[i]) "CPU reference $name has non-finite entry at $i = $(A[i])"
    end
end

# Build the column arrays for the compaction call at precision FT.
# nc columns: 1 no-snow (snl=0), then snl=-1, -2, -3.
function build(::Type{FT}) where {FT}
    nc = 4
    JOFF = NLEVSNO
    nlev_tot = NLEVSNO + 8   # snow layers + some soil layers (padded)

    snl          = Int[0, -1, -2, -3]
    col_gridcell = Int[1, 2, 3, 4]
    col_landunit = Int[1, 2, 3, 4]
    lakpoi       = Bool[false, false, false, false]
    urbpoi       = Bool[false, false, false, false]

    dz          = zeros(FT, nc, nlev_tot)
    t_soisno    = fill(FT(270.0), nc, nlev_tot)
    h2osoi_ice  = zeros(FT, nc, nlev_tot)
    h2osoi_liq  = zeros(FT, nc, nlev_tot)
    imelt       = zeros(Int, nc, nlev_tot)
    swe_old     = zeros(FT, nc, nlev_tot)
    frac_iceold = fill(FT(0.9), nc, nlev_tot)

    frac_sno    = FT[0.0, 0.8, 0.9, 0.95]
    frac_h2osfc = FT[0.0, 0.05, 0.05, 0.05]
    int_snow    = FT[0.0, 50.0, 120.0, 300.0]
    forc_wind   = FT[0.0, 6.0, 8.0, 10.0]

    # Populate snow layers [snl(c)+1 .. 0] for each snow column.
    for c in 2:nc
        for j in (snl[c] + 1):0
            jj = j + JOFF
            dz[c, jj]          = FT(0.05) + FT(0.02) * (j + NLEVSNO)  # thin-ish layers (m)
            t_soisno[c, jj]    = FT(268.0) + FT(1.0) * (j + NLEVSNO)  # below freezing, varies
            h2osoi_ice[c, jj]  = FT(8.0) + FT(2.0) * (j + NLEVSNO)    # > 0.1 → compaction active
            h2osoi_liq[c, jj]  = FT(1.0)                              # some liquid water term
            imelt[c, jj]       = 1                                    # melt-compaction branch
            swe_old[c, jj]     = h2osoi_ice[c, jj] + h2osoi_liq[c, jj] + FT(3.0)  # delta > 0
            frac_iceold[c, jj] = FT(0.92)
        end
    end

    return (; nc, snl, dz, t_soisno, h2osoi_ice, h2osoi_liq, imelt, frac_sno,
            frac_h2osfc, swe_old, int_snow, frac_iceold, forc_wind,
            col_gridcell, col_landunit, lakpoi, urbpoi)
end

function run_compaction!(B, mask, dtime)
    CLM.snow_compaction!(
        B.dz, dtime, B.snl, B.t_soisno, B.h2osoi_ice, B.h2osoi_liq, B.imelt,
        B.frac_sno, B.frac_h2osfc, B.swe_old, B.int_snow, B.frac_iceold,
        B.forc_wind, B.col_gridcell, B.col_landunit, B.lakpoi, B.urbpoi,
        mask, 1:B.nc, NLEVSNO)
end

function run_regime(name, backend, ob_method, wind_dep)
    name_b, to, FT = backend
    vp = CLM.varpar
    saved = vp.nlevsno
    vp.nlevsno = NLEVSNO
    sv_obm = CLM.OVERBURDEN_COMPACTION_METHOD[]
    sv_wds = CLM.WIND_DEPENDENT_SNOW_DENSITY[]
    CLM.OVERBURDEN_COMPACTION_METHOD[] = ob_method
    CLM.WIND_DEPENDENT_SNOW_DENSITY[] = wind_dep
    try
        B = build(FT)
        mask = falses(B.nc)
        for c in 2:B.nc; mask[c] = true; end   # all snow columns active
        snl_before = copy(B.snl)
        dt = FT(1800)

        # Device snapshot BEFORE the CPU run mutates dz (dz is the only output, so
        # only dz needs an independent copy; the inputs are read-only).
        D = (; nc = B.nc,
            snl = to(B.snl),
            dz = to(copy(B.dz)),
            t_soisno = to(B.t_soisno),
            h2osoi_ice = to(B.h2osoi_ice),
            h2osoi_liq = to(B.h2osoi_liq),
            imelt = to(B.imelt),
            frac_sno = to(B.frac_sno),
            frac_h2osfc = to(B.frac_h2osfc),
            swe_old = to(B.swe_old),
            int_snow = to(B.int_snow),
            frac_iceold = to(B.frac_iceold),
            forc_wind = to(B.forc_wind),
            col_gridcell = to(B.col_gridcell),
            col_landunit = to(B.col_landunit),
            lakpoi = to(collect(Bool, B.lakpoi)),
            urbpoi = to(collect(Bool, B.urbpoi)))
        dmask = to(collect(Bool, mask))

        # CPU reference run.
        run_compaction!(B, mask, dt)
        assert_finite("dz", B.dz)
        # snl must be unchanged by compaction.
        @assert B.snl == snl_before "snl changed during compaction (should not)"

        # Device run.
        run_compaction!(D, dmask, dt)

        d = reldiff(B.dz, D.dz)
        ok = d < 5e-5
        @printf("  [%s] regime %-28s  ob=%d wind=%-5s  max reldiff(dz) = %.3e\n",
                ok ? "PASS" : "FAIL", name, ob_method, string(wind_dep), d)
        return ok
    finally
        vp.nlevsno = saved
        CLM.OVERBURDEN_COMPACTION_METHOD[] = sv_obm
        CLM.WIND_DEPENDENT_SNOW_DENSITY[] = sv_wds
    end
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for snow_compaction! (whole per-column kernel)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)
    @printf("  Columns: snl = 0 (no snow), -1, -2, -3   nlevsno=%d\n\n", NLEVSNO)

    nfail = 0
    # destructive + overburden(Anderson) + melt + wind-drift
    nfail += run_regime("anderson + wind-drift", backend,
                        CLM.OVERBURDEN_COMPACTION_ANDERSON1976, true) ? 0 : 1
    # destructive + overburden(Vionnet) + melt + wind-drift
    nfail += run_regime("vionnet + wind-drift", backend,
                        CLM.OVERBURDEN_COMPACTION_VIONNET2012, true) ? 0 : 1
    # overburden(Anderson), no wind
    nfail += run_regime("anderson, no wind", backend,
                        CLM.OVERBURDEN_COMPACTION_ANDERSON1976, false) ? 0 : 1

    println()
    println(nfail == 0 ? "  WHOLE snow_compaction! MATCHES CPU ON $name ($FT) across snl=-1/-2/-3 ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
