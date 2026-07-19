# ==========================================================================
# gpu_validate_combinedivide_e2e.jl — end-to-end GPU parity for the WHOLE
# combine_snow_layers! + divide_snow_layers! drivers (A7, the hardest mechanical
# snow cluster): one thread per column running the full sequential, loop-carried
# snow-layer restructuring in-thread on its own fixed-max PADDED slice —
#   combine: thin-layer removal (water/aerosol transfer + shift), all-snow-gone
#            collapse to h2osno_no_layers, and pairwise thin-layer merging (combo!).
#   divide:  top-to-bottom split of layers thicker than their dzmax (combo_scalar +
#            mass_weighted_snow_radius), via per-column device-resident scratch rows.
#
# Builds a small Float32 instance with columns at DIFFERENT snl so every layer-count
# branch executes on BOTH backends:
#   col 1: snl=-1 thin single layer (combine-away candidate)
#   col 2: snl=-2 poised to combine (thin lower layer merges)
#   col 3: snl=-1 thick single layer (divide candidate -> snl=-2 or -3)
#   col 4: snl=-3 deep multi-layer
#   col 5: snl=0  no snow (mask false — untouched)
#
# Runs combine then divide on the CPU, adapts every array (+ aerosol + mask) to Metal,
# runs the SAME calls on the device, and compares snl AND every layer array field-by-
# field. CPU-reference fields are asserted finite so a both-NaN false PASS cannot slip.
#
#   julia --project=scripts scripts/gpu_validate_combinedivide_e2e.jl
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
cpu_all_finite(a) = all(isfinite, Array(a))

const NS = 5   # nlevsno (fixed max padded snow layers)
const NMAX = NS + 10   # total padded columns (snow + soil), generous

# Initialize the module-level snow-thickness threshold arrays (host Float64 Refs),
# mirroring init_snow_layers!.
function set_snow_thresholds!()
    nlevsno = NS
    dzmin  = zeros(Float64, nlevsno)
    dzmax_l = zeros(Float64, nlevsno)
    dzmax_u = zeros(Float64, nlevsno)
    dzmin[1] = 0.010; dzmax_l[1] = 0.03; dzmax_u[1] = 0.02
    dzmin[2] = 0.015; dzmax_l[2] = 0.07; dzmax_u[2] = 0.05
    for j in 3:nlevsno
        dzmin[j]  = dzmax_u[j-1] * 0.5
        dzmax_u[j] = 2.0 * dzmax_u[j-1] + 0.01
        dzmax_l[j] = dzmax_u[j] + dzmax_l[j-1]
        if j == nlevsno
            dzmax_u[j] = floatmax(Float64)
            dzmax_l[j] = floatmax(Float64)
        end
    end
    CLM.SNOW_DZMIN[]   = dzmin
    CLM.SNOW_DZMAX_L[] = dzmax_l
    CLM.SNOW_DZMAX_U[] = dzmax_u
    return nothing
end

# Build the column/layer arrays + aerosol at precision FT.
function build(::Type{FT}) where {FT}
    nc = 5
    snl   = zeros(Int, nc)
    dz    = zeros(FT, nc, NMAX)
    zi    = zeros(FT, nc, NMAX + 1)
    z     = zeros(FT, nc, NMAX)
    t     = zeros(FT, nc, NMAX)
    ice   = zeros(FT, nc, NMAX)
    liq   = zeros(FT, nc, NMAX)
    snwrds = zeros(FT, nc, NMAX)
    h2osno_no = zeros(FT, nc)
    snow_depth = zeros(FT, nc)
    frac_sno = ones(FT, nc)
    frac_sno_eff = ones(FT, nc)
    int_snow = zeros(FT, nc)

    aer = CLM.AerosolData{FT}()
    aer.mss_bcphi_col = zeros(FT, nc, NMAX); aer.mss_bcpho_col = zeros(FT, nc, NMAX)
    aer.mss_ocphi_col = zeros(FT, nc, NMAX); aer.mss_ocpho_col = zeros(FT, nc, NMAX)
    aer.mss_dst1_col = zeros(FT, nc, NMAX); aer.mss_dst2_col = zeros(FT, nc, NMAX)
    aer.mss_dst3_col = zeros(FT, nc, NMAX); aer.mss_dst4_col = zeros(FT, nc, NMAX)

    lun_itype = fill(CLM.ISTSOIL, 1)
    urbpoi = fill(false, 1)
    col_landunit = ones(Int, nc)

    # Fill a snow column: snl layers [snl+1..0], plus give the top soil layer (idx NS+1)
    # nonzero water so combine's "transfer liquid to soil layer 1" path is well-defined.
    function fill_col!(c, snl_c, dzs::Vector{Float64}, ices::Vector{Float64}, liqs::Vector{Float64})
        snl[c] = snl_c
        n = -snl_c
        for k in 1:n
            jj = (k + snl_c) + NS  # layer index k counts from top
            dz[c, jj]  = FT(dzs[k])
            ice[c, jj] = FT(ices[k])
            liq[c, jj] = FT(liqs[k])
            t[c, jj]   = FT(270.0)
            snwrds[c, jj] = FT(100.0)
            aer.mss_bcphi_col[c, jj] = FT(1.0e-9 * k)
            aer.mss_bcpho_col[c, jj] = FT(2.0e-9 * k)
            aer.mss_ocphi_col[c, jj] = FT(3.0e-9 * k)
            aer.mss_ocpho_col[c, jj] = FT(4.0e-9 * k)
            aer.mss_dst1_col[c, jj] = FT(5.0e-9 * k)
            aer.mss_dst2_col[c, jj] = FT(6.0e-9 * k)
            aer.mss_dst3_col[c, jj] = FT(7.0e-9 * k)
            aer.mss_dst4_col[c, jj] = FT(8.0e-9 * k)
        end
        # top soil layer (Fortran layer 1) liquid sink
        liq[c, NS + 1] = FT(10.0)
        # surface interface depth = 0; build zi downward for snow + into soil
        zi[c, NS + 1] = FT(0.0)
        for j in 1:NMAX
            zi[c, j + 1] = zi[c, j] + dz[c, j]
        end
        for j in 1:NMAX
            z[c, j] = zi[c, j + 1] - FT(0.5) * dz[c, j]
        end
    end

    # col 1: snl=-1 thin single layer -> combine-away (very little ice, frac_sno_eff*depth tiny)
    fill_col!(1, -1, [0.005], [0.001], [0.0])
    # col 2: snl=-2 poised to combine (lower layer thin/light -> merges to -1)
    fill_col!(2, -2, [0.05, 0.004], [10.0, 0.5], [1.0, 0.05])
    # col 3: snl=-1 thick single layer -> divide (exceeds dzmax_l[1]=0.03)
    fill_col!(3, -1, [0.20], [40.0], [4.0])
    # col 4: snl=-3 deep multi-layer (all healthy)
    fill_col!(4, -3, [0.02, 0.05, 0.15], [4.0, 10.0, 30.0], [0.4, 1.0, 3.0])
    # col 5: snl=0 no snow (mask false)
    snl[5] = 0

    return (; nc, snl, dz, zi, z, t, ice, liq, snwrds, h2osno_no, snow_depth,
            frac_sno, frac_sno_eff, int_snow, aer, lun_itype, urbpoi, col_landunit,
            qflx_sl_top_soil = zeros(nc))
end

function run_cpu!(B, mask, bounds)
    CLM.combine_snow_layers!(B.snl, B.dz, B.zi, B.z, B.t, B.ice, B.liq, B.h2osno_no,
        B.snow_depth, B.frac_sno, B.frac_sno_eff, B.int_snow, B.snwrds, B.aer,
        B.lun_itype, B.urbpoi, B.col_landunit, mask, bounds, NS,
        B.qflx_sl_top_soil, 1800.0)
    CLM.divide_snow_layers!(B.snl, B.dz, B.zi, B.z, B.t, B.ice, B.liq,
        B.frac_sno, B.snwrds, B.aer, false, mask, bounds, NS)
end

function run_dev!(D, mask, bounds)
    CLM.combine_snow_layers!(D.snl, D.dz, D.zi, D.z, D.t, D.ice, D.liq, D.h2osno_no,
        D.snow_depth, D.frac_sno, D.frac_sno_eff, D.int_snow, D.snwrds, D.aer,
        D.lun_itype, D.urbpoi, D.col_landunit, mask, bounds, NS,
        D.qflx_sl_top_soil, 1800.0)
    CLM.divide_snow_layers!(D.snl, D.dz, D.zi, D.z, D.t, D.ice, D.liq,
        D.frac_sno, D.snwrds, D.aer, false, mask, bounds, NS)
end

function main(backend)
    println("=" ^ 72)
    println("END-TO-END GPU parity for combine_snow_layers! + divide_snow_layers!")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    saved = (CLM.varpar.nlevsno,)
    CLM.varpar.nlevsno = NS
    set_snow_thresholds!()
    try
        B = build(FT)
        nc = B.nc
        bounds = 1:nc
        m_cpu = falses(nc); for c in 1:4; m_cpu[c] = true; end   # col 5 (snl=0) excluded

        # Device snapshot of the populated initial state BEFORE the CPU run mutates B.
        ad(x) = CLM.Adapt.adapt(device_array_type(), x)
        D = (; nc,
            snl = ad(B.snl), dz = ad(B.dz), zi = ad(B.zi), z = ad(B.z), t = ad(B.t),
            ice = ad(B.ice), liq = ad(B.liq), snwrds = ad(B.snwrds),
            h2osno_no = ad(B.h2osno_no), snow_depth = ad(B.snow_depth),
            frac_sno = ad(B.frac_sno), frac_sno_eff = ad(B.frac_sno_eff),
            int_snow = ad(B.int_snow), aer = ad(B.aer),
            lun_itype = ad(B.lun_itype), urbpoi = ad(B.urbpoi),
            col_landunit = ad(B.col_landunit),
            qflx_sl_top_soil = ad(B.qflx_sl_top_soil))
        m_dev = to(collect(Bool, m_cpu))

        # sanity: arrays must actually be on the device
        if !(D.dz isa device_array_type())
            println("  BLOCKED: dz did not move to the device under adapt.")
            return 2
        end

        # CPU reference run
        run_cpu!(B, m_cpu, bounds)
        # Assert CPU reference fields finite (no unset-input NaN that would mask a divergence)
        for (nm, a) in (("dz", B.dz), ("zi", B.zi), ("z", B.z), ("t", B.t),
                        ("ice", B.ice), ("liq", B.liq), ("snwrds", B.snwrds),
                        ("snow_depth", B.snow_depth), ("h2osno_no", B.h2osno_no))
            if !cpu_all_finite(a)
                @printf("  WARN: CPU reference field %s contains non-finite values!\n", nm)
            end
        end
        println("  post-call snl (CPU): ", B.snl)

        # Device run
        run_dev!(D, m_dev, bounds)
        println("  post-call snl (dev): ", Array(D.snl))

        checks = [
            ("snl",        B.snl,        D.snl),
            ("dz",         B.dz,         D.dz),
            ("zi",         B.zi,         D.zi),
            ("z",          B.z,          D.z),
            ("t_soisno",   B.t,          D.t),
            ("h2osoi_ice", B.ice,        D.ice),
            ("h2osoi_liq", B.liq,        D.liq),
            ("snw_rds",    B.snwrds,     D.snwrds),
            ("snow_depth", B.snow_depth, D.snow_depth),
            ("h2osno_no",  B.h2osno_no,  D.h2osno_no),
            ("frac_sno",   B.frac_sno,   D.frac_sno),
            ("frac_sno_eff", B.frac_sno_eff, D.frac_sno_eff),
            ("int_snow",   B.int_snow,   D.int_snow),
            ("mss_bcphi",  B.aer.mss_bcphi_col, D.aer.mss_bcphi_col),
            ("mss_bcpho",  B.aer.mss_bcpho_col, D.aer.mss_bcpho_col),
            ("mss_ocphi",  B.aer.mss_ocphi_col, D.aer.mss_ocphi_col),
            ("mss_ocpho",  B.aer.mss_ocpho_col, D.aer.mss_ocpho_col),
            ("mss_dst1",   B.aer.mss_dst1_col,  D.aer.mss_dst1_col),
            ("mss_dst2",   B.aer.mss_dst2_col,  D.aer.mss_dst2_col),
            ("mss_dst3",   B.aer.mss_dst3_col,  D.aer.mss_dst3_col),
            ("mss_dst4",   B.aer.mss_dst4_col,  D.aer.mss_dst4_col),
        ]
        nfail = 0
        for (nm, a, b) in checks
            d = reldiff(a, b)
            ok = d < 5e-3
            @printf("  [%s] %-13s  reldiff = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
            ok || (nfail += 1)
        end
        println()
        println(nfail == 0 ?
            "  WHOLE combine_snow_layers! + divide_snow_layers! MATCH CPU ON $name ($FT)" :
            "  DIVERGENCE — investigate.")
        return nfail == 0 ? 0 : 1
    finally
        CLM.varpar.nlevsno, = saved
    end
end

const BACKEND = detect_backend()
exit(main(BACKEND))
