# ==========================================================================
# gpu_validate_snowpercolation_e2e.jl — end-to-end Metal parity for the WHOLE
# snow-percolation chain (A7), run on the CPU and on the device and compared.
#
# Functions exercised (all whole-function, in the order they fire in the driver):
#   update_state_top_layer_fluxes!         (top snow layer dew/evap/rain input)
#   bulk_flux_snow_percolation!            (per-column sequential layer loop)
#   update_state_snow_percolation!         (per-(c,j) percolation transport)
#   calc_and_apply_aerosol_fluxes!         (per-column sequential, loop-carried
#                                           qin for BC/OC/dust species)
#   post_percolation_adjust_layer_thicknesses!
#   bulkdiag_snow_water_accumulated_snow!  (two-mask per-column)
#   sum_flux_add_snow_percolation!         (two-mask per-column)
#
# Columns span snl = 0 (no snow), -1, -2, -3 so the variable-layer-count branches
# (top-only, interior+bottom, deep multi-layer) all execute, plus a no-snow column
# to drive the mask_nosnow branches in the two diagnostic functions. Each snow
# column carries non-zero BC/OC/dust mass so aerosol transport down the column is
# genuinely exercised.
#
#   julia --project=scripts scripts/gpu_validate_snowpercolation_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

const NS = 5  # nlevsno

# NaN-aware reldiff: untouched padded layers stay NaN on both sides → agreement.
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        x = Float64(A[i]); y = Float64(B[i])
        (isnan(x) && isnan(y)) && continue
        d = abs(x - y); s = max(abs(x), abs(y))
        m = max(m, s > 1e-30 ? d / s : d)
    end
    return m
end

# assert the CPU reference is finite where it must be (per the banked A1 lesson:
# a both-NaN reldiff PASSES silently). Returns count of unexpected NaNs.
function count_bad(a, mask, snl, jrange)
    A = Array(a); n = 0
    for c in eachindex(mask)
        mask[c] || continue
        for j in (snl[c] + 1):0
            jj = j + NS
            isfinite(A[c, jj]) || (n += 1)
        end
    end
    return n
end

function build(::Type{FT}) where {FT}
    nc = 4
    snl = [0, -1, -2, -3]          # no-snow, 1, 2, 3 snow layers
    ntot = NS + 15
    dz          = zeros(FT, nc, ntot)
    h2osoi_ice  = zeros(FT, nc, ntot)
    h2osoi_liq  = zeros(FT, nc, ntot)
    qflx_snow_percolation = zeros(FT, nc, NS)
    frac_sno_eff = FT[0.0, 0.9, 0.85, 0.8]

    # Per-(layer) snow mass/thickness for the snow columns.
    for c in 2:nc
        for j in (snl[c] + 1):0
            jj = j + NS
            dz[c, jj]         = FT(0.02 + 0.01 * (j + 4))
            h2osoi_ice[c, jj] = FT(4.0 + 0.5 * (j + 4))
            h2osoi_liq[c, jj] = FT(2.0 + 0.3 * (j + 4))
        end
    end

    aer = CLM.AerosolData{FT}()
    CLM.aerosol_init!(aer, nc)
    # Zero the snow-layer aerosol mass everywhere, then seed the snow columns.
    for f in (:mss_bcphi_col, :mss_bcpho_col, :mss_ocphi_col, :mss_ocpho_col,
              :mss_dst1_col, :mss_dst2_col, :mss_dst3_col, :mss_dst4_col)
        fill!(getfield(aer, f), zero(FT))
    end
    for c in 2:nc
        for j in (snl[c] + 1):0
            jj = j + NS
            aer.mss_bcphi_col[c, jj] = FT(1.0e-7)
            aer.mss_bcpho_col[c, jj] = FT(2.0e-7)
            aer.mss_ocphi_col[c, jj] = FT(3.0e-7)
            aer.mss_ocpho_col[c, jj] = FT(4.0e-7)
            aer.mss_dst1_col[c, jj]  = FT(5.0e-7)
            aer.mss_dst2_col[c, jj]  = FT(6.0e-7)
            aer.mss_dst3_col[c, jj]  = FT(7.0e-7)
            aer.mss_dst4_col[c, jj]  = FT(8.0e-7)
        end
    end

    # Top-layer flux inputs (per column).
    qflx_soliddew_to_top_layer  = FT[0.0, 1.0e-5, 2.0e-5, 1.5e-5]
    qflx_solidevap_from_top_layer = FT[0.0, 5.0e-6, 4.0e-6, 3.0e-6]
    qflx_liq_grnd               = FT[1.0e-4, 2.0e-4, 1.5e-4, 1.2e-4]
    qflx_liqdew_to_top_layer    = FT[0.0, 3.0e-6, 2.0e-6, 1.0e-6]
    qflx_liqevap_from_top_layer = FT[0.0, 1.0e-6, 1.5e-6, 0.5e-6]

    # Diagnostic-function state.
    int_snow      = FT[10.0, 5.0, 8.0, 12.0]
    frac_sno      = FT[0.0, 0.9, 0.85, 0.8]
    snow_depth    = FT[0.0, 0.1, 0.2, 0.3]
    h2osno_no_layers = FT[0.0, 0.0, 0.0, 0.0]  # no-snow col → reset branch
    qflx_snow_drain        = FT[0.0, 1.0e-4, 2.0e-4, 3.0e-4]
    qflx_rain_plus_snomelt = zeros(FT, nc)
    qflx_snomelt           = FT[1.0e-5, 0.0, 0.0, 0.0]

    return (; nc, snl, dz, h2osoi_ice, h2osoi_liq, qflx_snow_percolation,
            frac_sno_eff, aer, qflx_soliddew_to_top_layer, qflx_solidevap_from_top_layer,
            qflx_liq_grnd, qflx_liqdew_to_top_layer, qflx_liqevap_from_top_layer,
            int_snow, frac_sno, snow_depth, h2osno_no_layers, qflx_snow_drain,
            qflx_rain_plus_snomelt, qflx_snomelt)
end

# Run the whole chain on a state bundle B with the given mask types.
function run_chain!(B, snl, mask_snow, mask_nosnow, dtime)
    bounds = 1:B.nc
    CLM.update_state_top_layer_fluxes!(B.h2osoi_ice, B.h2osoi_liq, dtime, snl,
        B.frac_sno_eff, B.qflx_soliddew_to_top_layer, B.qflx_solidevap_from_top_layer,
        B.qflx_liq_grnd, B.qflx_liqdew_to_top_layer, B.qflx_liqevap_from_top_layer,
        mask_snow, bounds, NS)
    CLM.bulk_flux_snow_percolation!(B.qflx_snow_percolation, dtime, snl, B.dz,
        B.frac_sno_eff, B.h2osoi_ice, B.h2osoi_liq, mask_snow, bounds, NS)
    CLM.update_state_snow_percolation!(B.h2osoi_liq, dtime, snl,
        B.qflx_snow_percolation, mask_snow, bounds, NS)
    CLM.calc_and_apply_aerosol_fluxes!(B.aer, dtime, snl, B.h2osoi_ice, B.h2osoi_liq,
        B.qflx_snow_percolation, mask_snow, bounds, NS)
    CLM.post_percolation_adjust_layer_thicknesses!(B.dz, snl, B.h2osoi_ice,
        B.h2osoi_liq, mask_snow, bounds, NS)
    CLM.bulkdiag_snow_water_accumulated_snow!(B.int_snow, B.frac_sno, B.snow_depth,
        dtime, B.frac_sno_eff, B.qflx_soliddew_to_top_layer, B.qflx_liqdew_to_top_layer,
        B.qflx_liq_grnd, B.h2osno_no_layers, mask_snow, mask_nosnow, bounds)
    # bottom-layer percolation (Fortran layer 0 → Julia NS) into the summed flux.
    qflx_perc_bottom = B.qflx_snow_percolation[:, NS]
    CLM.sum_flux_add_snow_percolation!(B.qflx_snow_drain, B.qflx_rain_plus_snomelt,
        B.frac_sno_eff, qflx_perc_bottom, B.qflx_liq_grnd, B.qflx_snomelt,
        mask_snow, mask_nosnow, bounds)
    return nothing
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for the snow-percolation chain (A7)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU chain exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    vp = CLM.varpar
    saved = vp.nlevsno
    vp.nlevsno = NS
    try
        B = build(FT)
        snl = B.snl
        mask_snow   = [snl[c] < 0  for c in 1:B.nc]   # has explicit snow layers
        mask_nosnow = [snl[c] == 0 for c in 1:B.nc]   # no snow layers

        # Device snapshot BEFORE the CPU run mutates B.
        ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
        D = (; nc = B.nc,
            dz = ad(B.dz), h2osoi_ice = ad(B.h2osoi_ice), h2osoi_liq = ad(B.h2osoi_liq),
            qflx_snow_percolation = ad(B.qflx_snow_percolation),
            frac_sno_eff = ad(B.frac_sno_eff), aer = ad(B.aer),
            qflx_soliddew_to_top_layer = ad(B.qflx_soliddew_to_top_layer),
            qflx_solidevap_from_top_layer = ad(B.qflx_solidevap_from_top_layer),
            qflx_liq_grnd = ad(B.qflx_liq_grnd),
            qflx_liqdew_to_top_layer = ad(B.qflx_liqdew_to_top_layer),
            qflx_liqevap_from_top_layer = ad(B.qflx_liqevap_from_top_layer),
            int_snow = ad(B.int_snow), frac_sno = ad(B.frac_sno),
            snow_depth = ad(B.snow_depth), h2osno_no_layers = ad(B.h2osno_no_layers),
            qflx_snow_drain = ad(B.qflx_snow_drain),
            qflx_rain_plus_snomelt = ad(B.qflx_rain_plus_snomelt),
            qflx_snomelt = ad(B.qflx_snomelt))

        if !(D.aer.mss_bcphi_col isa Metal.MtlArray)
            println("  BLOCKED: AerosolData did not move to the device under adapt.")
            return 2
        end

        dmask(m) = to(collect(Bool, m))
        dt = FT(1800)

        # CPU reference.
        run_chain!(B, snl, BitVector(mask_snow), BitVector(mask_nosnow), dt)

        # Assert the CPU reference fields are finite where snow exists (catch a
        # silent both-NaN false PASS before trusting any reldiff).
        bad = 0
        bad += count_bad(B.h2osoi_liq, mask_snow, snl, nothing)
        bad += count_bad(B.h2osoi_ice, mask_snow, snl, nothing)
        bad += count_bad(B.aer.mss_bcphi_col, mask_snow, snl, nothing)
        bad += count_bad(B.aer.mss_dst4_col, mask_snow, snl, nothing)
        if bad > 0
            println("  BLOCKED: CPU reference has $bad unexpected NaN(s) in snow layers; ",
                    "harness inputs are under-specified — reldiff would falsely PASS.")
            return 2
        end

        # Device run — snl and masks live on the device too (read inside kernels).
        run_chain!(D, to(snl), dmask(mask_snow), dmask(mask_nosnow), dt)

        # Verify snl unchanged (these functions do not restructure layers).
        @assert snl == B.snl

        checks = [
            ("h2osoi_ice", B.h2osoi_ice, D.h2osoi_ice),
            ("h2osoi_liq", B.h2osoi_liq, D.h2osoi_liq),
            ("qflx_snow_percolation", B.qflx_snow_percolation, D.qflx_snow_percolation),
            ("dz",          B.dz,          D.dz),
            ("mss_bcphi",   B.aer.mss_bcphi_col, D.aer.mss_bcphi_col),
            ("mss_bcpho",   B.aer.mss_bcpho_col, D.aer.mss_bcpho_col),
            ("mss_ocphi",   B.aer.mss_ocphi_col, D.aer.mss_ocphi_col),
            ("mss_ocpho",   B.aer.mss_ocpho_col, D.aer.mss_ocpho_col),
            ("mss_dst1",    B.aer.mss_dst1_col,  D.aer.mss_dst1_col),
            ("mss_dst2",    B.aer.mss_dst2_col,  D.aer.mss_dst2_col),
            ("mss_dst3",    B.aer.mss_dst3_col,  D.aer.mss_dst3_col),
            ("mss_dst4",    B.aer.mss_dst4_col,  D.aer.mss_dst4_col),
            ("int_snow",    B.int_snow,    D.int_snow),
            ("frac_sno",    B.frac_sno,    D.frac_sno),
            ("snow_depth",  B.snow_depth,  D.snow_depth),
            ("qflx_snow_drain", B.qflx_snow_drain, D.qflx_snow_drain),
            ("qflx_rain_plus_snomelt", B.qflx_rain_plus_snomelt, D.qflx_rain_plus_snomelt),
        ]
        nfail = 0
        for (nm, a, b) in checks
            d = reldiff(a, b)
            ok = d < 1f-4
            @printf("  [%s] %-24s  reldiff = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
            ok || (nfail += 1)
        end
        println()
        println(nfail == 0 ? "  WHOLE snow-percolation chain MATCHES CPU ON $name ($FT) ✓" :
                             "  DIVERGENCE — investigate.")
        return nfail == 0 ? 0 : 1
    finally
        vp.nlevsno = saved
    end
end

const BACKEND = detect_backend()
exit(main(BACKEND))
