# ==========================================================================
# gpu_validate_watertablevariants_e2e.jl — end-to-end GPU parity for the
# WHOLE perched_water_table! AND theta_based_water_table! drivers: one thread
# per hydrology column running the full nested, loop-carried water-table-depth
# search in-thread —
#   • perched: frost-table search (2→nlevsoi) + bottom-up saturation search
#     (k_frz→1) + interpolation branch (zwt_perched, frost_table).
#   • theta:   nbedrock→1 saturation search + k_zwt branch (zwt).
#
# Builds a small Float32 instance with SEVERAL columns deliberately placing the
# (perched / theta) water table at different depths so every branch runs:
#   col 1: shallow / fully-saturated upper column  -> water table near top
#   col 2: mid-column partially-saturated          -> interpolation branch
#   col 3: deep / mostly-dry column                -> water table deep
#   col 4: at/below bedrock (nbedrock < nlevsoi for theta; frozen upper layers
#          for the perched frost-table search)     -> bedrock + perched branches
# Runs each whole function on the CPU (Float64-equivalent reference at FT),
# adapts every state struct (+ mask + geometry) to the GPU, runs the SAME call
# on the device, and compares the mutated outputs field-by-field. CPU-reference
# fields are asserted finite first so a both-NaN false PASS cannot slip through.
#
#   julia --project=scripts scripts/gpu_validate_watertablevariants_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

# reldiff: NaN-aware (both-NaN agrees; one-sided NaN flags divergence).
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

# Build a small instance for the perched/theta water-table variants. Geometry
# uses direct [c,k] indexing (nc x nlevsoi), matching the host loops / tests.
function build(::Type{FT}) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevsoi = CLM.varpar.nlevsoi
    nc = 4; np = 4; nl = 1; ng = 1

    sh = CLM.SoilHydrologyData{FT}(); CLM.soilhydrology_init!(sh, nc)
    ss = CLM.SoilStateData{FT}();     CLM.soilstate_init!(ss, np, nc)
    ws = CLM.WaterStateData{FT}();    CLM.waterstate_init!(ws, nc, np, nl, ng)
    temp = CLM.TemperatureData{FT}(); CLM.temperature_init!(temp, np, nc, nl, ng)

    # --- geometry: direct [c,k] soil-layer thickness/depth/interface (m) ---
    col_dz = fill(FT(NaN), nc, nlevsoi)
    col_z  = fill(FT(NaN), nc, nlevsoi)
    col_zi = fill(FT(NaN), nc, nlevsoi)
    for c in 1:nc
        zi_prev = FT(0.0)
        for k in 1:nlevsoi
            dz = FT(0.05) * FT(1.4)^(k - 1)
            col_dz[c, k] = dz
            zi = zi_prev + dz
            col_zi[c, k] = zi
            col_z[c, k]  = FT(0.5) * (zi_prev + zi)
            zi_prev = zi
        end
    end
    zi_bot = Array(col_zi)[1, nlevsoi]

    # --- soil hydraulic properties (per column,layer) ---
    for c in 1:nc, k in 1:nlevsoi
        ss.watsat_col[c, k] = FT(0.45)
    end

    # --- temperature: warm everywhere by default ---
    for c in 1:nc, k in 1:nlevsoi
        temp.t_soisno_col[c, k] = FT(285.0)
    end

    # --- per-column water table scalars ---
    for c in 1:nc
        sh.zwt_col[c]         = FT(0.5 * zi_bot)
        sh.zwt_perched_col[c] = FT(0.0)
        sh.frost_table_col[c] = FT(0.0)
    end

    # Helper to set a column's moisture by per-layer saturation fraction.
    setcol! = (c, fracs) -> begin
        for k in 1:nlevsoi
            f = fracs(k)
            ws.h2osoi_liq_col[c, k] = f * FT(0.45) * Array(col_dz)[c, k] * FT(1000.0)
            ws.h2osoi_ice_col[c, k] = FT(0.0)
            ws.h2osoi_vol_col[c, k] = FT(NaN)
        end
    end

    # col 1: fully saturated whole column -> k_zwt==1, water table near surface.
    setcol!(1, k -> FT(1.0))
    # col 2: saturated bottom-up, sharp drop near layer 8 -> theta INTERPOLATION
    #        branch (k_zwt < nbedrock) with s1 < s2.
    setcol!(2, k -> k >= 8 ? FT(0.99) : FT(0.4))
    # col 3: mostly dry -> bottom-up search exits at nbedrock immediately
    #        -> water table deep (else branch / shallow bedrock).
    setcol!(3, k -> FT(0.4))
    # col 4: saturated upper, dry just below midcol -> perched interpolation.
    setcol!(4, k -> k <= 5 ? FT(0.99) : FT(0.3))

    # nbedrock: col 3 has a SHALLOW bedrock so theta hits the bedrock else-branch.
    col_nbedrock = fill(Int(nlevsoi), nc)
    col_nbedrock[3] = max(2, nlevsoi - 4)

    # Perched: col 4 has frozen layers below k=6 (so k_frz>=2, no index-0 OOB),
    # with a saturated→dry transition above the frost table for the perched search.
    for k in 1:nlevsoi
        temp.t_soisno_col[4, k] = k <= 5 ? FT(285.0) : FT(270.0)
    end
    # Keep zwt below the frost table for col 4 so the perched search runs.
    sh.zwt_col[4] = FT(0.95 * zi_bot)

    S = (; sh, ss, ws, temp, col_dz, col_z, col_zi, col_nbedrock)
    return (; nc, nlevsoi, S)
end

run_perched!(S, mask, nc, nlevsoi) = CLM.perched_water_table!(
    S.sh, S.ss, S.temp.t_soisno_col, S.ws,
    S.col_dz, S.col_z, S.col_zi, mask, 1:nc, nlevsoi)

run_theta!(S, mask, nc, nlevsoi) = CLM.theta_based_water_table!(
    S.sh, S.ss, S.ws, S.col_dz, S.col_z, S.col_zi,
    S.col_nbedrock, mask, 1:nc, nlevsoi)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for perched_water_table! + theta_based_water_table!")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU drivers exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    nfail = 0
    ad(x) = CLM.Adapt.adapt(device_array_type(), x)
    dmask(m) = to(collect(Bool, m))

    # ============================ perched_water_table! ====================
    begin
        H = build(FT); B = build(FT)
        Sd = (; sh = ad(B.S.sh), ss = ad(B.S.ss), ws = ad(B.S.ws), temp = ad(B.S.temp),
                col_dz = to(B.S.col_dz), col_z = to(B.S.col_z), col_zi = to(B.S.col_zi),
                col_nbedrock = to(B.S.col_nbedrock))
        if !(Sd.sh.zwt_perched_col isa device_array_type())
            println("  BLOCKED: perched state did not move to the device under adapt.")
            return 2
        end
        mask = trues(H.nc)
        run_perched!(H.S, mask, H.nc, H.nlevsoi)
        run_perched!(Sd, dmask(mask), B.nc, B.nlevsoi)

        checks = [
            ("zwt_perched", H.S.sh.zwt_perched_col, Sd.sh.zwt_perched_col),
            ("frost_table", H.S.sh.frost_table_col, Sd.sh.frost_table_col),
            ("h2osoi_vol",  H.S.ws.h2osoi_vol_col,  Sd.ws.h2osoi_vol_col),
        ]
        println("  --- perched_water_table! ---")
        for (nm, a, b) in checks
            @assert cpu_has_finite(a) "CPU reference for $nm is all-NaN/Inf — fix inputs"
            d = reldiff(a, b); ok = d < 1f-3
            @printf("  [%s] %-12s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
            ok || (nfail += 1)
        end
        zp_c = Array(H.S.sh.zwt_perched_col); zp_d = Array(Sd.sh.zwt_perched_col)
        ft_c = Array(H.S.sh.frost_table_col)
        for c in 1:H.nc
            @printf("    col %d: frost_table = %9.5f  zwt_perched_cpu = %9.5f  zwt_perched_dev = %9.5f\n",
                    c, ft_c[c], zp_c[c], zp_d[c])
        end
        println()
    end

    # ============================ theta_based_water_table! ================
    begin
        H = build(FT); B = build(FT)
        Sd = (; sh = ad(B.S.sh), ss = ad(B.S.ss), ws = ad(B.S.ws), temp = ad(B.S.temp),
                col_dz = to(B.S.col_dz), col_z = to(B.S.col_z), col_zi = to(B.S.col_zi),
                col_nbedrock = to(B.S.col_nbedrock))
        mask = trues(H.nc)
        run_theta!(H.S, mask, H.nc, H.nlevsoi)
        run_theta!(Sd, dmask(mask), B.nc, B.nlevsoi)

        checks = [
            ("zwt",        H.S.sh.zwt_col,        Sd.sh.zwt_col),
            ("h2osoi_vol", H.S.ws.h2osoi_vol_col, Sd.ws.h2osoi_vol_col),
        ]
        println("  --- theta_based_water_table! ---")
        for (nm, a, b) in checks
            @assert cpu_has_finite(a) "CPU reference for $nm is all-NaN/Inf — fix inputs"
            d = reldiff(a, b); ok = d < 1f-3
            @printf("  [%s] %-12s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
            ok || (nfail += 1)
        end
        zc = Array(H.S.sh.zwt_col); zd = Array(Sd.sh.zwt_col)
        nb = Array(H.S.col_nbedrock)
        for c in 1:H.nc
            @printf("    col %d: nbedrock = %2d  zwt_cpu = %9.5f  zwt_dev = %9.5f\n",
                    c, nb[c], zc[c], zd[c])
        end
        println()
    end

    println(nfail == 0 ? "  BOTH water-table variants MATCH CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
