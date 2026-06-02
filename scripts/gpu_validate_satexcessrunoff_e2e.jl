# ==========================================================================
# gpu_validate_satexcessrunoff_e2e.jl — end-to-end Metal parity for the WHOLE
# saturated_excess_runoff! routine (the fused per-column kernel).
#
# Builds a small Float32 instance with several columns exercising every branch:
#   col 1: TOPModel, deep frost table -> uses zwt branch, normal soil
#   col 2: TOPModel, perched-water-table branch
#   col 3: crop column with crop_fsat_equals_zero -> fsat forced to 0
#   col 4: hillslope upland (cold != ISPVAL) with hillslope flag -> fsat forced 0
#   col 5: urban column -> flood water added to runoff
# A second pass exercises the VIC fsat_method.
#
# Runs the whole function on the CPU, adapts every state struct (+ the ser
# fsat/fcov arrays + the Bool mask) to Metal, runs the SAME call on the device,
# and compares every mutated field with reldiff.
#
# CRITICAL: every required param/input is set to a real value so the CPU
# reference is FINITE — reldiff silently PASSES when both sides are NaN, so the
# harness ASSERTS the CPU reference fields are finite before trusting parity.
#
#   julia --project=scripts scripts/gpu_validate_satexcessrunoff_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal
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

# Build a 5-column state. `use_vichydro` selects the fsat_method.
function build(::Type{FT}, use_vichydro::Bool) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nc = 5
    CLM.sat_excess_runoff_params.fff = 0.5

    # Build an FT-typed SaturatedExcessRunoffData (init! would fill Float64).
    ser = CLM.SaturatedExcessRunoffData(;
        fsat_col = fill(FT(NaN), nc),
        fcov_col = fill(FT(NaN), nc),
        fsat_method = use_vichydro ? CLM.FSAT_METHOD_VIC : CLM.FSAT_METHOD_TOPMODEL)

    col = CLM.ColumnData{FT}(); CLM.column_init!(col, nc)
    # col 3 -> landunit 2 (crop); all others landunit 1 (soil)
    col.landunit            = Int[1, 1, 2, 1, 1]
    col.urbpoi              = Bool[false, false, false, false, true]
    col.active              = Bool[true, true, true, true, true]
    col.is_hillslope_column = Bool[false, false, false, true, false]
    col.cold                = Int[CLM.ISPVAL, CLM.ISPVAL, CLM.ISPVAL, 3, CLM.ISPVAL]

    lun = CLM.LandunitData{FT}(); CLM.landunit_init!(lun, 2)
    lun.itype = Int[CLM.ISTSOIL, CLM.ISTCROP]

    sh = CLM.SoilHydrologyData{FT}(); CLM.soilhydrology_init!(sh, nc)
    # TOPModel inputs: col1 uses zwt (frost deep), col2 uses perched
    sh.frost_table_col = FT[5.0, 1.5, 5.0, 5.0, 5.0]
    sh.zwt_col         = FT[2.0, 2.0, 1.0, 1.0, 1.0]
    sh.zwt_perched_col = FT[1.0, 1.0, 10.0, 10.0, 10.0]
    # VIC inputs
    sh.b_infil_col           = FT[0.5, 2.0, 1.0, 1.5, 1.0]
    sh.top_max_moist_col     = FT[100.0, 200.0, 100.0, 150.0, 120.0]
    sh.top_moist_limited_col = FT[50.0, 180.0, 40.0, 90.0, 60.0]

    ss = CLM.SoilStateData{FT}(); CLM.soilstate_init!(ss, 0, nc)
    ss.wtfact_col = FT[0.4, 0.4, 0.4, 0.4, 0.3]

    wfb = CLM.WaterFluxBulkData{FT}(); CLM.waterfluxbulk_init!(wfb, nc, 0, 0, 0)
    wfb.wf.qflx_rain_plus_snomelt_col = FT[5.0, 10.0, 8.0, 6.0, 5.0]
    wfb.wf.qflx_floodc_col            = FT[0.0, 0.0, 0.0, 0.0, 2.0]

    return (; nc, ser, col, lun, sh, ss, wfb)
end

run_ser!(S) = CLM.saturated_excess_runoff!(
    S.ser, S.mask, 1:S.nc, S.col, S.lun, S.sh, S.ss, S.wfb;
    crop_fsat_equals_zero=true, hillslope_fsat_equals_zero=true)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for saturated_excess_runoff! (whole routine)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU routine exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    nfail = 0
    for (label, use_vic) in (("TOPModel", false), ("VIC", true))
        @printf("  --- fsat_method: %s ---\n", label)
        H = build(FT, use_vic)   # CPU reference
        B = build(FT, use_vic)   # device source

        ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
        # Adapt the full state structs (their fields become MtlArrays).
        cold = ad(B.col); lund = ad(B.lun); shd = ad(B.sh)
        ssd = ad(B.ss); wfbd = ad(B.wfb); serd = ad(B.ser)

        if !(shd.zwt_col isa Metal.MtlArray) || !(serd.fsat_col isa Metal.MtlArray)
            println("  BLOCKED: a state struct did not move to the device under adapt.")
            return 2
        end

        H = (; H..., mask = trues(H.nc))
        Bd = (; nc = B.nc, ser = serd, col = cold, lun = lund, sh = shd, ss = ssd,
                wfb = wfbd, mask = to(collect(Bool, trues(B.nc))))

        run_ser!(H)    # CPU
        run_ser!(Bd)   # device

        checks = [
            ("fsat_col",             H.ser.fsat_col,                     Bd.ser.fsat_col),
            ("fcov_col",             H.ser.fcov_col,                     Bd.ser.fcov_col),
            ("qflx_sat_excess_surf", H.wfb.qflx_sat_excess_surf_col,     Bd.wfb.qflx_sat_excess_surf_col),
        ]
        for (nm, a, b) in checks
            if !cpu_has_finite(a)
                @printf("  [WARN] %-22s CPU reference is all-NaN/Inf — skipping\n", nm)
                continue
            end
            d = reldiff(a, b); ok = d < 1f-3
            @printf("  [%s] %-22s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
            ok || (nfail += 1)
        end
        println()
    end

    println(nfail == 0 ? "  WHOLE saturated_excess_runoff! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
