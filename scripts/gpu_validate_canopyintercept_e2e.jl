# ==========================================================================
# gpu_validate_canopyintercept_e2e.jl — end-to-end GPU parity for the WHOLE
# canopy_interception_and_throughfall! driver (all 9 steps: top-of-canopy
# inputs, interception/throughfall, canopy-excess, snow unloading, the
# patch->column scatter of fluxes onto the ground via _scatter_add!, and the
# wet/dry fraction diagnostic).
#
# There is no whole-function test in the suite (test_canopy_hydrology.jl only
# covers the sub-pieces), so this builds a small Float32 instance directly: two
# vegetated patches sharing ONE column (so the patch->column weighted scatter is
# actually exercised), runs the function on the CPU, adapts every state struct
# (and the BitVector masks -> device Vector{Bool}, the forcing/irrig arrays) to
# Metal, runs the SAME call on the device, and compares the mutated outputs.
#
#   julia --project=scripts scripts/gpu_validate_canopyintercept_e2e.jl
#
# NOTE (banked A1 lesson): reldiff PASSES when both sides are NaN. The setup
# below sets every input/param to a real finite value, and main() asserts the
# CPU reference fields are finite before trusting the parity numbers.
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end

allfinite(a) = all(isfinite, Array(a))

# Set the module-global params/control to real CLM defaults (host-resident; the
# wrappers resolve them to scalars before each launch).
function setconfig!()
    CLM.canopy_hydrology_read_nml!(use_clm5_fpi = false)
    CLM.canopy_hydrology_read_params!()   # CLM defaults (0.04, 6.0, 1.87e5, 1.56e5, 1.0, 1.0)
end

# Two vegetated patches on one column; np=2, nc=1, ng=1, nl=1.
function build(::Type{FT}) where {FT}
    np = 2; nc = 1; nl = 1; ng = 1
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    patch = CLM.PatchData{FT}(); CLM.patch_init!(patch, np)
    patch.column   .= 1
    patch.gridcell .= 1
    patch.landunit .= 1
    patch.itype    .= 1
    patch.active   .= true
    patch.wtcol[1] = FT(0.6); patch.wtcol[2] = FT(0.4)   # weights sum to 1 on the shared column

    col = CLM.ColumnData{FT}(); CLM.column_init!(col, nc)
    col.itype[1] = 1                                     # soil column (not sun/shade wall)

    cs = CLM.CanopyStateData{FT}(); CLM.canopystate_init!(cs, np)
    cs.frac_veg_nosno_patch .= 1
    cs.elai_patch[1] = FT(2.0); cs.elai_patch[2] = FT(1.5)
    cs.esai_patch[1] = FT(0.5); cs.esai_patch[2] = FT(0.3)

    water = CLM.WaterData{FT}()
    CLM.water_init_for_testing!(water, nc, np, nl, ng)
    # Seed canopy storage so canopy-excess / snow-unloading / frac-wet branches are live.
    water.waterstatebulk_inst.ws.snocan_patch[1] = FT(0.5)
    water.waterstatebulk_inst.ws.snocan_patch[2] = FT(0.2)
    water.waterstatebulk_inst.ws.liqcan_patch[1] = FT(0.3)
    water.waterstatebulk_inst.ws.liqcan_patch[2] = FT(0.1)

    masks = (; soilp  = BitVector([true, true]),
               nolakep = BitVector([true, true]),
               nolakec = BitVector([true]))
    forc = (; rain  = fill(FT(3.0e-4), nc),
              snow  = fill(FT(1.0e-4), nc),
              t     = fill(FT(272.0), nc),     # column-level air T (above 270.15 -> temp unload active)
              wind  = fill(FT(5.0), ng))       # gridcell-level wind
    irrig = (; sprinkler = fill(FT(0.5e-4), np),
               drip      = fill(FT(0.2e-4), np))

    return (; np, nc, ng, nl, patch, col, cs, water, masks, forc, irrig,
              dtime = FT(1800.0))
end

run_cit!(S, masks, forc, irrig, dtime, np, nc, ng) =
    CLM.canopy_interception_and_throughfall!(
        S.patch, S.col, S.cs, S.water, dtime,
        masks.soilp, masks.nolakep, masks.nolakec,
        1:np, 1:nc, 1:ng,
        forc.rain, forc.snow, forc.t, forc.wind,
        irrig.sprinkler, irrig.drip)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for canopy_interception_and_throughfall! (whole driver)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    setconfig!()
    H = build(FT)             # CPU reference
    setconfig!()
    B = build(FT)             # to be adapted to device (fresh initial state)

    ad(x) = CLM.Adapt.adapt(device_array_type(), x)
    Sd = (; patch = ad(B.patch), col = ad(B.col), cs = ad(B.cs), water = ad(B.water))
    forc_d  = (; (k => ad(getfield(B.forc,  k)) for k in keys(B.forc))...)
    irrig_d = (; (k => ad(getfield(B.irrig, k)) for k in keys(B.irrig))...)
    # BitVector masks -> device Vector{Bool}
    masks_d = (; soilp   = to(collect(Bool, B.masks.soilp)),
                 nolakep = to(collect(Bool, B.masks.nolakep)),
                 nolakec = to(collect(Bool, B.masks.nolakec)))

    if !(Sd.water.waterstatebulk_inst.ws.snocan_patch isa device_array_type())
        println("  BLOCKED: WaterData did not move to the device under adapt.")
        return 2
    end

    run_cit!(H, H.masks, H.forc, H.irrig, H.dtime, H.np, H.nc, H.ng)        # CPU
    run_cit!(Sd, masks_d, forc_d, irrig_d, B.dtime, B.np, B.nc, B.ng)       # device

    Hwf = H.water.waterfluxbulk_inst; Hws = H.water.waterstatebulk_inst; Hwd = H.water.waterdiagnosticbulk_inst
    Dwf = Sd.water.waterfluxbulk_inst; Dws = Sd.water.waterstatebulk_inst; Dwd = Sd.water.waterdiagnosticbulk_inst

    # FIRST verify the CPU reference is finite (reldiff PASSES on NaN==NaN).
    cpu_fields = [
        ("qflx_through_snow", Hwf.wf.qflx_through_snow_patch),
        ("qflx_through_liq",  Hwf.wf.qflx_through_liq_patch),
        ("qflx_intercept_sn", Hwf.wf.qflx_intercepted_snow_patch),
        ("qflx_intercept_lq", Hwf.wf.qflx_intercepted_liq_patch),
        ("qflx_snocanfall",   Hwf.wf.qflx_snocanfall_patch),
        ("qflx_liqcanfall",   Hwf.wf.qflx_liqcanfall_patch),
        ("qflx_snow_unload",  Hwf.wf.qflx_snow_unload_patch),
        ("qflx_snow_grnd_col",Hwf.wf.qflx_snow_grnd_col),
        ("qflx_liq_grnd_col", Hwf.wf.qflx_liq_grnd_col),
        ("qflx_snow_h2osfc",  Hwf.wf.qflx_snow_h2osfc_col),
        ("snocan",            Hws.ws.snocan_patch),
        ("liqcan",            Hws.ws.liqcan_patch),
        ("fwet",              Hwd.fwet_patch),
        ("fdry",              Hwd.fdry_patch),
        ("fcansno",           Hwd.fcansno_patch),
    ]
    nbad = 0
    for (nm, a) in cpu_fields
        if !allfinite(a)
            @printf("  [NONFINITE] CPU reference %-18s = %s\n", nm, string(Array(a)))
            nbad += 1
        end
    end
    if nbad > 0
        println("\n  BLOCKED: CPU reference has non-finite fields — parity numbers untrustworthy.")
        return 2
    end

    checks = [
        ("qflx_through_snow", Hwf.wf.qflx_through_snow_patch,     Dwf.wf.qflx_through_snow_patch),
        ("qflx_through_liq",  Hwf.wf.qflx_through_liq_patch,      Dwf.wf.qflx_through_liq_patch),
        ("qflx_intercept_sn", Hwf.wf.qflx_intercepted_snow_patch, Dwf.wf.qflx_intercepted_snow_patch),
        ("qflx_intercept_lq", Hwf.wf.qflx_intercepted_liq_patch,  Dwf.wf.qflx_intercepted_liq_patch),
        ("qflx_snocanfall",   Hwf.wf.qflx_snocanfall_patch,       Dwf.wf.qflx_snocanfall_patch),
        ("qflx_liqcanfall",   Hwf.wf.qflx_liqcanfall_patch,       Dwf.wf.qflx_liqcanfall_patch),
        ("qflx_snow_unload",  Hwf.wf.qflx_snow_unload_patch,      Dwf.wf.qflx_snow_unload_patch),
        ("qflx_snow_grnd_col",Hwf.wf.qflx_snow_grnd_col,          Dwf.wf.qflx_snow_grnd_col),
        ("qflx_liq_grnd_col", Hwf.wf.qflx_liq_grnd_col,           Dwf.wf.qflx_liq_grnd_col),
        ("qflx_snow_h2osfc",  Hwf.wf.qflx_snow_h2osfc_col,        Dwf.wf.qflx_snow_h2osfc_col),
        ("snocan",            Hws.ws.snocan_patch,                Dws.ws.snocan_patch),
        ("liqcan",            Hws.ws.liqcan_patch,                Dws.ws.liqcan_patch),
        ("fwet",              Hwd.fwet_patch,                     Dwd.fwet_patch),
        ("fdry",              Hwd.fdry_patch,                     Dwd.fdry_patch),
        ("fcansno",           Hwd.fcansno_patch,                  Dwd.fcansno_patch),
    ]
    nfail = 0
    for (nm, a, b) in checks
        d = reldiff(a, b); ok = d < 1f-2
        @printf("  [%s] %-18s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ?
        "  WHOLE canopy_interception_and_throughfall! MATCHES CPU ON $name ($FT) ✓" :
        "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
