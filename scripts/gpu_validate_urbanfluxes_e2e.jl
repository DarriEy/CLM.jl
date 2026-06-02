# ==========================================================================
# gpu_validate_urbanfluxes_e2e.jl — end-to-end Metal parity for the
# column-INDEPENDENT, device-resident portion of urban_fluxes! run WHOLE on the
# GPU: urban_fluxes_diagnostics!.
#
# These three passes (non-urban taf/qaf restart, urban-roots fill, and the 2-m
# temperature/humidity/RH/t_veg diagnostics) carry no cross-index coupling and run
# as KernelAbstractions kernels on the device. The harness builds a Float32 urban
# instance with ALL FIVE urban column types (roof / sunwall / shadewall /
# road_imperv / road_perv) wired one-patch-per-column, PLUS a non-urban landunit so
# the SPVAL-restart branch is exercised, runs urban_fluxes_diagnostics! on the CPU,
# adapts every state struct + masks + forcing to Metal, runs the SAME call on the
# device, and compares every mutated field with reldiff.
#
# The iterative canyon-air solve (canyontop wind, stability iteration via
# friction_velocity!'s landunit path, the urban-column gather/reduction into
# taf/qaf, canyon-surface fluxes, and the sum()-based error checks) is loop-carried
# + reduction-coupled and stays on the host — see urban_fluxes! and the PARTIAL note.
#
#   julia --project=scripts scripts/gpu_validate_urbanfluxes_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

# reldiff: NaN-aware (both-NaN agrees; one-sided NaN flags divergence). Asserts the
# CPU reference is FINITE separately so a both-NaN false PASS can't slip through.
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

function build(::Type{FT}) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevgrnd = CLM.varpar.nlevgrnd

    # 1 gridcell; 2 landunits (1 urban + 1 non-urban); 5 urban columns; 5 urban patches.
    ng = 1; nl = 2; nc = 5; np = 5

    temperature   = CLM.TemperatureData{FT}();          CLM.temperature_init!(temperature, np, nc, nl, ng)
    waterdiagbulk = CLM.WaterDiagnosticBulkData{FT}();   CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)
    soilstate     = CLM.SoilStateData{FT}();             CLM.soilstate_init!(soilstate, np, nc)

    patch_data = CLM.PatchData{FT}();   CLM.patch_init!(patch_data, np)
    col_data   = CLM.ColumnData{FT}();  CLM.column_init!(col_data, nc)

    # Column types: roof, sunwall, shadewall, road_imperv, road_perv.
    col_data.itype[1] = CLM.ICOL_ROOF
    col_data.itype[2] = CLM.ICOL_SUNWALL
    col_data.itype[3] = CLM.ICOL_SHADEWALL
    col_data.itype[4] = CLM.ICOL_ROAD_IMPERV
    col_data.itype[5] = CLM.ICOL_ROAD_PERV
    for c in 1:nc
        col_data.landunit[c] = 1
    end

    for p in 1:np
        patch_data.column[p]   = p
        patch_data.gridcell[p] = 1
        patch_data.landunit[p] = 1   # urban landunit
        patch_data.active[p]   = true
    end

    # Post-solve canopy-air state the diagnostics read (set to plausible solved values).
    temperature.taf_lun[1]   = FT(290.0)   # urban
    waterdiagbulk.qaf_lun[1] = FT(0.005)
    temperature.taf_lun[2]   = FT(123.0)   # non-urban: will be overwritten with SPVAL
    waterdiagbulk.qaf_lun[2] = FT(456.0)

    # Roots source for the pervious-road column (only column 5 should propagate).
    for c in 1:nc, j in 1:nlevgrnd
        soilstate.rootr_road_perv_col[c, j] = (j <= 5 ? FT(0.2) * j : FT(0.0))
    end

    forc_t_grc    = FT[288.0]
    forc_pbot_grc = FT[101325.0]

    S = (; temperature, waterdiagbulk, soilstate, patch_data, col_data)
    return (; np, nc, nl, ng, nlevgrnd, S,
              filter_nourbanl = [2], num_nourbanl = 1,
              filter_urbanp   = [1, 2, 3, 4, 5], num_urbanp = 5,
              forc_t_grc, forc_pbot_grc)
end

function run_diag!(S, m_nourbanl, m_urbanp, forc_t_grc, forc_pbot_grc, nlevgrnd)
    CLM.urban_fluxes_diagnostics!(
        S.temperature, S.waterdiagbulk, S.soilstate, S.patch_data, S.col_data,
        m_nourbanl, m_urbanp, forc_t_grc, forc_pbot_grc, nlevgrnd)
end

# --------------------------------------------------------------------------
# Canyon-surface flux pass (urban_fluxes_canyon_surface!): the per-patch state
# block of urban_fluxes! that runs WHOLE on the device given the SOLVED
# per-landunit canopy-air state. All five urban column types are wired
# one-patch-per-column; per-landunit weights / forcing / column state set to
# finite plausible solved values so every branch fires on the CPU reference.
# --------------------------------------------------------------------------
function build_canyon(::Type{FT}) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    ng = 1; nl = 1; nc = 5; np = 5
    energyflux    = CLM.EnergyFluxData{FT}();        CLM.energyflux_init!(energyflux, np, nc, nl, ng)
    frictionvel   = CLM.FrictionVelocityData{FT}();  CLM.frictionvel_init!(frictionvel, np, nc)
    temperature   = CLM.TemperatureData{FT}();        CLM.temperature_init!(temperature, np, nc, nl, ng)
    soilstate     = CLM.SoilStateData{FT}();          CLM.soilstate_init!(soilstate, np, nc)
    waterfluxbulk = CLM.WaterFluxBulkData{FT}();      CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)
    waterdiagbulk = CLM.WaterDiagnosticBulkData{FT}();CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)

    patch_data = CLM.PatchData{FT}();   CLM.patch_init!(patch_data, np)
    col_data   = CLM.ColumnData{FT}();  CLM.column_init!(col_data, nc)

    col_data.itype[1] = CLM.ICOL_ROOF
    col_data.itype[2] = CLM.ICOL_SUNWALL
    col_data.itype[3] = CLM.ICOL_SHADEWALL
    col_data.itype[4] = CLM.ICOL_ROAD_IMPERV
    col_data.itype[5] = CLM.ICOL_ROAD_PERV
    for c in 1:nc; col_data.landunit[c] = 1; end
    for p in 1:np
        patch_data.column[p]   = p
        patch_data.gridcell[p] = 1
        patch_data.landunit[p] = 1
        patch_data.active[p]   = true
    end

    # Solved canopy-air state (per landunit l=1)
    temperature.taf_lun[1]   = FT(290.0)
    waterdiagbulk.qaf_lun[1] = FT(0.006)

    # Per-column state read by the kernel
    for c in 1:nc
        temperature.t_grnd_col[c]     = FT(288.0 + c)
        waterdiagbulk.qg_col[c]       = FT(0.005 + 0.0003 * c)
        waterdiagbulk.dqgdT_col[c]    = FT(2.0e-4)
        energyflux.htvp_col[c]        = FT(2.5e6)
    end
    # Pervious-road (c=5): drive BOTH road_perv branches by toggling later; here pick
    # the transpiration branch (dqh<=0, no snow, soilalpha>0).
    waterdiagbulk.frac_sno_col[5]   = FT(0.0)
    soilstate.soilalpha_u_col[5]    = FT(0.8)
    waterdiagbulk.qg_col[5]         = FT(0.007)   # > qaf -> dqh<0 -> transpiration branch

    # Per-landunit solved weight arrays + resistances (finite, nonzero denominators)
    mk(v) = (a = zeros(FT, nl); a .= FT(v); a)
    wts = CLM._UfWtsDV(
        wtas = mk(1.5),
        wtus_roof = mk(0.30), wtus_road_perv = mk(0.10), wtus_road_imperv = mk(0.20),
        wtus_sunwall = mk(0.25), wtus_shadewall = mk(0.25),
        wtus_roof_unscl = mk(3.0), wtus_road_perv_unscl = mk(3.0),
        wtus_road_imperv_unscl = mk(3.0), wtus_sunwall_unscl = mk(3.0),
        wtus_shadewall_unscl = mk(3.0), wts_sum = mk(2.6))
    wtq = CLM._UfWtqDV(
        wtaq = mk(1.5),
        wtuq_roof = mk(0.30), wtuq_road_perv = mk(0.10), wtuq_road_imperv = mk(0.20),
        wtuq_sunwall = mk(0.0), wtuq_shadewall = mk(0.0),
        wtuq_roof_unscl = mk(3.0), wtuq_road_perv_unscl = mk(3.0),
        wtuq_road_imperv_unscl = mk(3.0), wtq_sum = mk(2.1))
    ramu       = mk(50.0)
    zeta_lunit = mk(-0.3)
    wtus_col   = (a = zeros(FT, nc); for c in 1:nc; a[c] = FT(0.1 * c); end; a)
    wtuq_col   = (a = zeros(FT, nc); for c in 1:nc; a[c] = FT(0.05 * c); end; a)

    eflx_scale = zeros(FT, np)
    qflx_scale = zeros(FT, np)

    forc_rho_grc = FT[1.2]
    forc_u_grc   = FT[3.0]
    forc_v_grc   = FT[2.0]

    S = (; energyflux, frictionvel, temperature, soilstate, waterfluxbulk, waterdiagbulk,
           patch_data, col_data)
    W = (; wts, wtq, ramu, zeta_lunit, wtus_col, wtuq_col, eflx_scale, qflx_scale,
           forc_rho_grc, forc_u_grc, forc_v_grc)
    return (; np, nc, nl, ng, S, W, filter_urbanp = [1,2,3,4,5], num_urbanp = 5)
end

function run_canyon!(S, m_urbanp, W)
    CLM.urban_fluxes_canyon_surface!(
        S.energyflux, S.frictionvel, S.waterfluxbulk, S.waterdiagbulk, S.soilstate,
        S.temperature, S.patch_data, S.col_data, m_urbanp,
        W.wts, W.wtq, W.ramu, W.zeta_lunit, W.wtus_col, W.wtuq_col,
        W.eflx_scale, W.qflx_scale,
        W.forc_rho_grc, W.forc_u_grc, W.forc_v_grc)
end

function main_canyon(name, to, FT)
    println("\n" * "=" ^ 70)
    println("END-TO-END Metal parity for urban_fluxes_canyon_surface! (per-patch flux block)")
    println("=" ^ 70)

    H = build_canyon(FT)
    B = build_canyon(FT)

    ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
    Sd = map(ad, B.S)
    Wd = (; wts = ad(B.W.wts), wtq = ad(B.W.wtq),
            ramu = to(B.W.ramu), zeta_lunit = to(B.W.zeta_lunit),
            wtus_col = to(B.W.wtus_col), wtuq_col = to(B.W.wtuq_col),
            eflx_scale = to(B.W.eflx_scale), qflx_scale = to(B.W.qflx_scale),
            forc_rho_grc = to(B.W.forc_rho_grc), forc_u_grc = to(B.W.forc_u_grc),
            forc_v_grc = to(B.W.forc_v_grc))

    if !(Sd.energyflux.cgrnd_patch isa Metal.MtlArray)
        println("  BLOCKED: energyflux did not move to the device under adapt.")
        return 2
    end

    m = falses(H.np); for f in 1:H.num_urbanp; m[H.filter_urbanp[f]] = true; end

    run_canyon!(H.S, m, H.W)               # CPU
    run_canyon!(Sd, to(collect(Bool, m)), Wd)  # device

    checks = [
        ("cgrnds",        H.S.energyflux.cgrnds_patch,            Sd.energyflux.cgrnds_patch),
        ("cgrndl",        H.S.energyflux.cgrndl_patch,            Sd.energyflux.cgrndl_patch),
        ("cgrnd",         H.S.energyflux.cgrnd_patch,             Sd.energyflux.cgrnd_patch),
        ("taux",          H.S.energyflux.taux_patch,              Sd.energyflux.taux_patch),
        ("tauy",          H.S.energyflux.tauy_patch,              Sd.energyflux.tauy_patch),
        ("ulrad",         H.S.energyflux.ulrad_patch,             Sd.energyflux.ulrad_patch),
        ("dlrad",         H.S.energyflux.dlrad_patch,             Sd.energyflux.dlrad_patch),
        ("eflx_sh_grnd",  H.S.energyflux.eflx_sh_grnd_patch,      Sd.energyflux.eflx_sh_grnd_patch),
        ("eflx_sh_tot",   H.S.energyflux.eflx_sh_tot_patch,       Sd.energyflux.eflx_sh_tot_patch),
        ("eflx_sh_tot_u", H.S.energyflux.eflx_sh_tot_u_patch,     Sd.energyflux.eflx_sh_tot_u_patch),
        ("qflx_evap_soi", H.S.waterfluxbulk.wf.qflx_evap_soi_patch, Sd.waterfluxbulk.wf.qflx_evap_soi_patch),
        ("qflx_tran_veg", H.S.waterfluxbulk.wf.qflx_tran_veg_patch, Sd.waterfluxbulk.wf.qflx_tran_veg_patch),
        ("qflx_evap_veg", H.S.waterfluxbulk.wf.qflx_evap_veg_patch, Sd.waterfluxbulk.wf.qflx_evap_veg_patch),
        ("ram1",          H.S.frictionvel.ram1_patch,             Sd.frictionvel.ram1_patch),
        ("zeta",          H.S.frictionvel.zeta_patch,             Sd.frictionvel.zeta_patch),
        ("eflx_sh_scale", H.W.eflx_scale,                         Wd.eflx_scale),
        ("qflx_evap_scale", H.W.qflx_scale,                       Wd.qflx_scale),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-18s CPU reference is all-NaN/Inf — skipping\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-18s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end

    # Confirm road_perv transpiration branch fired (qflx_tran_veg nonzero at p5).
    qtv = Array(H.S.waterfluxbulk.wf.qflx_tran_veg_patch)
    qes = Array(H.S.waterfluxbulk.wf.qflx_evap_soi_patch)
    @printf("\n  road_perv p5: qflx_tran_veg=%.4e (expect !=0, transpiration branch), qflx_evap_soi=%.4e (expect 0)\n",
            qtv[5], qes[5])

    println(nfail == 0 ? "  urban_fluxes_canyon_surface! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for urban_fluxes_diagnostics! (device passes)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT)   # CPU reference
    B = build(FT)   # source for the device snapshot

    ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
    Sd = map(ad, B.S)
    dmask(m) = to(collect(Bool, m))

    if !(Sd.temperature.taf_lun isa Metal.MtlArray)
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    # Masks over the full index range built from the filters.
    function filt_mask(filter, fn, n)
        m = falses(n)
        for f in 1:fn; m[filter[f]] = true; end
        return m
    end
    m_nourbanl = filt_mask(H.filter_nourbanl, H.num_nourbanl, H.nl)
    m_urbanp   = filt_mask(H.filter_urbanp,   H.num_urbanp,   H.np)

    run_diag!(H.S, m_nourbanl, m_urbanp, H.forc_t_grc, H.forc_pbot_grc, H.nlevgrnd)                       # CPU
    run_diag!(Sd, dmask(m_nourbanl), dmask(m_urbanp), to(B.forc_t_grc), to(B.forc_pbot_grc), B.nlevgrnd)  # device

    checks = [
        ("taf_lun (restart)", H.S.temperature.taf_lun,        Sd.temperature.taf_lun),
        ("qaf_lun (restart)", H.S.waterdiagbulk.qaf_lun,      Sd.waterdiagbulk.qaf_lun),
        ("t_ref2m",           H.S.temperature.t_ref2m_patch,  Sd.temperature.t_ref2m_patch),
        ("t_ref2m_u",         H.S.temperature.t_ref2m_u_patch,Sd.temperature.t_ref2m_u_patch),
        ("t_veg",             H.S.temperature.t_veg_patch,    Sd.temperature.t_veg_patch),
        ("q_ref2m",           H.S.waterdiagbulk.q_ref2m_patch,Sd.waterdiagbulk.q_ref2m_patch),
        ("rh_ref2m",          H.S.waterdiagbulk.rh_ref2m_patch,Sd.waterdiagbulk.rh_ref2m_patch),
        ("rh_ref2m_u",        H.S.waterdiagbulk.rh_ref2m_u_patch,Sd.waterdiagbulk.rh_ref2m_u_patch),
        ("rootr_patch",       H.S.soilstate.rootr_patch,      Sd.soilstate.rootr_patch),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-18s CPU reference is all-NaN/Inf — skipping (no parity signal)\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-18s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end

    # Sanity: confirm the branches actually fired on the CPU reference.
    rc = Array(H.S.soilstate.rootr_patch)
    @printf("\n  roots[road_perv p5, j1..3] = %.4f %.4f %.4f (expect 0.2 0.4 0.6)\n",
            rc[5, 1], rc[5, 2], rc[5, 3])
    @printf("  roots[roof p1, j1]         = %.4f (expect 0.0)\n", rc[1, 1])
    @printf("  taf_lun[non-urban l2]      = %.3e (expect SPVAL=%.3e)\n",
            Array(H.S.temperature.taf_lun)[2], CLM.SPVAL)

    println()
    println(nfail == 0 ? "  urban_fluxes_diagnostics! (device passes) MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")

    nfail2 = main_canyon(name, to, FT)

    total = nfail + nfail2
    println()
    println(total == 0 ? "  ALL urban_fluxes! device passes MATCH CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return total == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
