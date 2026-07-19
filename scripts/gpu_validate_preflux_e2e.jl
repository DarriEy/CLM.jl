# ==========================================================================
# gpu_validate_preflux_e2e.jl — end-to-end GPU parity for the kernelized
# portions of biogeophys_pre_flux_calcs!: set_z0m_displa! (step 1) and
# calc_initial_temp_energy! (step 3). Both run WHOLE on Metal here.
#
# NOTE on scope: the orchestrator's step 2, set_roughness_and_forc_heights_nonlake!,
# lives in src/types/friction_velocity.jl (a NON-assigned file owned elsewhere) and
# is still a host loop, so calling the full biogeophys_pre_flux_calcs! on device
# would break inside that peer function (scalar indexing of MtlArrays). This harness
# therefore validates the two steps that pre_flux_calcs.jl owns and kernelizes —
# every @kernel in pre_flux_calcs.jl is exercised here: z0m/displa, tssbef, t_grnd,
# emg, htvp, thv/beta/zii, zero-fluxes, emv, thm.
#
#   julia --project=scripts scripts/gpu_validate_preflux_e2e.jl
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
cpu_has_finite(a) = any(isfinite, Array(a))

function build(::Type{FT}) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    np = 1; nc = 1; nl = 1; ng = 1
    nlevsno  = CLM.varpar.nlevsno
    nlev_soisno = nlevsno + CLM.varpar.nlevmaxurbgrnd

    canopystate = CLM.CanopyStateData{FT}(); CLM.canopystate_init!(canopystate, np)
    frictionvel = CLM.FrictionVelocityData{FT}(); CLM.frictionvel_init!(frictionvel, np, nc)
    temperature = CLM.TemperatureData{FT}(); CLM.temperature_init!(temperature, np, nc, nl, ng)
    energyflux = CLM.EnergyFluxData{FT}(); CLM.energyflux_init!(energyflux, np, nc, nl, ng)
    waterstatebulk = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)
    waterdiagbulk = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)

    patch_data = CLM.PatchData{FT}(); CLM.patch_init!(patch_data, np)
    patch_data.column[1] = 1; patch_data.gridcell[1] = 1; patch_data.landunit[1] = 1
    patch_data.itype[1] = 1

    col_data = CLM.ColumnData{FT}(); CLM.column_init!(col_data, nc)
    col_data.snl[1] = -1                 # one snow layer -> exercise the snow branch
    col_data.landunit[1] = 1

    lun_data = CLM.LandunitData{FT}(); CLM.landunit_init!(lun_data, nl)
    lun_data.itype[1] = CLM.ISTSOIL; lun_data.urbpoi[1] = false

    # canopy / vegetation
    canopystate.htop_patch[1] = FT(15.0)
    canopystate.elai_patch[1] = FT(2.0)
    canopystate.esai_patch[1] = FT(0.5)

    # temperatures (combined snow+soil indexing): set snow + soil layers
    for j in 1:nlev_soisno
        temperature.t_soisno_col[1, j] = FT(278.0)
        temperature.t_ssbef_col[1, j]  = FT(0.0)
    end
    temperature.t_soisno_col[1, nlevsno] = FT(272.0)        # top snow layer (snl=-1 -> jtop=nlevsno)
    temperature.t_h2osfc_col[1] = FT(275.0)

    # water states (combined indexing). Top layer (jtop) ice>0,liq=0 -> sublimation.
    jtop = col_data.snl[1] + 1 + nlevsno
    for j in 1:nlev_soisno
        waterstatebulk.ws.h2osoi_liq_col[1, j] = FT(10.0)
        waterstatebulk.ws.h2osoi_ice_col[1, j] = FT(0.0)
    end
    waterstatebulk.ws.h2osoi_liq_col[1, jtop] = FT(0.0)
    waterstatebulk.ws.h2osoi_ice_col[1, jtop] = FT(5.0)

    waterdiagbulk.frac_sno_col[1]     = FT(0.4)
    waterdiagbulk.frac_sno_eff_col[1] = FT(0.4)
    waterdiagbulk.frac_h2osfc_col[1]  = FT(0.1)

    frictionvel.forc_hgt_t_patch[1] = FT(10.0)

    forc_t  = FT[280.0]
    forc_th = FT[281.0]
    forc_q  = FT[0.006]

    # PFT roughness/displacement ratios (avoid the fallback so the branch is real)
    npft = 14
    z0mr   = fill(FT(0.055), npft)
    displar = fill(FT(0.67), npft)

    S = (; canopystate, frictionvel, temperature, energyflux,
           waterstatebulk, waterdiagbulk, patch_data, col_data, lun_data)
    return (; np, nc, nl, ng, S, forc_t, forc_th, forc_q, z0mr, displar)
end

function run_preflux!(S, forc_t, forc_th, forc_q, z0mr, displar, m_c, m_p)
    # Step 1: z0m / displa (whole function, on device)
    CLM.set_z0m_displa!(S.canopystate, S.frictionvel, S.patch_data, S.lun_data,
                        m_p, 1:length(m_p);
                        z0param_method="ZengWang2007", z0mr=z0mr, displar=displar)
    # Step 3: initialize temperature/energy variables (whole function, on device)
    CLM.calc_initial_temp_energy!(S.temperature, S.energyflux, S.canopystate,
                                  S.frictionvel, S.waterstatebulk, S.waterdiagbulk,
                                  S.col_data, S.lun_data, S.patch_data,
                                  forc_t, forc_th, forc_q,
                                  m_c, m_p, 1:length(m_c), 1:length(m_p))
    return nothing
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for pre_flux_calcs.jl kernels")
    println("(set_z0m_displa! + calc_initial_temp_energy!, whole on device)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT)
    B = build(FT)

    ad(x) = CLM.Adapt.adapt(device_array_type(), x)
    Sd = map(ad, B.S)
    dmask(m) = to(collect(Bool, m))

    if !(Sd.temperature.t_grnd_col isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    m_c = trues(H.nc); m_p = trues(H.np)

    run_preflux!(H.S, H.forc_t, H.forc_th, H.forc_q, H.z0mr, H.displar, m_c, m_p)
    run_preflux!(Sd, to(B.forc_t), to(B.forc_th), to(B.forc_q),
                 to(B.z0mr), to(B.displar), dmask(m_c), dmask(m_p))

    checks = [
        ("z0m_patch",   H.S.canopystate.z0m_patch,        Sd.canopystate.z0m_patch),
        ("displa_patch",H.S.canopystate.displa_patch,     Sd.canopystate.displa_patch),
        ("t_ssbef",     H.S.temperature.t_ssbef_col,      Sd.temperature.t_ssbef_col),
        ("t_grnd",      H.S.temperature.t_grnd_col,       Sd.temperature.t_grnd_col),
        ("emg",         H.S.temperature.emg_col,          Sd.temperature.emg_col),
        ("htvp",        H.S.energyflux.htvp_col,          Sd.energyflux.htvp_col),
        ("thv",         H.S.temperature.thv_col,          Sd.temperature.thv_col),
        ("beta",        H.S.temperature.beta_col,         Sd.temperature.beta_col),
        ("zii",         H.S.col_data.zii,                 Sd.col_data.zii),
        ("emv",         H.S.temperature.emv_patch,        Sd.temperature.emv_patch),
        ("thm",         H.S.temperature.thm_patch,        Sd.temperature.thm_patch),
        ("eflx_sh_tot", H.S.energyflux.eflx_sh_tot_patch, Sd.energyflux.eflx_sh_tot_patch),
        ("cgrnd",       H.S.energyflux.cgrnd_patch,       Sd.energyflux.cgrnd_patch),
        ("cgrndl",      H.S.energyflux.cgrndl_patch,      Sd.energyflux.cgrndl_patch),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-13s CPU reference is all-NaN/Inf — skipping (no parity signal)\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-13s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end

    # sanity: t_grnd snow-branch must be a weighted blend (frac_sno_eff>0)
    @printf("\n  CPU t_grnd[1] = %.4f   htvp[1] = %.3e (HSUB=%.3e HVAP=%.3e)\n",
            Array(H.S.temperature.t_grnd_col)[1], Array(H.S.energyflux.htvp_col)[1],
            Float64(CLM.HSUB), Float64(CLM.HVAP))

    println()
    println(nfail == 0 ? "  pre_flux kernels MATCH CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
