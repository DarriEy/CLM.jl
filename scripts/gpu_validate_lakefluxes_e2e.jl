# ==========================================================================
# gpu_validate_lakefluxes_e2e.jl — end-to-end GPU parity for the WHOLE
# lake_fluxes! driver (the single per-patch _lake_fluxes_kernel! with its
# internal sequential Monin-Obukhov stability iteration).
#
# Builds a small Float32 lake instance, runs lake_fluxes! on the CPU, adapts
# every state struct to Metal (via the LakeFluxDV @adapt_structure device view),
# runs the SAME call on device, and compares the mutated outputs with reldiff.
# Exercises the full chain: fetch/roughness init, qsat + stability functions,
# the loop-carried Newton/stability iteration, the temperature corrections, and
# the final flux / momentum-stress / mixing-parameter writes.
#
# Two columns span both branches: c1 = open unfrozen water (snl=0, t_grnd>TFRZ),
# c2 = frozen lake with a snow layer (snl<0, t_grnd<TFRZ) so the frozen-roughness
# and snow-subsurface paths fire too.
#
# IMPORTANT (A1 lesson): reldiff silently PASSES when both sides are NaN, so this
# harness first asserts the CPU-F32 reference outputs are FINITE before trusting
# any parity number. Every required scalar input is set to a real finite value.
#
#   julia --project=scripts scripts/gpu_validate_lakefluxes_e2e.jl
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

function build(::Type{FT}) where {FT}
    vp = CLM.varpar
    vp.nlevsno = 5; vp.nlevsoi = 10; vp.nlevgrnd = 15; vp.nlevlak = 10
    vp.nlevurb = 5; vp.nlevmaxurbgrnd = max(vp.nlevurb, vp.nlevgrnd)
    CLM.varcon_init!()
    nlevsno, nlevgrnd, nlevlak = vp.nlevsno, vp.nlevgrnd, vp.nlevlak
    nlevtot = nlevsno + vp.nlevmaxurbgrnd
    nc = 2; np = 2; ng = 2; nl = 2
    dzlak = FT.(CLM.dzlak[]); zlak = FT.(CLM.zlak[])
    dzsoi = FT.(CLM.dzsoi[]); zsoi = FT.(CLM.zsoi[]); zisoi = FT.(CLM.zisoi[])

    # --- column ---
    col = CLM.ColumnData{FT}(); CLM.column_init!(col, nc)
    snls = (0, -2)  # c1 open water, c2 two snow layers
    for c in 1:nc
        col.gridcell[c] = c
        col.snl[c] = snls[c]
        col.lakedepth[c] = FT(30.0)
        for j in 1:nlevlak; col.dz_lake[c, j] = dzlak[j]; col.z_lake[c, j] = zlak[j]; end
        for j in 1:nlevgrnd
            jj = j + nlevsno
            col.dz[c, jj] = dzsoi[j]; col.z[c, jj] = zsoi[j]; col.zi[c, jj + 1] = zisoi[j + 1]
        end
        # snow layers: give the snow dz a real positive value so dzsur is finite
        for jj in 1:nlevsno; col.dz[c, jj] = FT(0.05); col.z[c, jj] = FT(-0.05 * (nlevsno - jj + 1)); end
        col.zi[c, nlevsno + 1] = FT(0.0)
    end

    # --- patch ---
    patch = CLM.PatchData{FT}(); CLM.patch_init!(patch, np)
    for p in 1:np
        patch.column[p] = p; patch.gridcell[p] = p; patch.landunit[p] = p
    end

    # --- landunit ---
    lun = CLM.LandunitData{FT}(); CLM.landunit_init!(lun, nl)

    # --- temperature ---
    temp = CLM.TemperatureData{FT}(); CLM.temperature_init!(temp, np, nc, 1, 1)
    tgs = (FT(285.0), FT(265.0))   # c1 unfrozen, c2 frozen
    tlk = (FT(283.0), FT(272.0))
    for c in 1:nc
        temp.t_grnd_col[c] = tgs[c]
        for j in 1:nlevlak; temp.t_lake_col[c, j] = tlk[c]; end
        for jj in 1:nlevtot; temp.t_soisno_col[c, jj] = (c == 1 ? FT(283.0) : FT(265.0)); end
    end

    # --- solar absorbed ---
    sa = CLM.SolarAbsorbedData{FT}(); CLM.solarabs_init!(sa, np, 1)
    for p in 1:np
        sa.sabg_patch[p] = FT(120.0)
        for j in 1:(nlevsno + 1); sa.sabg_lyr_patch[p, j] = FT(120.0 / (nlevsno + 1)); end
    end

    # --- water state ---
    wsb = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(wsb, nc, np, 1, 1)
    for c in 1:nc
        wsb.ws.h2osno_no_layers_col[c] = FT(0.0)
        for jj in 1:nlevtot; wsb.ws.h2osoi_liq_col[c, jj] = FT(0.0); wsb.ws.h2osoi_ice_col[c, jj] = FT(0.0); end
        # surface layer water: c1 liquid, c2 ice (so htvp picks HSUB for the frozen one)
        jtop = col.snl[c] + 1 + nlevsno
        if c == 1
            wsb.ws.h2osoi_liq_col[c, jtop] = FT(10.0)
        else
            wsb.ws.h2osoi_ice_col[c, jtop] = FT(10.0)
        end
    end
    wdb = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(wdb, nc, np, 1, 1)
    wfb = CLM.WaterFluxBulkData{FT}(); CLM.waterfluxbulk_init!(wfb, nc, np, 1, 1)

    # --- energy flux ---
    ef = CLM.EnergyFluxData{FT}(); CLM.energyflux_init!(ef, np, nc, 1, 1)

    # --- friction velocity ---
    fv = CLM.FrictionVelocityData{FT}(); CLM.frictionvel_init!(fv, np, nc)

    # --- lake state ---
    ls = CLM.LakeStateData{FT}(); CLM.lakestate_init!(ls, nc, np)
    for c in 1:nc
        ls.savedtke1_col[c] = FT(CLM.TKWAT)
        ls.ws_col[c] = FT(0.05); ls.ks_col[c] = FT(0.1)
        ls.betaprime_col[c] = FT(0.4)
    end

    # --- forcings (column-indexed) ---
    forc_t = fill(FT(0.0), nc); forc_th = fill(FT(0.0), nc); forc_q = fill(FT(0.0), nc)
    forc_pbot = fill(FT(0.0), nc); forc_rho = fill(FT(0.0), nc); forc_lwrad = fill(FT(0.0), nc)
    for c in 1:nc
        forc_t[c]     = (c == 1 ? FT(284.0) : FT(264.0))
        forc_th[c]    = (c == 1 ? FT(285.0) : FT(265.0))
        forc_q[c]     = FT(0.006)
        forc_pbot[c]  = FT(101325.0)
        forc_rho[c]   = FT(1.2)
        forc_lwrad[c] = FT(300.0)
    end
    # --- forcings (gridcell-indexed) ---
    forc_u = fill(FT(3.0), ng); forc_v = fill(FT(2.0), ng)
    forc_hgt_u = fill(FT(30.0), ng); forc_hgt_t = fill(FT(30.0), ng); forc_hgt_q = fill(FT(30.0), ng)
    grc_lat = fill(FT(0.893), ng)   # 51.2 deg N in radians — see run_lf! below

    return (; nc, np, col, patch, lun, temp, sa, wsb, wdb, wfb, ef, fv, ls,
        forc_t, forc_th, forc_q, forc_pbot, forc_rho, forc_lwrad,
        forc_u, forc_v, forc_hgt_u, forc_hgt_t, forc_hgt_q, grc_lat)
end

run_lf!(d, m_c, m_p) = CLM.lake_fluxes!(d.temp, d.ef, d.fv, d.sa, d.ls, d.wsb, d.wdb, d.wfb,
    d.col, d.patch, d.lun,
    d.forc_t, d.forc_th, d.forc_q, d.forc_pbot, d.forc_rho, d.forc_lwrad,
    d.forc_u, d.forc_v, d.forc_hgt_u, d.forc_hgt_t, d.forc_hgt_q,
    m_c, m_p, 1:d.nc, 1:d.np; dtime = 1800.0,
    # Non-zero latitude on purpose: ks = 6.6*sqrt(|sin(lat)|)*u2m^-1.84 is exactly
    # zero at lat = 0, so validating on the default would leave the ks branch (and
    # the sqrt/pow it lowers to on device) UNCOVERED. 0.893 rad = 51.2 deg, the
    # lake parity reference's own latitude.
    grc_lat = d.grc_lat)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for lake_fluxes! (whole per-patch driver)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)
    sv = (CLM.varpar.nlevsno, CLM.varpar.nlevsoi, CLM.varpar.nlevgrnd,
          CLM.varpar.nlevlak, CLM.varpar.nlevurb, CLM.varpar.nlevmaxurbgrnd)
    try
        B = build(FT)
        m_c = fill(true, B.nc); m_p = fill(true, B.np)
        ad(x) = CLM.Adapt.adapt(device_array_type(), x)
        D = (; nc = B.nc, np = B.np,
            col = ad(B.col), patch = ad(B.patch), lun = ad(B.lun), temp = ad(B.temp),
            sa = ad(B.sa), wsb = ad(B.wsb), wdb = ad(B.wdb), wfb = ad(B.wfb),
            ef = ad(B.ef), fv = ad(B.fv), ls = ad(B.ls),
            forc_t = to(B.forc_t), forc_th = to(B.forc_th), forc_q = to(B.forc_q),
            forc_pbot = to(B.forc_pbot), forc_rho = to(B.forc_rho), forc_lwrad = to(B.forc_lwrad),
            forc_u = to(B.forc_u), forc_v = to(B.forc_v),
            forc_hgt_u = to(B.forc_hgt_u), forc_hgt_t = to(B.forc_hgt_t), forc_hgt_q = to(B.forc_hgt_q),
            grc_lat = to(B.grc_lat))
        if !(D.temp.t_grnd_col isa device_array_type() && D.fv.ustar_patch isa device_array_type())
            println("  BLOCKED: a lake state struct did not move to the device under adapt.")
            return 2
        end
        dmask(m) = to(collect(Bool, m))
        run_lf!(B, m_c, m_p)
        run_lf!(D, dmask(m_c), dmask(m_p))

        # A1 lesson: reldiff(NaN, NaN) PASSES silently. Assert CPU-F32 outputs are
        # finite first so a false NaN==NaN pass cannot masquerade as parity.
        cpu_outs = [
            ("eflx_sh_tot",    B.ef.eflx_sh_tot_patch),
            ("eflx_sh_grnd",   B.ef.eflx_sh_grnd_patch),
            ("eflx_lh_tot",    B.ef.eflx_lh_tot_patch),
            ("eflx_lh_grnd",   B.ef.eflx_lh_grnd_patch),
            ("eflx_lwrad_out", B.ef.eflx_lwrad_out_patch),
            ("eflx_lwrad_net", B.ef.eflx_lwrad_net_patch),
            ("eflx_soil_grnd", B.ef.eflx_soil_grnd_patch),
            ("taux",           B.ef.taux_patch),
            ("tauy",           B.ef.tauy_patch),
            ("htvp",           B.ef.htvp_col),
            ("qflx_evap_soi",  B.wfb.wf.qflx_evap_soi_patch),
            ("qflx_evap_tot",  B.wfb.wf.qflx_evap_tot_patch),
            ("ustar",          B.fv.ustar_patch),
            ("z0mg",           B.fv.z0mg_patch),
            ("z0hg",           B.fv.z0hg_patch),
            ("z0qg",           B.fv.z0qg_patch),
            ("t_grnd",         B.temp.t_grnd_col),
            ("ws_col",         B.ls.ws_col),
            ("ks_col",         B.ls.ks_col),
        ]
        for (nm, a) in cpu_outs
            if !allfinite(a)
                @printf("  BLOCKED: CPU-F32 reference output %s is non-finite %s — parity numbers untrustworthy.\n",
                        nm, string(Array(a)))
                return 2
            end
        end
        println("  CPU-F32 reference outputs all finite ✓\n")

        checks = [
            ("eflx_sh_tot",    B.ef.eflx_sh_tot_patch,    D.ef.eflx_sh_tot_patch),
            ("eflx_sh_grnd",   B.ef.eflx_sh_grnd_patch,   D.ef.eflx_sh_grnd_patch),
            ("eflx_lh_tot",    B.ef.eflx_lh_tot_patch,    D.ef.eflx_lh_tot_patch),
            ("eflx_lh_grnd",   B.ef.eflx_lh_grnd_patch,   D.ef.eflx_lh_grnd_patch),
            ("eflx_lwrad_out", B.ef.eflx_lwrad_out_patch, D.ef.eflx_lwrad_out_patch),
            ("eflx_lwrad_net", B.ef.eflx_lwrad_net_patch, D.ef.eflx_lwrad_net_patch),
            ("eflx_soil_grnd", B.ef.eflx_soil_grnd_patch, D.ef.eflx_soil_grnd_patch),
            ("taux",           B.ef.taux_patch,           D.ef.taux_patch),
            ("tauy",           B.ef.tauy_patch,           D.ef.tauy_patch),
            ("htvp",           B.ef.htvp_col,             D.ef.htvp_col),
            ("qflx_evap_soi",  B.wfb.wf.qflx_evap_soi_patch, D.wfb.wf.qflx_evap_soi_patch),
            ("qflx_evap_tot",  B.wfb.wf.qflx_evap_tot_patch, D.wfb.wf.qflx_evap_tot_patch),
            ("ustar",          B.fv.ustar_patch,          D.fv.ustar_patch),
            ("z0mg",           B.fv.z0mg_patch,           D.fv.z0mg_patch),
            ("z0hg",           B.fv.z0hg_patch,           D.fv.z0hg_patch),
            ("z0qg",           B.fv.z0qg_patch,           D.fv.z0qg_patch),
            ("t_grnd",         B.temp.t_grnd_col,         D.temp.t_grnd_col),
            ("ws_col",         B.ls.ws_col,               D.ls.ws_col),
            ("ks_col",         B.ls.ks_col,               D.ls.ks_col),
        ]
        nfail = 0
        for (nm, a, b) in checks
            d = reldiff(a, b); ok = d < 1f-2
            @printf("  [%s] %-15s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
            ok || (nfail += 1)
        end
        println()
        println(nfail == 0 ? "  WHOLE lake_fluxes! MATCHES CPU ON $name ($FT) ✓" :
                             "  DIVERGENCE — investigate.")
        return nfail == 0 ? 0 : 1
    finally
        CLM.varpar.nlevsno, CLM.varpar.nlevsoi, CLM.varpar.nlevgrnd,
            CLM.varpar.nlevlak, CLM.varpar.nlevurb, CLM.varpar.nlevmaxurbgrnd = sv
    end
end

const BACKEND = detect_backend()
exit(main(BACKEND))
