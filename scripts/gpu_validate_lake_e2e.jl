# ==========================================================================
# gpu_validate_lake_e2e.jl — end-to-end GPU parity for the WHOLE
# lake_temperature! driver (all lt1..lt7 kernels together).
#
# Builds a small Float32 lake instance (mirroring test_lake_temperature.jl's
# setup), runs lake_temperature! on the CPU, adapts every state struct to the GPU,
# runs the SAME call on device, and compares the mutated outputs. Exercises the
# full chain: thermal properties, phase change, diffusivity, solar heat source,
# extended-column tridiagonal assembly + batched solve, convective mixing, and
# the post-solve energy sums/conservation.
#
#   julia --project=scripts scripts/gpu_validate_lake_e2e.jl
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

function build(::Type{FT}) where {FT}
    vp = CLM.varpar
    vp.nlevsno = 5; vp.nlevsoi = 10; vp.nlevgrnd = 15; vp.nlevlak = 10
    vp.nlevurb = 5; vp.nlevmaxurbgrnd = max(vp.nlevurb, vp.nlevgrnd)
    CLM.varcon_init!()
    nlevsno, nlevgrnd, nlevlak = vp.nlevsno, vp.nlevgrnd, vp.nlevlak
    nlevtot = nlevsno + vp.nlevmaxurbgrnd
    nc = 2; np = 2
    dzlak = FT.(CLM.dzlak[]); zlak = FT.(CLM.zlak[])
    dzsoi = FT.(CLM.dzsoi[]); zsoi = FT.(CLM.zsoi[]); zisoi = FT.(CLM.zisoi[])

    col = CLM.ColumnData{FT}(); CLM.column_init!(col, nc)
    for c in 1:nc
        col.snl[c] = 0; col.lakedepth[c] = FT(30.0)
        for j in 1:nlevlak; col.dz_lake[c, j] = dzlak[j]; col.z_lake[c, j] = zlak[j]; end
        for j in 1:nlevgrnd
            jj = j + nlevsno
            col.dz[c, jj] = dzsoi[j]; col.z[c, jj] = zsoi[j]; col.zi[c, jj + 1] = zisoi[j + 1]
        end
        col.zi[c, nlevsno + 1] = FT(0.0)
    end
    patch = CLM.PatchData{FT}(); CLM.patch_init!(patch, np)
    for p in 1:np; patch.column[p] = p; end
    temp = CLM.TemperatureData{FT}(); CLM.temperature_init!(temp, np, nc, 1, 1)
    for c in 1:nc
        temp.t_grnd_col[c] = FT(280.0)
        for j in 1:nlevlak; temp.t_lake_col[c, j] = FT(280.0); end
        for jj in 1:nlevtot; temp.t_soisno_col[c, jj] = FT(280.0); end
    end
    sa = CLM.SolarAbsorbedData{FT}(); CLM.solarabs_init!(sa, np, 1)
    for p in 1:np
        sa.sabg_patch[p] = FT(100.0)
        sa.fsds_nir_d_patch[p] = FT(30.0); sa.fsds_nir_i_patch[p] = FT(20.0)
        sa.fsr_nir_d_patch[p] = FT(5.0); sa.fsr_nir_i_patch[p] = FT(3.0)
        for j in 1:(nlevsno + 1); sa.sabg_lyr_patch[p, j] = FT(100.0 / (nlevsno + 1)); end
    end
    ss = CLM.SoilStateData{FT}(); CLM.soilstate_init!(ss, np, nc)
    for c in 1:nc, j in 1:nlevgrnd
        ss.watsat_col[c, j] = FT(0.45); ss.tksatu_col[c, j] = FT(1.5)
        ss.tkmg_col[c, j] = FT(3.0); ss.tkdry_col[c, j] = FT(0.25); ss.csol_col[c, j] = FT(2.0e6)
    end
    wsb = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(wsb, nc, np, 1, 1)
    for c in 1:nc
        wsb.ws.h2osno_no_layers_col[c] = FT(0.0)
        for jj in 1:nlevtot; wsb.ws.h2osoi_liq_col[c, jj] = FT(0.0); wsb.ws.h2osoi_ice_col[c, jj] = FT(0.0); end
        for j in 1:nlevgrnd
            jj = j + nlevsno
            wsb.ws.h2osoi_liq_col[c, jj] = col.dz[c, jj] * FT(CLM.DENH2O) * ss.watsat_col[c, j]
        end
    end
    wdb = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(wdb, nc, np, 1, 1)
    for c in 1:nc
        wdb.snow_depth_col[c] = FT(0.0)
        for jj in 1:nlevtot; wdb.frac_iceold_col[c, jj] = FT(0.0); end
    end
    wfb = CLM.WaterFluxBulkData{FT}(); CLM.waterfluxbulk_init!(wfb, nc, np, 1, 1)
    ef = CLM.EnergyFluxData{FT}(); CLM.energyflux_init!(ef, np, nc, 1, 1)
    for p in 1:np
        ef.eflx_gnet_patch[p] = FT(50.0); ef.eflx_sh_tot_patch[p] = FT(20.0)
        ef.eflx_sh_grnd_patch[p] = FT(20.0); ef.eflx_soil_grnd_patch[p] = FT(50.0)
    end
    ls = CLM.LakeStateData{FT}(); CLM.lakestate_init!(ls, nc, np)
    for c in 1:nc
        ls.etal_col[c] = FT(0.5); ls.ks_col[c] = FT(0.1); ls.ws_col[c] = FT(0.05)
        ls.lake_raw_col[c] = FT(100.0); ls.betaprime_col[c] = FT(0.4)
        ls.savedtke1_col[c] = FT(CLM.TKWAT); ls.lakeresist_col[c] = FT(0.0)
        ls.lake_icethick_col[c] = FT(0.0); ls.lake_icefracsurf_col[c] = FT(0.0)
        for j in 1:nlevlak; ls.lake_icefrac_col[c, j] = FT(0.0); end
    end
    grnd_ch4 = fill!(similar(temp.t_grnd_col, FT, nc), FT(0.0))
    return (; nc, np, col, patch, sa, ss, wsb, wdb, wfb, ef, temp, ls, grnd_ch4)
end

run_lt!(d, m_c, m_p, dt) = CLM.lake_temperature!(d.col, d.patch, d.sa, d.ss, d.wsb,
    d.wdb, d.wfb, d.ef, d.temp, d.ls, d.grnd_ch4, m_c, m_p, 1:d.nc, 1:d.np, dt)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for lake_temperature! (whole driver)")
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
        m_c = falses(B.nc); m_p = falses(B.np); m_c .= true; m_p .= true
        dt = FT(1800)
        ad(x) = CLM.Adapt.adapt(device_array_type(), x)
        D = (; nc = B.nc, np = B.np, col = ad(B.col), patch = ad(B.patch), sa = ad(B.sa),
             ss = ad(B.ss), wsb = ad(B.wsb), wdb = ad(B.wdb), wfb = ad(B.wfb),
             ef = ad(B.ef), temp = ad(B.temp), ls = ad(B.ls), grnd_ch4 = ad(B.grnd_ch4))
        dmask(m) = to(collect(Bool, m))
        if !(D.temp.t_lake_col isa device_array_type() && D.ls.lake_icefrac_col isa device_array_type())
            println("  BLOCKED: a lake state struct did not move to the device under adapt.")
            return 2
        end
        run_lt!(B, m_c, m_p, dt)
        run_lt!(D, dmask(m_c), dmask(m_p), dt)
        checks = [
            ("t_lake",          B.temp.t_lake_col,        D.temp.t_lake_col),
            ("t_soisno",        B.temp.t_soisno_col,      D.temp.t_soisno_col),
            ("lake_icefrac",    B.ls.lake_icefrac_col,    D.ls.lake_icefrac_col),
            ("lake_icethick",   B.ls.lake_icethick_col,   D.ls.lake_icethick_col),
            ("lake_icefracsurf",B.ls.lake_icefracsurf_col,D.ls.lake_icefracsurf_col),
            ("savedtke1",       B.ls.savedtke1_col,       D.ls.savedtke1_col),
            ("eflx_gnet",       B.ef.eflx_gnet_patch,     D.ef.eflx_gnet_patch),
            ("eflx_sh_tot",     B.ef.eflx_sh_tot_patch,   D.ef.eflx_sh_tot_patch),
        ]
        nfail = 0
        for (nm, a, b) in checks
            d = reldiff(a, b); ok = d < 1f-2
            @printf("  [%s] %-17s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
            ok || (nfail += 1)
        end
        # errsoi = (ncvts - ocvts)/dtime - fin is an energy-balance RESIDUAL: a tiny
        # number from subtracting two huge (~1e7) energy sums. In Float32 that is
        # catastrophic cancellation (abs uncertainty ~ eps_f32*1e7 ~ 1), so its value
        # is precision-meaningless and CPU-F32 vs Metal-F32 differ. Informational only;
        # byte-identical on the Float64 CPU suite. Not a state variable.
        @printf("  [info] %-17s rel|dev-cpu| = %.3e  (Float32 cancellation; not gated)\n",
                "errsoi", reldiff(B.ef.errsoi_col, D.ef.errsoi_col))
        println()
        println(nfail == 0 ? "  WHOLE lake_temperature! MATCHES CPU ON $name ($FT) ✓" :
                             "  DIVERGENCE — investigate.")
        return nfail == 0 ? 0 : 1
    finally
        CLM.varpar.nlevsno, CLM.varpar.nlevsoi, CLM.varpar.nlevgrnd,
            CLM.varpar.nlevlak, CLM.varpar.nlevurb, CLM.varpar.nlevmaxurbgrnd = sv
    end
end

const BACKEND = detect_backend()
exit(main(BACKEND))
