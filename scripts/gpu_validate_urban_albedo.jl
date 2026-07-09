# ==========================================================================
# gpu_validate_urban_albedo.jl — WHOLE-function Metal parity for the urban
# canyon radiative-transfer kernels: incident_direct!, incident_diffuse!,
# net_solar! (iterative multiple-reflection solve), and wasteheat!.
#
#   julia --project=scripts scripts/gpu_validate_urban_albedo.jl
# ==========================================================================
using CLM, Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
mf(x) = mf(Metal.MtlArray, x)

function reldiff(H, D)
    A = Array(H); B = Array(D); m = 0.0; n = 0
    for i in eachindex(A, B)
        (isfinite(A[i]) && isfinite(B[i])) || continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + abs(Float64(A[i])))); n += 1
    end
    return m, n
end

function main(backend)
    println("="^64)
    println("Urban albedo canyon-RT kernels — host vs device parity")
    println("="^64)
    backend === nothing && (println("  No GPU backend."); return 0)
    name, _, FT = backend
    @printf("  Backend: %s (%s)\n\n", name, FT)

    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nl = 2; numrad = CLM.NUMRAD
    coszen     = [0.5, 0.7]
    canyon_hwr = [1.0, 1.5]
    zen        = acos.(coszen)
    wtroad_perv = [0.4, 0.3]
    mask = trues(nl)

    # View factors from canyon_hwr (same formula as urbanparams_populate!)
    up = CLM.UrbanParamsData(); CLM.urbanparams_init!(up, nl; numrad=numrad)
    for l in 1:nl
        hwr = canyon_hwr[l]
        up.vf_sr[l] = sqrt(hwr^2 + 1.0) - hwr
        up.vf_wr[l] = 0.5 * (1.0 - up.vf_sr[l])
        up.vf_sw[l] = 0.5 * (hwr + 1.0 - sqrt(hwr^2 + 1.0)) / hwr
        up.vf_rw[l] = up.vf_sw[l]
        up.vf_ww[l] = 1.0 - up.vf_sw[l] - up.vf_rw[l]
        for ib in 1:numrad
            up.alb_roof_dir[l, ib] = 0.20; up.alb_roof_dif[l, ib] = 0.20
            up.alb_improad_dir[l, ib] = 0.10; up.alb_improad_dif[l, ib] = 0.10
            up.alb_perroad_dir[l, ib] = 0.15; up.alb_perroad_dif[l, ib] = 0.15
            up.alb_wall_dir[l, ib] = 0.25; up.alb_wall_dif[l, ib] = 0.25
        end
    end

    sdir = ones(nl, numrad); sdif = ones(nl, numrad)
    mk_out() = zeros(nl, numrad)
    nfail = 0; ncmp = 0

    # ---- incident_direct! ----
    hr, hs, hh = mk_out(), mk_out(), mk_out()
    CLM.incident_direct!(mask, canyon_hwr, coszen, zen, sdir, hr, hs, hh)
    dr, ds, dh = mf(mk_out()), mf(mk_out()), mf(mk_out())
    CLM.incident_direct!(mf(mask), mf(canyon_hwr), mf(coszen), mf(zen), mf(sdir), dr, ds, dh)
    Metal.synchronize()
    for (nm, H, D) in (("id_road", hr, dr), ("id_sunwall", hs, ds), ("id_shadewall", hh, dh))
        r, n = reldiff(H, D); ncmp += n; ok = r < 1f-3; ok || (nfail += 1)
        @printf("  [%s] %-14s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", nm, r, n)
    end

    # ---- incident_diffuse! ----
    fr, fs, fh = mk_out(), mk_out(), mk_out()
    CLM.incident_diffuse!(mask, canyon_hwr, sdif, fr, fs, fh, up)
    dfr, dfs, dfh = mf(mk_out()), mf(mk_out()), mf(mk_out())
    CLM.incident_diffuse!(mf(mask), mf(canyon_hwr), mf(sdif), dfr, dfs, dfh, mf(up))
    Metal.synchronize()
    for (nm, H, D) in (("if_road", fr, dfr), ("if_sunwall", fs, dfs), ("if_shadewall", fh, dfh))
        r, n = reldiff(H, D); ncmp += n; ok = r < 1f-3; ok || (nfail += 1)
        @printf("  [%s] %-14s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", nm, r, n)
    end

    # ---- net_solar! (uses incident outputs; sref outputs + solarabs sabs) ----
    sref_h = [ones(nl, numrad) for _ in 1:2]; append!(sref_h, [zeros(nl, numrad) for _ in 1:8])
    # order: improad_dir, perroad_dir, sunwall_dir, shadewall_dir, roof_dir, improad_dif...
    srh = [zeros(nl, numrad) for _ in 1:10]
    saH = CLM.SolarAbsorbedData(); CLM.solarabs_init!(saH, nl, nl); CLM.solarabs_init_cold!(saH, 1:nl)
    CLM.net_solar!(mask, coszen, canyon_hwr, wtroad_perv, sdir, sdif,
        up.alb_improad_dir, up.alb_perroad_dir, up.alb_wall_dir, up.alb_roof_dir,
        up.alb_improad_dif, up.alb_perroad_dif, up.alb_wall_dif, up.alb_roof_dif,
        hr, hs, hh, fr, fs, fh,
        srh[1], srh[2], srh[3], srh[4], srh[5], srh[6], srh[7], srh[8], srh[9], srh[10],
        up, saH)
    srd = [mf(zeros(nl, numrad)) for _ in 1:10]
    saD = mf(saH); upd = mf(up)
    CLM.net_solar!(mf(mask), mf(coszen), mf(canyon_hwr), mf(wtroad_perv), mf(sdir), mf(sdif),
        upd.alb_improad_dir, upd.alb_perroad_dir, upd.alb_wall_dir, upd.alb_roof_dir,
        upd.alb_improad_dif, upd.alb_perroad_dif, upd.alb_wall_dif, upd.alb_roof_dif,
        mf(hr), mf(hs), mf(hh), mf(fr), mf(fs), mf(fh),
        srd[1], srd[2], srd[3], srd[4], srd[5], srd[6], srd[7], srd[8], srd[9], srd[10],
        upd, saD)
    Metal.synchronize()
    snames = ["sref_improad_dir","sref_perroad_dir","sref_sunwall_dir","sref_shadewall_dir","sref_roof_dir",
              "sref_improad_dif","sref_perroad_dif","sref_sunwall_dif","sref_shadewall_dif","sref_roof_dif"]
    for i in 1:10
        r, n = reldiff(srh[i], srd[i]); ncmp += n; ok = r < 1f-3; ok || (nfail += 1)
        @printf("  [%s] %-20s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", snames[i], r, n)
    end
    for (nm, hf, df) in (("sabs_roof_dir", saH.sabs_roof_dir_lun, saD.sabs_roof_dir_lun),
                         ("sabs_sunwall_dir", saH.sabs_sunwall_dir_lun, saD.sabs_sunwall_dir_lun),
                         ("sabs_improad_dif", saH.sabs_improad_dif_lun, saD.sabs_improad_dif_lun))
        r, n = reldiff(hf, df); ncmp += n; ok = r < 1f-3; ok || (nfail += 1)
        @printf("  [%s] %-20s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", nm, r, n)
    end

    # ---- wasteheat! ----
    CLM.urban_ctrl.read_namelist = true
    CLM.urban_ctrl.building_temp_method = CLM.BUILDING_TEMP_METHOD_PROG
    CLM.urban_ctrl.urban_hac = CLM.URBAN_WASTEHEAT_ON
    lun = CLM.LandunitData(); CLM.landunit_init!(lun, nl)
    for l in 1:nl
        lun.gridcell[l] = 1; lun.wtlunit_roof[l] = 0.5; lun.canyon_hwr[l] = canyon_hwr[l]
    end
    filt = collect(1:nl)
    whr = [10.0, 12.0]; whs = [5.0, 6.0]; whh = [4.0, 3.0]
    hcr = [2.0, 1.0]; hcs = [1.0, 2.0]; hch = [0.5, 0.5]
    mk_ef() = (ef = CLM.EnergyFluxData();
               ef.eflx_wasteheat_lun = zeros(nl); ef.eflx_heat_from_ac_lun = zeros(nl);
               ef.eflx_urban_ac_lun = [3.0, 4.0]; ef.eflx_urban_heat_lun = [2.0, 1.0]; ef)
    efH = mk_ef()
    CLM.wasteheat!(efH, lun, nl, filt, whr, whs, whh, hcr, hcs, hch)
    efB = mk_ef(); efD = mf(efB); lunD = mf(lun)
    CLM.wasteheat!(efD, lunD, nl, mf(filt), mf(whr), mf(whs), mf(whh), mf(hcr), mf(hcs), mf(hch))
    Metal.synchronize()
    for (nm, hf, df) in (("eflx_wasteheat", efH.eflx_wasteheat_lun, efD.eflx_wasteheat_lun),
                         ("eflx_heat_from_ac", efH.eflx_heat_from_ac_lun, efD.eflx_heat_from_ac_lun))
        r, n = reldiff(hf, df); ncmp += n; ok = r < 1f-3; ok || (nfail += 1)
        @printf("  [%s] %-20s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", nm, r, n)
    end

    println()
    println(nfail == 0 ? "  URBAN CANYON-RT + wasteheat kernels MATCH host on $name over $ncmp finite outputs" :
                          "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
exit(main(detect_backend()))
