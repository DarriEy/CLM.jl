# ==========================================================================
# gpu_validate_building_temp.jl — whole-routine Metal parity for the prognostic
# urban interior building-temperature solve (building_temperature!,
# UrbBuildTempOleson2015). Builds a single urban landunit (roof + sunwall +
# shadewall columns), runs building_temperature! on the host and on Metal, and
# asserts the 5 interior temperatures + 4 building energy fluxes match.
#
#   julia --project=scripts scripts/gpu_validate_building_temp.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
mf(x) = mf(Metal.MtlArray, x)

# ---- set module globals for the prognostic-building-temp scenario -----------
function set_globals!()
    vp = CLM.varpar
    vp.nlevsno = 5; vp.nlevgrnd = 10; vp.nlevurb = 5
    vp.nlevmaxurbgrnd = 10; vp.nlevsoi = 8
    CLM.urban_ctrl.read_namelist        = true
    CLM.urban_ctrl.building_temp_method = CLM.BUILDING_TEMP_METHOD_PROG
    CLM.urban_ctrl.urban_hac            = CLM.URBAN_HAC_ON
    CLM.urban_ctrl.urban_explicit_ac    = true
    return nothing
end

# ---- fixture: one urban landunit, columns roof/sunwall/shadewall ------------
function build_fixture(tbmax_val::Float64 = 297.0)
    vp = CLM.varpar
    nlevsno = vp.nlevsno; nlevurb = vp.nlevurb
    nlevurb_jj = nlevurb + nlevsno
    nc = 3; nl = 1; ng = 1; np = 1
    col = CLM.ColumnData{Float64}();     CLM.column_init!(col, nc)
    lun = CLM.LandunitData{Float64}();   CLM.landunit_init!(lun, nl)
    temp = CLM.TemperatureData{Float64}(); CLM.temperature_init!(temp, np, nc, nl, ng)
    ef  = CLM.EnergyFluxData{Float64}(); CLM.energyflux_init!(ef, np, nc, nl, ng)
    up  = CLM.UrbanParamsData{Float64}(); CLM.urbanparams_init!(up, nl; nlevurb=nlevurb)

    lun.urbpoi[1]       = true
    lun.itype[1]        = CLM.ISTURB_MD
    lun.canyon_hwr[1]   = 1.0
    lun.wtlunit_roof[1] = 0.5
    lun.ht_roof[1]      = 10.0
    up.t_building_min[1] = 290.0

    col.landunit .= 1
    col.itype[1] = CLM.ICOL_ROOF
    col.itype[2] = CLM.ICOL_SUNWALL
    col.itype[3] = CLM.ICOL_SHADEWALL
    col.snl     .= 0
    for c in 1:nc
        col.z[c, nlevurb_jj]      = 0.20
        col.zi[c, nlevurb_jj + 1] = 0.25
        temp.t_soisno_col[c, nlevurb_jj] = 291.0
        temp.t_ssbef_col[c, nlevurb_jj]  = 290.5
    end
    temp.t_roof_inner_lun[1] = 292.0
    temp.t_sunw_inner_lun[1] = 292.5
    temp.t_shdw_inner_lun[1] = 291.5
    temp.t_floor_lun[1]      = 293.0
    temp.t_building_lun[1]   = 294.0
    temp.taf_lun[1]          = 300.0

    tk = zeros(Float64, nc, nlevsno + vp.nlevmaxurbgrnd)
    for c in 1:nc
        tk[c, nlevurb_jj] = 1.5
    end
    tbmax = fill(tbmax_val, nl)
    p_ac  = fill(1.0, nl)
    mask_urbanl  = trues(nl)
    mask_nolakec = trues(nc)
    return col, lun, temp, ef, up, tk, tbmax, p_ac, mask_urbanl, mask_nolakec
end

outs(temp, ef) = (
    t_roof_inner = temp.t_roof_inner_lun, t_sunw_inner = temp.t_sunw_inner_lun,
    t_shdw_inner = temp.t_shdw_inner_lun, t_floor = temp.t_floor_lun,
    t_building = temp.t_building_lun, eflx_building = ef.eflx_building_lun,
    eflx_urban_ac = ef.eflx_urban_ac_lun, eflx_urban_heat = ef.eflx_urban_heat_lun,
    eflx_ventilation = ef.eflx_ventilation_lun)

function main(backend)
    println("="^72)
    println("Metal parity — building_temperature! (urban prognostic 5x5 solve)")
    println("="^72)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s (%s)\n\n", name, FT)
    set_globals!()
    dtime = 1800.0
    nfail = 0

    # Two scenarios: "cool" (t_building < setpoint → no AC clamp) and "hot"
    # (low setpoint → explicit-AC upper clamp branch fires: t_building reset +
    # eflx_urban_ac > 0). Exercises both the raw 5x5 solve and the HAC clamp.
    for (label, tbmax_val) in (("cool (no AC clamp)", 297.0), ("hot (explicit AC clamp)", 291.0))
        println("  --- scenario: $label ---")
        # ---- host reference ----
        col, lun, temp, ef, up, tk, tbmax, p_ac, mu, mn = build_fixture(tbmax_val)
        CLM.building_temperature!(col, lun, temp, ef, up, tbmax, p_ac, tk, mu, mn, 1:1, dtime)
        H = outs(temp, ef)

        # ---- device ----
        col2, lun2, temp2, ef2, up2, tk2, tbmax2, p_ac2, mu2, mn2 = build_fixture(tbmax_val)
        colD = mf(col2); lunD = mf(lun2); tempD = mf(temp2); efD = mf(ef2); upD = mf(up2)
        tkD = mf(tk2); tbmaxD = mf(tbmax2); p_acD = mf(p_ac2); muD = mf(mu2); mnD = mf(mn2)
        tempD.t_soisno_col isa Metal.MtlArray || (println("  BLOCKED: adapt failed."); return 2)
        try
            CLM.building_temperature!(colD, lunD, tempD, efD, upD, tbmaxD, p_acD, tkD,
                                      muD, mnD, 1:1, dtime)
            Metal.synchronize()
        catch e
            println("\n  ✗ DEVICE FAILURE:")
            println("      ", first(split(sprint(showerror, e), "\n")))
            for fr in stacktrace(catch_backtrace())
                s = string(fr.file)
                if occursin("/src/", s) && occursin("CLM", s)
                    println("      at ", basename(s), ":", fr.line, "  ", fr.func)
                end
            end
            return 1
        end
        D = outs(tempD, efD)

        for k in keys(H)
            h = Float64(Array(H[k])[1]); d = Float64(Array(D[k])[1])
            # Fields not written in the taken branch (e.g. eflx_urban_heat under the
            # AC-clamp path) keep the fixture's uninitialized NaN on BOTH backends —
            # matching NaN is agreement, not divergence.
            if !isfinite(h) || !isfinite(d)
                ok = (isnan(h) == isnan(d))
                @printf("  [%s] %-18s host=%s dev=%s (untouched — both non-finite)\n",
                        ok ? "PASS" : "FAIL", k, string(h), string(d))
            else
                rel = abs(h - d) / (1.0 + abs(h))
                ok = rel < 1f-3
                @printf("  [%s] %-18s host=%.6f dev=%.6f rel=%.2e\n", ok ? "PASS" : "FAIL", k, h, d, rel)
            end
            ok || (nfail += 1)
        end
        println()
    end
    println(nfail == 0 ? "  building_temperature! MATCHES host on $name (both scenarios)" :
                         "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
exit(main(detect_backend()))
