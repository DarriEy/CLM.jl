# =============================================================================
# Coupled end-to-end validation of the LUNA ×1e5 unit fix (commit 56152a3).
#
# Free-runs the Bow CN+PHS+LUNA driver from the n1757845 IC dump and, at every
# step, compares Julia's LUNA-updated vcmx25_z (tree patch 2, grass patch 3) to
# the Fortran after_hydrologydrainage dump's vcmx25_z. The LUNA end-of-day
# acclimation update fires once in this window (~n1757856, local midnight); the
# fix changes that update's target, so this is the direct coupled test that the
# corrected reference-NUE keeps vcmx25_z tracking Fortran (not just the
# standalone 1-day update).
#
# Usage: julia +1.12 --project=. scripts/validate_luna_coupled.jl
# =============================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const SUM   = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const N0    = 1757845
const NST   = 28
const BASE  = 1753153

start_date = DateTime(2002, 1, 1) + Hour(N0 - BASE)
ffile = replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc")
(inst, bounds, filt, tm) = build_bow_inst(; dtime = 3600, start_date = start_date,
                                            use_cn = true, use_luna = true)
if isempty(inst.photosyns.vcmx25_z_patch)
    inst.photosyns.vcmx25_z_patch = fill(30.0, bounds.endp, CLM.NLEVCAN)
    inst.photosyns.jmx25_z_patch  = fill(60.0, bounds.endp, CLM.NLEVCAN)
end
icdump = joinpath(SUM, "pdump_before_step_n$(N0).nc")
inject_dump!(inst, bounds, icdump)
let ds = NCDataset(icdump, "r")
    if haskey(ds, "vegwp")
        vw = ds["vegwp"][:, :]
        for pd in 1:size(vw, 2), seg in 1:4
            inst.canopystate.vegwp_patch[pd, seg] = Float64(vw[seg, pd])
        end
    end
    close(ds)
end

config = CLM.CLMDriverConfig(use_cn = true, use_aquifer_layer = false,
                             use_hydrstress = true, use_luna = true)
filt_ia = CLM.clump_filter_inactive_and_active
ng, nc, np = bounds.endg, bounds.endc, bounds.endp
fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, ffile)
tf = replace(ffile, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
if isfile(tf)
    dt = NCDataset(tf, "r")
    if haskey(dt, "TOPO")
        ft = Float64(dt["TOPO"][1])
        for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
        for c in 1:nc; inst.topo.topo_col[c] = ft; end
    end
    close(dt)
end

ps = inst.photosyns
@printf("%-5s | %-26s | %-26s | %s\n", "step", "vcmx25_z tree (Jl/Ft/ratio)",
        "vcmx25_z grass (Jl/Ft/ratio)", "EOD?")
worst = 0.0
for i in 0:(NST - 1)
    cur = start_date + Hour(i)
    calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
    nextsw = calday + 3600.0 / CLM.SECSPDAY
    (declinm1, _) = CLM.compute_orbital(calday - 3600.0 / CLM.SECSPDAY)
    CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:bounds.endg)
    CLM.advance_timestep!(tm)
    CLM.read_forcing_step!(fr, inst.atm2lnd, cur, ng, nc)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    (yr, mon, d, tod) = CLM.get_curr_date(tm)
    eod = CLM.is_end_curr_day(tm)
    CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
        CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
        nstep = tm.nstep, is_first_step = false,
        is_beg_curr_day = CLM.is_beg_curr_day(tm), is_end_curr_day = eod,
        is_beg_curr_year = CLM.is_beg_curr_year(tm), dtime = 3600.0, mon = mon, day = d,
        photosyns = inst.photosyns)
    n = N0 + i
    da = NCDataset(joinpath(SUM, "pdump_after_hydrologydrainage_n$(n).nc"), "r")
    fv = da["vcmx25_z"][:, :]    # dims (levcan, pft) on file
    close(da)
    ft_tree  = Float64(fv[1, 2]); ft_grass = Float64(fv[1, 3])
    jl_tree  = Float64(ps.vcmx25_z_patch[2, 1]); jl_grass = Float64(ps.vcmx25_z_patch[3, 1])
    rt = jl_tree / ft_tree; rg = jl_grass / ft_grass
    global worst = max(worst, abs(rt - 1), abs(rg - 1))
    @printf("n%-4d | %7.3f /%7.3f /%6.4f | %7.3f /%7.3f /%6.4f | %s\n",
            i + 1, jl_tree, ft_tree, rt, jl_grass, ft_grass, rg, eod ? "EOD" : "")
end
CLM.forcing_reader_close!(fr)
@printf("\nworst |ratio-1| over window = %.4f\n", worst)
