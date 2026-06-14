# ==========================================================================
# spring_spinup_test.jl — does the summer BTRAN inversion survive a SPUN-UP init?
#
# verify_vs_fortran.jl COLD-STARTS Julia (pedotransfer soil water) and compares
# to a Fortran reference from a calibrated SPUN-UP run. To separate "cold-start
# initialization artifact" from "per-step physics bug", here we instead inject
# the Fortran spun-up restart (clm2.r.2003-01-01) as the IC and free-run Julia
# Jan 1 -> mid-July 2003, tracking btran / root-zone soil water / snow. If the
# summer btran stays LOW (~0.005, matching the known-dry Fortran July state),
# the annual "inversion" was the cold-start. If it climbs back to ~0.8, the
# soil re-wets and there is a real physics bug.
#
#   julia +1.12 --project=. scripts/spring_spinup_test.jl [endnstep]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using CLM, NCDatasets, Printf, Dates

const N_START = 8761                                   # 2003-01-01 00:00
const N_END   = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 13461   # mid-July
const RESTART = joinpath(DUMPDIR, "Bow_at_Banff_lumped.clm2.r.2003-01-01-00000.nc")

println("="^70)
println("  SPUN-UP free-run: inject Fortran restart, run n$(N_START)..$(N_END)")
println("="^70)
isfile(RESTART) || error("missing restart: $RESTART")

(inst, bounds, filt, tm) = build_bow_inst(; dtime=3600, start_date=DateTime(2003,1,1))
# Inject the Fortran spun-up prognostic state (soil water, temps, snow, ...)
CLM.read_fortran_restart!(RESTART, inst, bounds)

config  = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=false)
filt_ia = CLM.clump_filter_inactive_and_active
ng, nc, np = bounds.endg, bounds.endc, bounds.endp

fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, FFORCING)
topo_file = replace(FFORCING, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
if isfile(topo_file)
    ds_topo = NCDataset(topo_file, "r")
    if haskey(ds_topo, "TOPO")
        ft = Float64(ds_topo["TOPO"][1])
        for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
        for c in 1:nc; inst.topo.topo_col[c] = ft; end
    end
    close(ds_topo)
end

bt()   = inst.energyflux.btran_patch
hsno() = sum(inst.water.waterstatebulk_inst.ws.h2osoi_liq_col[1, 1:12]) +
         sum(inst.water.waterstatebulk_inst.ws.h2osoi_ice_col[1, 1:12])
# root-zone (top 10 soil layers) volumetric-ish liquid water
rzliq() = sum(inst.water.waterstatebulk_inst.ws.h2osoi_liq_col[1, 13:22])

@printf("%-18s %8s %8s %8s %10s %10s\n", "date", "btran2", "btran3", "snow_mm", "rootliq", "zwt")
first_decl = true
for n in N_START:N_END
    step_start = tm.current_date
    calday = CLM.get_curr_calday(tm)
    (declin, _) = CLM.compute_orbital(calday)
    obliqr = CLM.ORB_OBLIQR_DEFAULT
    nextsw_cday = calday + 3600.0 / CLM.SECSPDAY
    if first_decl
        (declinm1, _) = CLM.compute_orbital(calday - 3600.0/CLM.SECSPDAY)
        CLM.init_daylength!(inst.gridcell, declin, declinm1, obliqr, 1:bounds.endg)
        global first_decl = false
    end
    CLM.advance_timestep!(tm)
    CLM.read_forcing_step!(fr, inst.atm2lnd, step_start, ng, nc)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    (yr, mon, d, tod) = CLM.get_curr_date(tm)
    CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
                 true, nextsw_cday, declin, declin, obliqr, false, false, "", false;
                 nstep=tm.nstep, is_first_step=false,
                 is_beg_curr_day=CLM.is_beg_curr_day(tm),
                 is_end_curr_day=CLM.is_end_curr_day(tm),
                 is_beg_curr_year=CLM.is_beg_curr_year(tm),
                 dtime=3600.0, mon=mon, day=d, photosyns=inst.photosyns)
    if n % 480 == 0 || n == N_END   # ~every 20 days + the end
        b = bt()
        @printf("%-18s %8.4f %8.4f %8.1f %10.2f %10.4f\n",
                string(step_start), b[2], b[3], hsno(), rzliq(), inst.soilhydrology.zwt_col[1])
    end
end
CLM.forcing_reader_close!(fr)

# Compare the end state to the Fortran dump at N_END (known-dry July)
fdump = joinpath(DUMPDIR, "pdump_after_hydrologydrainage_n$(N_END).nc")
if isfile(fdump)
    ds = NCDataset(fdump, "r"); rd(v) = ismissing(v) ? NaN : Float64(v)
    println("\n--- end-state vs Fortran n$N_END dump ---")
    @printf("  btran (Julia free-run)   p2=%.4f p3=%.4f\n", bt()[2], bt()[3])
    @printf("  Fortran summer btran is ~0.005-0.009 (known dry); cold-start Julia annual was ~0.82\n")
    fliq = sum(rd(ds["H2OSOI_LIQ"][k,1]) for k in 13:22)
    @printf("  root-zone liq  Julia=%.2f  Fortran=%.2f mm\n", rzliq(), fliq)
    close(ds)
end
