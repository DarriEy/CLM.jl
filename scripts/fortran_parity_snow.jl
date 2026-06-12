# ==========================================================================
# fortran_parity_snow.jl — pinpoint the INT_SNOW single-step divergence.
# Prints snow-pack quantities for col 1 before and after one Julia step, with
# the Fortran before/after values for comparison.
#
#   julia +1.12 --project=. scripts/fortran_parity_snow.jl [nstep]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const NSTEP       = length(ARGS) >= 1 ? ARGS[1] : "8761"
const DUMP_BEFORE = joinpath(DUMPDIR, "pdump_before_step_n$(NSTEP).nc")
const DUMP_AFTER  = joinpath(DUMPDIR, "pdump_after_hydrologydrainage_n$(NSTEP).nc")

fval(f, v) = begin
    ds = NCDataset(f, "r"); x = haskey(ds, v) ? Float64(ds[v][:][1]) : NaN; close(ds); x
end

println("="^64); println("  INT_SNOW snow-accumulation divergence  (nstep=$NSTEP)"); println("="^64)

(inst, bounds, filt, tm) = build_bow_inst(; dtime=3600, start_date=DateTime(2003,1,1))
inject_dump!(inst, bounds, DUMP_BEFORE)

ws = inst.water.waterstatebulk_inst
wd = inst.water.waterdiagnosticbulk_inst
snowsum() = sum(ws.ws.h2osoi_liq_col[1, j] + ws.ws.h2osoi_ice_col[1, j] for j in 1:CLM.varpar.nlevsno)

@printf("\n  %-22s %14s %14s\n", "quantity", "Julia", "Fortran")
@printf("  %s\n", "-"^52)
@printf("  %-22s %14.4f %14s\n", "n_melt[1]", inst.scf_method.n_melt[1], "—")
@printf("  %-22s %14.4f %14.4f\n", "int_snow (before)", ws.int_snow_col[1], fval(DUMP_BEFORE, "INT_SNOW"))
@printf("  %-22s %14.4f %14.4f\n", "h2osno_layers (before)", snowsum(), NaN)
@printf("  %-22s %14d %14.0f\n", "snl (before)", inst.column.snl[1], fval(DUMP_BEFORE, "SNLSNO"))
@printf("  %-22s %14.4f %14.4f\n", "snow_depth (before)", wd.snow_depth_col[1], fval(DUMP_BEFORE, "SNOW_DEPTH"))
@printf("  %-22s %14.4f %14.4f\n", "frac_sno (before)", wd.frac_sno_col[1], fval(DUMP_BEFORE, "frac_sno"))

# ---- run one step (mirror clm_run.jl loop body) ----
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
# Fortran step 8761 uses the forcing at the STEP START (PRECTmms[0]); Julia's
# clm_run loop advances first then reads, landing one interval late. Read at the
# step-start time to match Fortran.
step_start = tm.current_date
CLM.advance_timestep!(tm)
CLM.read_forcing_step!(fr, inst.atm2lnd, step_start, ng, nc)
CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
calday = CLM.get_curr_calday(tm); (declin, eccf) = CLM.compute_orbital(calday)
(yr, mon, d, tod) = CLM.get_curr_date(tm)

# snowfall forcing this step (downscaled, as the snow model sees it)
qsnow = inst.atm2lnd.forc_snow_downscaled_col[1]
qrain = inst.atm2lnd.forc_rain_downscaled_col[1]
@printf("  %-22s %14.6e %14s\n", "forc_snow (mm/s)", qsnow, "—")
@printf("  %-22s %14.6e %14s\n", "forc_rain (mm/s)", qrain, "—")
@printf("  %-22s %14.6e %14s\n", "newsnow=qsnow*dt (mm)", qsnow*3600.0, "Fortran≈0")

CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true,
             calday + 3600.0/CLM.SECSPDAY, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
             false, false, "", false; nstep=tm.nstep, is_first_step=(tm.nstep==1),
             is_beg_curr_day=CLM.is_beg_curr_day(tm), is_end_curr_day=CLM.is_end_curr_day(tm),
             is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=3600.0, mon=mon, day=d,
             photosyns=inst.photosyns)
CLM.forcing_reader_close!(fr)

println("  " * "-"^52)
@printf("  %-22s %14.4f %14.4f\n", "int_snow (AFTER)", ws.int_snow_col[1], fval(DUMP_AFTER, "INT_SNOW"))
@printf("  %-22s %14.4f %14.4f\n", "h2osno_layers (AFTER)", snowsum(), NaN)
@printf("  %-22s %14d %14.0f\n", "snl (AFTER)", inst.column.snl[1], fval(DUMP_AFTER, "SNLSNO"))
@printf("  %-22s %14.4f %14.4f\n", "snow_depth (AFTER)", wd.snow_depth_col[1], fval(DUMP_AFTER, "SNOW_DEPTH"))
@printf("  %-22s %14.4f %14.4f\n", "frac_sno (AFTER)", wd.frac_sno_col[1], fval(DUMP_AFTER, "frac_sno"))
@printf("  %-22s %14.6e %14s\n", "qflx_snow_grnd (AFTER)", inst.water.waterfluxbulk_inst.wf.qflx_snow_grnd_col[1], "—")
println()
@printf("  Δint_snow: Julia %+.3f  Fortran %+.3f\n",
        ws.int_snow_col[1]-fval(DUMP_BEFORE,"INT_SNOW"), fval(DUMP_AFTER,"INT_SNOW")-fval(DUMP_BEFORE,"INT_SNOW"))
