using NCDatasets, Dates, Printf, Statistics, CLM

# Quick diagnostic: run 2 days in July and print intermediate radiation values
const basedir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
fsurdat  = joinpath(basedir, "settings/CLM/parameters/surfdata_clm.nc")
paramfile = joinpath(basedir, "settings/CLM/parameters/clm5_params.nc")
fforcing = joinpath(basedir, "data/forcing/CLM_input/clmforc.2002.nc")

# Initialize
(inst, bounds, filt, tm) = CLM.clm_initialize!(fsurdat=fsurdat, paramfile=paramfile,
    start_date=DateTime(2002, 7, 1), dtime=1800)

# Read forcing
ff = NCDataset(fforcing)
# Find the step index for July 1 (day 182, half-hourly)
# Step 1 = Jan 1 00:00, so Jul 1 00:00 = step (181 * 48) + 1 = 8689
step_start = 181 * 48 + 1

# Set some state to reasonable summer values
nc = bounds.endc
np = bounds.endp
nlevsno = CLM.varpar.nlevsno

# Set initial temperatures to ~280K (summer ground)
inst.temperature.t_grnd_col .= 273.15
inst.temperature.t_soisno_col .= 273.15
inst.temperature.t_veg_patch .= 280.0

# Check what the forcing provides
forc_sw = ff["FSDS"]
println("Forcing FSDS dims: ", dimnames(forc_sw))
vals = Float64[]
for t in step_start:step_start+47
    v = forc_sw[1,1,t]
    push!(vals, ismissing(v) ? NaN : Float64(v))
end
println("FSDS range for day 182: ", minimum(vals), " to ", maximum(vals),
        " mean=", mean(vals))

# Check if direct/diffuse components exist
for vname in ["FSDSVD", "FSDSVN", "FSDSND", "FSDSNI", "FSDSdir", "FSDSdif"]
    if haskey(ff, vname)
        println("Has $vname")
    end
end

# Check what keys are in the forcing file
println("\nForcing file variables: ", join(keys(ff), ", "))

# Check the downscaled radiation used in the driver
println("\nDiagnostic: a2l forcing fields")
println("  forc_solad_downscaled_col size: ", size(inst.atm2lnd.forc_solad_downscaled_col))
println("  forc_solai_grc size: ", size(inst.atm2lnd.forc_solai_grc))

# Check albedo-related fields
sa = inst.surfalb
println("\nAlbedo fields:")
println("  albgrd_col size: ", size(sa.albgrd_col))
println("  albgri_col size: ", size(sa.albgri_col))
println("  fabd_patch size: ", size(sa.fabd_patch))
println("  fabi_patch size: ", size(sa.fabi_patch))

# Check patch/column structure
println("\nSubgrid structure:")
for p in 1:np
    @printf("  Patch %d: itype=%d, column=%d, wtgcell=%.4f\n",
            p, inst.patch.itype[p], inst.patch.column[p], inst.patch.wtgcell[p])
end
for c in 1:nc
    @printf("  Column %d: itype=%d, gridcell=%d, wtgcell=%.4f, snl=%d\n",
            c, inst.column.itype[c], inst.column.gridcell[c], inst.column.wtgcell[c], inst.column.snl[c])
end

# Check frac_sno
wdb = inst.water.waterdiagnosticbulk_inst
println("\nSnow state:")
for c in 1:nc
    @printf("  Column %d: frac_sno=%.4f, snow_depth=%.4f, h2osno=%.4f\n",
            c, wdb.frac_sno_col[c], wdb.snow_depth_col[c],
            inst.water.waterstatebulk_inst.ws.h2osno_no_layers_col[c])
end

# Now run one timestep and check radiation
println("\n--- Running first timestep at July 1 ---")
# Read forcing for this step
CLM.read_forcing_step!(inst.atm2lnd, ff, step_start + 12, inst.column, inst.gridcell)  # noon
CLM.downscale_forcings!(inst.atm2lnd, inst.column, inst.gridcell)

println("\nForcing after read+downscale (step $(step_start+12), noon Jul 1):")
@printf("  forc_solad_downscaled_col: %s\n", string(inst.atm2lnd.forc_solad_downscaled_col))
@printf("  forc_solai_grc: %s\n", string(inst.atm2lnd.forc_solai_grc))
@printf("  forc_t_downscaled_col: %s\n", string(inst.atm2lnd.forc_t_downscaled_col))

close(ff)
println("\nDone.")
