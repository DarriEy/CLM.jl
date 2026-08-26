#!/usr/bin/env julia

include(joinpath(dirname(@__DIR__), "fortran_parity_common.jl"))

const NSTEP = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2901
const RUN_DIR = length(ARGS) >= 2 ? abspath(ARGS[2]) : pwd()
const ORACLE_BEFORE = joinpath(RUN_DIR, "btran_oracle_before_step_n$(NSTEP).nc")
const ORACLE_AFTER = joinpath(RUN_DIR, "btran_oracle_after_hydrologydrainage_n$(NSTEP).nc")

function oracle_datetime(path)
    NCDataset(path) do ds
        ymd = Int(ds["timemgr_rst_curr_ymd"][])
        tod = Int(ds["timemgr_rst_curr_tod"][])
        DateTime(ymd ÷ 10000, (ymd ÷ 100) % 100, ymd % 100) + Second(tod)
    end
end

function inject_btran_oracle!(inst, bounds, _)
    NCDataset(ORACLE_BEFORE) do ds
        # The Fortran single-column bounds contain the three soil patches;
        # CLM.jl also materializes a fourth lake patch on column 2. Preserve
        # that unrelated patch and inject the shared soil subset explicitly.
        p = 1:size(ds["BTRAN"], 1)
        inst.canopystate.vegwp_patch[p, :] .= permutedims(ds["VEGWP"][:, :])
        inst.soilstate.smp_l_col[1, :] .= vec(ds["SMP"][:, 1])
        inst.soilstate.hk_l_col[1, :] .= vec(ds["HK"][:, 1])
        inst.soilstate.root_conductance_patch[p, :] .= permutedims(ds["ROOT_CONDUCTANCE"][:, :])
        inst.soilstate.k_soil_root_patch[p, :] .= permutedims(ds["K_SOIL_ROOT"][:, :])
        inst.soilstate.soil_conductance_patch[p, :] .= permutedims(ds["SOIL_CONDUCTANCE"][:, :])
        inst.soilstate.rootfr_patch[p, :] .= permutedims(ds["ROOTFR"][:, :])
        inst.soilstate.bsw_col[1, :] .= vec(ds["BSW"][:, 1])
        inst.soilstate.sucsat_col[1, :] .= vec(ds["SUCSAT"][:, 1])
        inst.atm2lnd.forc_pbot240_downscaled_patch[p] .= ds["PBOT240"][:]
        inst.atm2lnd.forc_pco2_240_patch[p] .= ds["PCO2_240"][:]
        inst.atm2lnd.forc_po2_240_patch[p] .= ds["PO2_240"][:]
        inst.canopystate.elai240_patch[p] .= ds["ELAI240"][:]
        inst.canopystate.laisun_patch[p] .= ds["LAISUN"][:]
        inst.canopystate.laisha_patch[p] .= ds["LAISHA"][:]
        inst.canopystate.elai_patch[p] .= ds["ELAI"][:]
        inst.canopystate.esai_patch[p] .= ds["ESAI"][:]
        inst.canopystate.tsai_patch[p] .= ds["TSAI"][:]
        inst.solarabs.par240d_z_patch[p, :] .= permutedims(ds["PAR240D_Z"][:, :])
        inst.solarabs.par240x_z_patch[p, :] .= permutedims(ds["PAR240X_Z"][:, :])
        inst.temperature.t_a10_patch[p] .= ds["T_A10"][:]
        inst.temperature.t_veg10_day_patch[p] .= ds["T_VEG10_DAY"][:]
        inst.temperature.t_veg10_night_patch[p] .= ds["T_VEG10_NIGHT"][:]
        inst.water.waterdiagnosticbulk_inst.rh10_af_patch[p] .= ds["RH10_AF"][:]
        inst.water.waterdiagnosticbulk_inst.fdry_patch[p] .= ds["FDRY"][:]
        inst.atm2lnd.forc_rho_downscaled_col[1] = ds["FORC_RHO"][1]
        inst.atm2lnd.forc_pbot_downscaled_col[1] = ds["FORC_PBOT"][1]
        inst.temperature.thm_patch[p] .= ds["THM"][:]
        inst.water.waterfluxbulk_inst.wf.qflx_tran_veg_patch[p] .= ds["QFLX_TRAN_VEG"][:]
        inst.frictionvel.rb10_patch[p] .= ds["RB10"][:]
        inst.photosyns.vcmx25_z_patch[p, :] .= permutedims(ds["VCMX25_Z"][:, :])
        inst.photosyns.jmx25_z_patch[p, :] .= permutedims(ds["JMX25_Z"][:, :])
        inst.photosyns.pnlc_z_patch[p, :] .= permutedims(ds["PNLC_Z"][:, :])
        inst.photosyns.enzs_z_patch[p, :] .= permutedims(ds["ENZS_Z"][:, :])
        inst.surfalb.nrad_patch[p] .= Int.(ds["NRAD"][:])
        inst.surfalb.tlai_z_patch[p, :] .= permutedims(ds["TLAI_Z"][:, :])
        inst.surfalb.fsun_z_patch[p, :] .= permutedims(ds["FSUN_Z"][:, :])
        inst.surfalb.fabd_sun_z_patch[p, :] .= permutedims(ds["FABD_SUN_Z"][:, :])
        inst.surfalb.fabd_sha_z_patch[p, :] .= permutedims(ds["FABD_SHA_Z"][:, :])
        inst.surfalb.fabi_sun_z_patch[p, :] .= permutedims(ds["FABI_SUN_Z"][:, :])
        inst.surfalb.fabi_sha_z_patch[p, :] .= permutedims(ds["FABI_SHA_Z"][:, :])
        inst.surfalb.vcmaxcintsun_patch[p] .= ds["VCMAXCINTSUN"][:]
        inst.surfalb.vcmaxcintsha_patch[p] .= ds["VCMAXCINTSHA"][:]
    end
    compare_inst_to_dump(inst, joinpath(RUN_DIR, "pdump_before_step_n$(NSTEP).nc");
                         label="injected initial state", tol=1e-9)
end

for path in (ORACLE_BEFORE, ORACLE_AFTER,
             joinpath(RUN_DIR, "pdump_before_step_n$(NSTEP).nc"),
             joinpath(RUN_DIR, "pdump_after_hydrologydrainage_n$(NSTEP).nc"))
    isfile(path) || error("missing shared-state oracle input: $path")
end

step_date = oracle_datetime(ORACLE_BEFORE)
println("BTRAN shared-state comparison: nstep=$NSTEP start=$step_date")
photo_trace = joinpath(RUN_DIR, "btran_photo_julia_n$(NSTEP).txt")
isfile(photo_trace) && rm(photo_trace)
have_photo_debug = isdefined(CLM, :PHS_PHOTO_DEBUG)
if have_photo_debug
    CLM.PHS_PHOTO_DEBUG[] = true
    CLM.PHS_PHOTO_DEBUG_PATH[] = photo_trace
end
(inst, bounds) = run_one_parity_step!(NSTEP; dumpdir=RUN_DIR, step_date=step_date,
                                      use_hydrstress=true, use_luna=true,
                                      forcing_offset_hours=-1,
                                      driver_is_first_step=false,
                                      pre_step_hook=inject_btran_oracle!)
if have_photo_debug
    CLM.PHS_PHOTO_DEBUG[] = false
    println("Julia PHS call trace: $photo_trace")
end

NCDataset(ORACLE_AFTER) do ds
    f = Float64.(ds["BTRAN"][:])
    j = Float64.(inst.energyflux.btran_patch[1:length(f)])
    println("patch,Julia,Fortran,abs_diff,rel_diff")
    for p in eachindex(f)
        a = abs(j[p] - f[p])
        r = a / (1 + max(abs(j[p]), abs(f[p])))
        @printf("%d,%.17g,%.17g,%.17g,%.17g\n", p, j[p], f[p], a, r)
    end
    checks = [
        ("K_SOIL_ROOT", inst.soilstate.k_soil_root_patch[1:length(f), :]),
        ("ROOT_CONDUCTANCE", inst.soilstate.root_conductance_patch[1:length(f), :]),
        ("SOIL_CONDUCTANCE", inst.soilstate.soil_conductance_patch[1:length(f), :]),
        ("VEGWP", inst.canopystate.vegwp_patch[1:length(f), :]),
        ("LAISUN", inst.canopystate.laisun_patch[1:length(f)]),
        ("LAISHA", inst.canopystate.laisha_patch[1:length(f)]),
        ("FDRY", inst.water.waterdiagnosticbulk_inst.fdry_patch[1:length(f)]),
    ]
    println("diagnostic,max_abs,max_rel")
    for (name, jraw) in checks
        fraw = Array(ds[name])
        fvals = vec(Float64.(fraw))
        jvals = ndims(fraw) == 2 ? vec(Float64.(permutedims(jraw))) : vec(Float64.(jraw))
        keep = [isfinite(fvals[i]) && abs(fvals[i]) < 1e35 for i in eachindex(fvals)]
        a = maximum(abs.(jvals[keep] .- fvals[keep]))
        r = maximum(abs.(jvals[keep] .- fvals[keep]) ./
                    (1 .+ max.(abs.(jvals[keep]), abs.(fvals[keep]))))
        @printf("%s,%.17g,%.17g\n", name, a, r)
    end
end

compare_inst_to_dump(inst, joinpath(RUN_DIR,
    "pdump_after_hydrologydrainage_n$(NSTEP).nc"); label="BTRAN shared-state", tol=1e-9)
