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
        inst.atm2lnd.forc_pbot240_downscaled_patch[p] .= ds["PBOT240"][:]
        inst.atm2lnd.forc_pco2_240_patch[p] .= ds["PCO2_240"][:]
        inst.atm2lnd.forc_po2_240_patch[p] .= ds["PO2_240"][:]
        inst.canopystate.elai240_patch[p] .= ds["ELAI240"][:]
        inst.solarabs.par240d_z_patch[p, :] .= permutedims(ds["PAR240D_Z"][:, :])
        inst.solarabs.par240x_z_patch[p, :] .= permutedims(ds["PAR240X_Z"][:, :])
        inst.temperature.t_a10_patch[p] .= ds["T_A10"][:]
        inst.temperature.t_veg10_day_patch[p] .= ds["T_VEG10_DAY"][:]
        inst.temperature.t_veg10_night_patch[p] .= ds["T_VEG10_NIGHT"][:]
        inst.water.waterdiagnosticbulk_inst.rh10_af_patch[p] .= ds["RH10_AF"][:]
        inst.frictionvel.rb10_patch[p] .= ds["RB10"][:]
        inst.photosyns.vcmx25_z_patch[p, :] .= permutedims(ds["VCMX25_Z"][:, :])
        inst.photosyns.jmx25_z_patch[p, :] .= permutedims(ds["JMX25_Z"][:, :])
        inst.photosyns.pnlc_z_patch[p, :] .= permutedims(ds["PNLC_Z"][:, :])
        inst.photosyns.enzs_z_patch[p, :] .= permutedims(ds["ENZS_Z"][:, :])
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
(inst, bounds) = run_one_parity_step!(NSTEP; dumpdir=RUN_DIR, step_date=step_date,
                                      use_hydrstress=true, use_luna=true,
                                      forcing_offset_hours=-1,
                                      pre_step_hook=inject_btran_oracle!)

NCDataset(ORACLE_AFTER) do ds
    f = Float64.(ds["BTRAN"][:])
    j = Float64.(inst.energyflux.btran_patch[1:length(f)])
    println("patch,Julia,Fortran,abs_diff,rel_diff")
    for p in eachindex(f)
        a = abs(j[p] - f[p])
        r = a / (1 + max(abs(j[p]), abs(f[p])))
        @printf("%d,%.17g,%.17g,%.17g,%.17g\n", p, j[p], f[p], a, r)
    end
end

compare_inst_to_dump(inst, joinpath(RUN_DIR,
    "pdump_after_hydrologydrainage_n$(NSTEP).nc"); label="BTRAN shared-state", tol=1e-9)
