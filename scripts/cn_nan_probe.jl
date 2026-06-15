# Localize the post-step CN NaN: run one use_cn step, then scan the CN flux
# chain (photosynthesis -> allocation -> decomposition) for NaN counts.
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const BGC_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_dumps_2202"
const NSTEP = 1753153

nanc(a) = count(isnan, Array(a))
function rep(name, a)
    A = Array(a); @printf("  %-34s nan=%d/%d  finite-sample=%s\n", name, count(isnan,A), length(A),
        string(round.(filter(isfinite, vec(A))[1:min(end,3)], sigdigits=4)))
end

inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=BGC_DUMPDIR, step_date=DateTime(2003,1,1))
ph  = inst.photosyns
cf  = inst.bgc_vegetation.cnveg_carbonflux_inst
cs  = inst.bgc_vegetation.cnveg_carbonstate_inst
scf = inst.soilbiogeochem_carbonflux
println("\n=== photosynthesis outputs ===")
for f in (:psnsun_patch,:psnsha_patch,:fpsn_patch,:rssun_patch,:rssha_patch)
    hasproperty(ph,f) && rep("photosyns.$f", getfield(ph,f))
end
println("=== cnveg carbon FLUX (veg input chain) ===")
for f in (:gpp_patch,:availc_patch,:psnsun_to_cpool_patch,:psnshade_to_cpool_patch,
          :leaf_mr_patch,:leafc_to_litter_patch,:cpool_to_leafc_patch)
    hasproperty(cf,f) && rep("cnvegcf.$f", getfield(cf,f))
end
println("=== cnveg carbon STATE (post-step) ===")
for f in (:leafc_patch,:frootc_patch,:deadstemc_patch,:cpool_patch,:xsmrpool_patch)
    hasproperty(cs,f) && rep("cnvegcs.$f", getfield(cs,f))
end
println("=== cnveg STATE (allometry / demand) ===")
cst = inst.bgc_vegetation.cnveg_state_inst
for f in (:c_allometry_patch,:n_allometry_patch)
    hasproperty(cst,f) && rep("cnvegst.$f", getfield(cst,f))
end
cnf = inst.bgc_vegetation.cnveg_nitrogenflux_inst
for f in (:plant_ndemand_patch,)
    hasproperty(cnf,f) && rep("cnvegnf.$f", getfield(cnf,f))
end
rep("soilbgc_state.fpg_col", inst.soilbiogeochem_state.fpg_col)
println("=== soil decomposition FLUX ===")
for f in (:phr_vr_col,:hr_vr_col,:decomp_cascade_hr_vr_col,:decomp_cpools_sourcesink_col)
    hasproperty(scf,f) && rep("soilcf.$f", getfield(scf,f))
end
println("=== soil decomposition STATE (post-step) ===")
rep("decomp_cpools_vr_col", inst.soilbiogeochem_carbonstate.decomp_cpools_vr_col)
rep("sminn_vr_col", inst.soilbiogeochem_nitrogenstate.sminn_vr_col)
