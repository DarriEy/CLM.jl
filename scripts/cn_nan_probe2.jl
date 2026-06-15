# Pinpoint the photosynthesis NaN input for the veg patches in the use_cn step.
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const BGC_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_dumps_2202"
const NSTEP = 1753153

pv(name, a) = (A=Array(a); @printf("  %-26s %s\n", name, string(round.(Float64.(vec(A)), sigdigits=4))))
inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=BGC_DUMPDIR, step_date=DateTime(2003,1,1))
cs = inst.canopystate; sa = inst.surfalb; so = inst.solarabs; ph = inst.photosyns
println("\n=== sun geometry / radiation ===")
pv("coszen_grc", sa.coszen_grc); pv("coszen_col", sa.coszen_col)
println("=== canopy state (per patch: bare,tree,grass,bare2) ===")
pv("elai", cs.elai_patch); pv("esai", cs.esai_patch)
pv("laisun", cs.laisun_patch); pv("laisha", cs.laisha_patch); pv("fsun", cs.fsun_patch)
println("=== absorbed radiation ===")
pv("sabv", so.sabv_patch); pv("sabg", so.sabg_patch)
hasproperty(so,:parsun_z_patch) && @printf("  parsun_z nan=%d size=%s\n", count(isnan,so.parsun_z_patch), string(size(so.parsun_z_patch)))
hasproperty(so,:parsun_z_patch) && pv("parsun_z[:,1]", so.parsun_z_patch[:,1])
println("=== photosynthesis outputs/inputs ===")
for f in (:psnsun_patch,:fpsn_patch,:rssun_patch,:lnca_patch,:vcmaxcint_patch)
    hasproperty(ph,f) && pv("ph.$f", getfield(ph,f))
end
println("=== veg leaf state feeding vcmax ===")
ccs=inst.bgc_vegetation.cnveg_carbonstate_inst; cns=inst.bgc_vegetation.cnveg_nitrogenstate_inst
pv("leafc", ccs.leafc_patch); pv("leafn", cns.leafn_patch)
pv("t_veg", inst.temperature.t_veg_patch); pv("t_a10", inst.temperature.t_a10_patch)
hasproperty(inst.gridcell,:dayl_grc) && pv("dayl_grc", inst.gridcell.dayl_grc)
