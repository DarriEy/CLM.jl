# Clean single-step translation error: Julia evolved after-state vs Fortran AFTER
# dump, ABSOLUTE per-layer relative error (not deltas vs stale before-dump).
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const DD = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852
inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=DD,
    use_hydrstress=true, use_luna=true,
    step_date=DateTime(2002,1,1)+Hour(NSTEP-1753153),
    forcing_file=replace(FFORCING, "clmforc.2003.nc"=>"clmforc.2002.nc"))
sns=inst.soilbiogeochem_nitrogenstate
ad=NCDataset(joinpath(DD,"pdump_after_hydrologydrainage_n$(NSTEP).nc"),"r")
fnh4=Float64.(ad["smin_nh4_vr"][:,1]); fno3=Float64.(ad["smin_no3_vr"][:,1]); fsm=Float64.(ad["sminn_vr"][:,1])
jnh4=sns.smin_nh4_vr_col; jno3=sns.smin_no3_vr_col; jsm=sns.sminn_vr_col
@printf("%-3s | %12s %12s %8s | %12s %12s %8s\n","lev","J nh4","F nh4","relerr","J sminn","F sminn","relerr")
worst=0.0
for j in 1:8
  rn=abs(Float64(jnh4[1,j])-fnh4[j])/max(abs(fnh4[j]),1e-12)
  rs=abs(Float64(jsm[1,j])-fsm[j])/max(abs(fsm[j]),1e-12)
  global worst=max(worst,rs)
  @printf("%-3d | %12.5e %12.5e %7.3f%% | %12.5e %12.5e %7.3f%%\n",j,
          Float64(jnh4[1,j]),fnh4[j],100rn,Float64(jsm[1,j]),fsm[j],100rs)
end
@printf("--- worst sminn rel err (top 8 layers) = %.3f%%\n",100worst)
close(ad)
