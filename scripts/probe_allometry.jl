include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const DD = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852
inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=DD,
    use_hydrstress=true, use_luna=true,
    step_date=DateTime(2002,1,1)+Hour(NSTEP-1753153),
    forcing_file=replace(FFORCING, "clmforc.2003.nc"=>"clmforc.2002.nc"))
cnf=inst.bgc_vegetation.cnveg_nitrogenflux_inst
ccf=inst.bgc_vegetation.cnveg_carbonflux_inst
cst=inst.bgc_vegetation.cnveg_state_inst
bds=NCDataset(joinpath(DD,"pdump_before_step_n$(NSTEP).nc"),"r")
fnd=Float64.(bds["plant_ndemand"][:]); fca=Float64.(bds["plant_calloc"][:])
println("p | availc_J  nallom_J  callom_J  ratioJ=na/ca | ndemand_recon  ndemand_field | F_calloc  F_ndemand  F_ratio")
for p in 2:3
  av=Float64(ccf.availc_patch[p]); na=Float64(cst.n_allometry_patch[p]); ca=Float64(cst.c_allometry_patch[p])
  rj=na/ca; recon=av*rj; fld=Float64(cnf.plant_ndemand_patch[p])
  frat=fnd[p]/fca[p]
  @printf("%d | %.4e %.4e %.4e %.5f | %.4e %.4e | %.4e %.4e %.5f\n",
          p, av, na, ca, rj, recon, fld, fca[p], fnd[p], frat)
  @printf("   ratios vs Fortran: availc(J/Fcalloc)=%.3f  allom_ratio(J/F)=%.3f  ndemand(J/F)=%.3f\n",
          av/fca[p], rj/frat, fld/fnd[p])
end
close(bds)
