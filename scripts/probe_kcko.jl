include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const SC = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup"
inst,bounds = run_one_parity_step!(1757873; use_cn=true, dumpdir=SC, use_hydrstress=true, use_luna=true,
  step_date=DateTime(2002,1,1)+Hour(1757873-1753153),
  forcing_file=replace(FFORCING,"clmforc.2003.nc"=>"clmforc.2002.nc"))
ps=inst.photosyns; a=inst.atm2lnd
for p in 2:3
  @printf("p%d  kc=%.5f ko=%.5f cp=%.5f oair=%.4f ci_sun=%.5f\n", p,
    Float64(ps.kc_patch[p]), Float64(ps.ko_patch[p]), Float64(ps.cp_patch[p]),
    Float64(a.forc_po2_grc[1]), Float64(ps.cisun_z_patch[p,1]))
end
