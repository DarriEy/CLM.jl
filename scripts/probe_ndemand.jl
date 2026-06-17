include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const DD = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852
inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=DD,
    use_hydrstress=true, use_luna=true,
    step_date=DateTime(2002,1,1)+Hour(NSTEP-1753153),
    forcing_file=replace(FFORCING, "clmforc.2003.nc"=>"clmforc.2002.nc"))
cnf = inst.bgc_vegetation.cnveg_nitrogenflux_inst
ccf = inst.bgc_vegetation.cnveg_carbonflux_inst
cst = inst.bgc_vegetation.cnveg_state_inst
bds = NCDataset(joinpath(DD,"pdump_before_step_n$(NSTEP).nc"),"r")
g(v)=haskey(bds,v) ? Float64.(bds[v][:]) : fill(NaN,3)
fnd=g("plant_ndemand"); fca=g("plant_calloc"); fcn=g("plantCN")
@printf("%-22s | %4s %14s %14s %7s\n","field","p","Julia","Fortran","ratio")
for p in 2:3
  @printf("%-22s | %4d %14.5e %14.5e %7.3f\n","plant_ndemand",p,Float64(cnf.plant_ndemand_patch[p]),fnd[p],Float64(cnf.plant_ndemand_patch[p])/fnd[p])
  @printf("%-22s | %4d %14.5e %14.5e %7.3f\n","availc(plant_calloc?)",p,Float64(ccf.availc_patch[p]),fca[p],Float64(ccf.availc_patch[p])/fca[p])
  @printf("%-22s | %4d %14.5e %14.5e %7.3f\n","plantCN",p,Float64(cst.plantCN_patch[p]),fcn[p],Float64(cst.plantCN_patch[p])/fcn[p])
end
close(bds)
