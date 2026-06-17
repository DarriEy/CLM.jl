# Split the single-step mineral-N over-drain into NH4 vs NO3 per layer, Julia vs
# Fortran, to localize whether the 1.148x sminn over-drain is plant-uptake/nitrif
# (NH4 side) or leaching/denit (NO3 side). Re-injects Fortran state before step.
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const DD = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852
const NL = 6

inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=DD,
    use_hydrstress=true, use_luna=true,
    step_date=DateTime(2002,1,1)+Hour(NSTEP-1753153),
    forcing_file=replace(FFORCING, "clmforc.2003.nc"=>"clmforc.2002.nc"))

nf  = inst.soilbiogeochem_nitrogenflux
sns = inst.soilbiogeochem_nitrogenstate
bds = NCDataset(joinpath(DD,"pdump_before_step_n$(NSTEP).nc"),"r")
ads = NCDataset(joinpath(DD,"pdump_after_hydrologydrainage_n$(NSTEP).nc"),"r")

f_nh4_b = Float64.(bds["smin_nh4_vr"][:,1]); f_nh4_a = Float64.(ads["smin_nh4_vr"][:,1])
f_no3_b = Float64.(bds["smin_no3_vr"][:,1]); f_no3_a = Float64.(ads["smin_no3_vr"][:,1])
j_nh4 = sns.smin_nh4_vr_col; j_no3 = sns.smin_no3_vr_col

@printf("%-3s | %11s %11s %11s | %11s %11s %11s\n", "lev",
        "F dNH4","J dNH4","d(J-F)NH4", "F dNO3","J dNO3","d(J-F)NO3")
for j in 1:NL
    dFh=f_nh4_a[j]-f_nh4_b[j]; dJh=Float64(j_nh4[1,j])-f_nh4_b[j]
    dFn=f_no3_a[j]-f_no3_b[j]; dJn=Float64(j_no3[1,j])-f_no3_b[j]
    @printf("%-3d | %11.3e %11.3e %11.3e | %11.3e %11.3e %11.3e\n",
            j, dFh,dJh,dJh-dFh, dFn,dJn,dJn-dFn)
end
sFh=sum(f_nh4_a[1:NL].-f_nh4_b[1:NL]); sJh=sum(Float64(j_nh4[1,j])-f_nh4_b[j] for j in 1:NL)
sFn=sum(f_no3_a[1:NL].-f_no3_b[1:NL]); sJn=sum(Float64(j_no3[1,j])-f_no3_b[j] for j in 1:NL)
@printf("--- top-%d  dNH4 F=%.3e J=%.3e ratio=%.3f | dNO3 F=%.3e J=%.3e ratio=%.3f\n",
        NL, sFh,sJh, sJh/sFh, sFn,sJn, sJn/sFn)

# Julia per-layer sink/source terms (gN/m3/s) to attribute the drain
@printf("\n%-3s | %11s %11s %11s %11s %11s\n", "lev",
        "gross_nmin","immob_nh4","f_nit","f_denit","nh4_to_plt")
for j in 1:NL
    @printf("%-3d | %11.3e %11.3e %11.3e %11.3e %11.3e\n", j,
        Float64(nf.gross_nmin_vr_col[1,j]), Float64(nf.actual_immob_nh4_vr_col[1,j]),
        Float64(nf.f_nit_vr_col[1,j]), Float64(nf.f_denit_vr_col[1,j]),
        Float64(nf.smin_nh4_to_plant_vr_col[1,j]))
end
# leaching/runoff if present
for nm in ("smin_no3_leached_vr_col","smin_no3_runoff_vr_col","sminn_leached_vr_col")
    if hasproperty(nf, Symbol(nm))
        v = getfield(nf, Symbol(nm))
        @printf("%-20s top-%d sum = %.3e\n", nm, NL, sum(Float64(v[1,j]) for j in 1:NL))
    end
end
close(bds); close(ads)
