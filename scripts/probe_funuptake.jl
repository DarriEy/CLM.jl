# Is the FUN NH4 over-drain a TOTAL-uptake error or a vertical-DISTRIBUTION error?
# Compare (a) Julia per-patch FUN N-uptake stream totals vs Fortran dump
# (Nactive/Nam/Necm/Nnonmyc/Nfix _nh4) and (b) full 25-layer column NH4 drain
# (Julia FUN nh4 uptake vs Fortran d smin_nh4). Re-injects Fortran state.
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const DD = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852

inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=DD,
    use_hydrstress=true, use_luna=true,
    step_date=DateTime(2002,1,1)+Hour(NSTEP-1753153),
    forcing_file=replace(FFORCING, "clmforc.2003.nc"=>"clmforc.2002.nc"))

cnf = inst.bgc_vegetation.cnveg_nitrogenflux_inst
snf = inst.soilbiogeochem_nitrogenflux
bds = NCDataset(joinpath(DD,"pdump_after_hydrologydrainage_n$(NSTEP).nc"),"r")
ads = NCDataset(joinpath(DD,"pdump_after_hydrologydrainage_n$(NSTEP).nc"),"r")
nl = CLM.varpar.nlevdecomp
dz = CLM.dzsoi_decomp[]

println("=== Per-patch FUN N-uptake stream totals (gN/m2/s): Julia vs Fortran ===")
@printf("%-18s | %4s %12s %12s %7s\n","stream","p","Julia","Fortran","ratio")
for (jl_field, f_var) in (("Nactive_nh4_patch","Nactive_nh4"),
                          ("Nnonmyc_nh4_patch","Nnonmyc_nh4"),
                          ("Nactive_no3_patch","Nactive_no3"),
                          ("Nnonmyc_no3_patch","Nnonmyc_no3"),
                          ("Nfix_patch","Nfix"))
    jv = hasproperty(cnf, Symbol(jl_field)) ? getfield(cnf, Symbol(jl_field)) : nothing
    fv = haskey(bds, f_var) ? Float64.(bds[f_var][:]) : nothing
    for p in 2:3   # patches 2,3 = tree,grass (slot 1 bare)
        jval = jv === nothing ? NaN : Float64(jv[p])
        fval = fv === nothing ? NaN : fv[p]
        @printf("%-18s | %4d %12.4e %12.4e %7.3f\n", jl_field, p, jval, fval, jval/fval)
    end
end
# Total NH4 vs NO3 vs (NH4+NO3) per patch — partitioning vs magnitude
println("\n=== Per-patch TOTAL uptake: NH4, NO3, sum (Julia vs Fortran) ===")
for p in 2:3
    jnh4 = Float64(cnf.Nactive_nh4_patch[p]) + Float64(cnf.Nnonmyc_nh4_patch[p])
    jno3 = Float64(cnf.Nactive_no3_patch[p]) + Float64(cnf.Nnonmyc_no3_patch[p])
    fnh4 = Float64(bds["Nactive_nh4"][p]) + Float64(bds["Nnonmyc_nh4"][p])
    fno3 = Float64(bds["Nactive_no3"][p]) + Float64(bds["Nnonmyc_no3"][p])
    @printf("p%d NH4: J=%.4e F=%.4e (%.3f) | NO3: J=%.4e F=%.4e (%.3f) | SUM: J=%.4e F=%.4e (%.3f)\n",
            p, jnh4, fnh4, jnh4/fnh4, jno3, fno3, jno3/fno3, jnh4+jno3, fnh4+fno3, (jnh4+jno3)/(fnh4+fno3))
end

println("\n=== Full-column NH4 drain (all $nl layers): Julia FUN-uptake vs Fortran d smin_nh4 ===")
f_nh4_b = Float64.(bds["smin_nh4_vr"][:,1]); f_nh4_a = Float64.(ads["smin_nh4_vr"][:,1])
# Julia FUN nh4 uptake col, gN/m3/s -> *dz*dt for gN/m2 removed; but compare as d(smin_nh4) gN/m3
jl_fun_nh4 = snf.sminn_to_plant_fun_nh4_vr_col  # gN/m3/s
dt = 3600.0
tot_jl_uptake = sum(Float64(jl_fun_nh4[1,j])*dz[j]*dt for j in 1:nl)   # gN/m2 removed by FUN
tot_F_drain   = -sum((f_nh4_a[j]-f_nh4_b[j])*dz[j] for j in 1:nl)      # gN/m2 net NH4 lost (Fortran)
@printf("Total FUN NH4 uptake (Julia, gN/m2)   = %.5e\n", tot_jl_uptake)
@printf("Total NH4 net drain (Fortran, gN/m2)  = %.5e  ratio(up/drain)=%.3f\n",
        tot_F_drain, tot_jl_uptake/tot_F_drain)
# vertical profile of uptake fraction — BOTH dz-weighted (mass per layer)
println("\nlev | dz |  Julia FUN-nh4 massfrac | Fortran dNH4 massfrac")
jl_frac_tot = sum(Float64(jl_fun_nh4[1,j])*dz[j] for j in 1:nl)
F_frac_tot  = -sum((f_nh4_a[j]-f_nh4_b[j])*dz[j] for j in 1:nl)
for j in 1:10
    jf = Float64(jl_fun_nh4[1,j])*dz[j]/jl_frac_tot
    ff = -(f_nh4_a[j]-f_nh4_b[j])*dz[j]/F_frac_tot
    @printf("%-3d | %.3f | %22.4f | %.4f\n", j, dz[j], jf, ff)
end
close(bds); close(ads)
