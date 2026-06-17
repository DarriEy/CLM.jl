# Scope the sminn drift: per-step sminn_vr change (Julia vs Fortran before/after
# dumps) + Julia soil-BGC N-flux terms, to localize which source term drifts.
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const DD = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852
const NL = 5   # show top 5 decomp layers

inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=DD,
    use_hydrstress=true, use_luna=true,
    step_date=DateTime(2002,1,1)+Hour(NSTEP-1753153),
    forcing_file=replace(FFORCING, "clmforc.2003.nc"=>"clmforc.2002.nc"))

nf = inst.soilbiogeochem_nitrogenflux
sns = inst.soilbiogeochem_nitrogenstate
# Fortran sminn_vr before and after the step
bds = NCDataset(joinpath(DD,"pdump_before_step_n$(NSTEP).nc"),"r")
ads = NCDataset(joinpath(DD,"pdump_after_hydrologydrainage_n$(NSTEP).nc"),"r")
fb(v) = haskey(bds,v) ? Float64.(bds[v][:,1]) : fill(NaN,NL)
fa(v) = haskey(ads,v) ? Float64.(ads[v][:,1]) : fill(NaN,NL)

sminn_b = fb("sminn_vr"); sminn_a = fa("sminn_vr")
jl_sminn = sns.sminn_vr_col

@printf("%-4s | %12s %12s %12s | %12s %12s %12s\n", "lev",
        "F Δsminn", "J Δsminn", "Δ(J-F)", "J gross_nmin", "J immob", "J f_denit")
for j in 1:NL
    dF = sminn_a[j] - sminn_b[j]
    dJ = Float64(jl_sminn[1,j]) - sminn_b[j]   # Julia injected = before, so after-before
    @printf("%-4d | %12.4e %12.4e %12.4e | %12.4e %12.4e %12.4e\n", j,
            dF, dJ, dJ-dF,
            Float64(nf.gross_nmin_vr_col[1,j]), Float64(nf.actual_immob_vr_col[1,j]),
            Float64(nf.f_denit_vr_col[1,j]))
end
# layer sums
sF = sum(sminn_a[1:NL].-sminn_b[1:NL]); sJ = sum(Float64(jl_sminn[1,j])-sminn_b[j] for j in 1:NL)
@printf("--- top-%d Δsminn: Fortran=%.4e Julia=%.4e ratio=%.3f\n", NL, sF, sJ, sJ/sF)
@printf("--- Julia flux sums (gN/m3/s) gross_nmin=%.3e immob=%.3e f_nit=%.3e f_denit=%.3e sminn_to_plant=%.3e\n",
        sum(Float64(nf.gross_nmin_vr_col[1,j]) for j in 1:NL),
        sum(Float64(nf.actual_immob_vr_col[1,j]) for j in 1:NL),
        sum(Float64(nf.f_nit_vr_col[1,j]) for j in 1:NL),
        sum(Float64(nf.f_denit_vr_col[1,j]) for j in 1:NL),
        sum(Float64(nf.sminn_to_plant_vr_col[1,j]) for j in 1:NL))
close(bds); close(ads)
