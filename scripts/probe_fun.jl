# Probe: single-step FUN N-flux parity vs the Fortran after_competition dump.
# FUN computes plant N acquisition (active/mycorrhizal/nonmyc uptake, fixation,
# retranslocation) — the smin-N-balance lever behind the multi-step drift.
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const DD = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852

inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=DD,
    use_hydrstress=true, use_luna=true,
    step_date=DateTime(2002,1,1)+Hour(NSTEP-1753153),
    forcing_file=replace(FFORCING, "clmforc.2003.nc"=>"clmforc.2002.nc"))

nf = inst.bgc_vegetation.cnveg_nitrogenflux_inst
ds = NCDataset(joinpath(DD,"pdump_after_competition_n$(NSTEP).nc"),"r")
f_pfts = Int.(ds["pfts1d_itypveg"][:]); jl_pfts = Int.(inst.patch.itype)
fv(v,i) = haskey(ds,v) ? (a=ds[v][:]; i<=length(a) ? Float64(a[i]) : NaN) : NaN

# (Fortran dump var, Julia field) pairs
pairs = [("plant_ndemand", :plant_ndemand_patch), ("Nuptake", :Nuptake_patch),
         ("Nfix", :Nfix_patch), ("Nretrans", :Nretrans_patch),
         ("Nactive_nh4", :Nactive_nh4_patch), ("Necm_nh4", :Necm_nh4_patch),
         ("Nnonmyc_nh4", :Nnonmyc_nh4_patch), ("Nnonmyc_no3", :Nnonmyc_no3_patch)]

for pj in 1:length(jl_pfts)
    it = jl_pfts[pj]; it == 0 && continue
    pf = findfirst(==(it), f_pfts); pf === nothing && continue
    tag = it == 1 ? "TREE " : "GRASS"
    @printf("\n%s  (Julia | Fortran | ratio)\n", tag)
    for (fvar, jfield) in pairs
        jv = hasproperty(nf, jfield) ? Float64(getfield(nf, jfield)[pj]) : NaN
        fval = fv(fvar, pf)
        r = (abs(fval) > 1e-30) ? jv/fval : NaN
        @printf("  %-14s %12.4e | %12.4e | %6.3f%s\n", fvar, jv, fval, r,
                abs(r-1) > 0.02 ? "  DIFF" : "")
    end
end
close(ds)

# leafc_to_litter_fun (drives retranslocation) — Julia vs Fortran
let inst2=inst
    cf = inst2.bgc_vegetation.cnveg_carbonflux_inst
    ds2 = NCDataset(joinpath(DD,"pdump_after_competition_n$(NSTEP).nc"),"r")
    f2 = Int.(ds2["pfts1d_itypveg"][:]); j2 = Int.(inst2.patch.itype)
    println("\n--- leafc_to_litter_fun (drives retrans) ---")
    for pj in 1:length(j2)
        it=j2[pj]; it==0 && continue; pf=findfirst(==(it),f2); pf===nothing && continue
        fval = haskey(ds2,"leafc_to_litter_fun") ? Float64(ds2["leafc_to_litter_fun"][pf]) : NaN
        @printf("  %-5s leafc_to_litter_fun J=%.4e F=%.4e\n", it==1 ? "TREE" : "PFT$it",
                cf.leafc_to_litter_fun_patch[pj], fval)
    end
    close(ds2)
end
