# Probe: decompose grass vs tree froot_mr inputs (frootn, crootfr profile, tcsoi)
# against the Fortran dump. froot_mr = frootn*br_root_acc*sum_j(tcsoi[j]*crootfr[j]).
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const DD = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852

inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=DD,
    use_hydrstress=true, use_luna=true,
    step_date=DateTime(2002,1,1)+Hour(NSTEP-1753153),
    forcing_file=replace(FFORCING, "clmforc.2003.nc"=>"clmforc.2002.nc"))

ds = NCDataset(joinpath(DD,"pdump_after_hydrologydrainage_n$(NSTEP).nc"),"r")
f_pfts = Int.(ds["pfts1d_itypveg"][:]); jl_pfts = Int.(inst.patch.itype)
frootn_f = Float64.(ds["frootn"][:])
cns = inst.bgc_vegetation.cnveg_nitrogenstate_inst
ss = inst.soilstate; cf = inst.bgc_vegetation.cnveg_carbonflux_inst
patch = inst.patch
nlevgrnd = size(ss.crootfr_patch, 2)

Q10 = inst.cn_shared_params.Q10
tsno = inst.temperature.t_soisno_col
nlevsno = CLM.varpar.nlevsno
TFRZ = CLM.TFRZ
for pj in 1:length(jl_pfts)
    it = jl_pfts[pj]; it == 0 && continue
    pf = findfirst(==(it), f_pfts); pf === nothing && continue
    tag = it == 1 ? "TREE " : "GRASS"
    c = patch.column[pj]
    fn_j = cns.frootn_patch[pj]; fn_f = frootn_f[pf]
    cr = ss.crootfr_patch[pj, :]
    tcsoi = [Q10^((tsno[c, nlevsno+j]-TFRZ-20.0)/10.0) for j in 1:nlevgrnd]
    wsum = sum(tcsoi[j]*cr[j] for j in 1:nlevgrnd)
    @printf("%s col=%d frootn J=%.5f F=%.5f (%.3f) sum(cr)=%.4f sum(tcsoi*cr)=%.4f froot_mr=%.4e\n",
            tag, c, fn_j, fn_f, fn_j/fn_f, sum(cr[1:nlevgrnd]), wsum, cf.froot_mr_patch[pj])
    @printf("   crootfr[1:8]= %s\n", join([@sprintf("%.4f",cr[j]) for j in 1:8], " "))
    @printf("   tcsoi[1:8]  = %s\n", join([@sprintf("%.4f",tcsoi[j]) for j in 1:8], " "))
    @printf("   tsoil[1:8]K = %s\n", join([@sprintf("%.2f",tsno[c,nlevsno+j]) for j in 1:8], " "))
end
close(ds)
