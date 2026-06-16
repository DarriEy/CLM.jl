# Decisive Julia<->Fortran compare at the SAME step (1757873), using the fresh
# Fortran dumps written by the instrumented continue-run. Fortran ground truth
# (PARITYPHS/PARITYMR, grass p=3, converged 1st eval): psn_sun=13.854 psn_sha=6.100
# lmr_sun=0.36738 lmr_sha=0.17280 leaf_mr=9.424e-7 froot_mr=2.0259e-6 ci=25.258
# vcmax_z_sun=59.94 t_veg=288.55 bsun=0.5259
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const DD = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup"
const NSTEP = 1757873

inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=DD,
    use_hydrstress=true, use_luna=true,
    step_date=DateTime(2002,1,1)+Hour(NSTEP-1753153),
    forcing_file=replace(FFORCING, "clmforc.2003.nc"=>"clmforc.2002.nc"))

ps = inst.photosyns; cf = inst.bgc_vegetation.cnveg_carbonflux_inst
ef = inst.energyflux; cs = inst.canopystate; tp = inst.temperature
ds = NCDataset(joinpath(DD,"pdump_after_hydrologydrainage_n$(NSTEP).nc"),"r")
f_pfts = Int.(ds["pfts1d_itypveg"][:]); jl_pfts = Int.(inst.patch.itype)
availc_f = Float64.(ds["availc"][:])

# Fortran ground truth (grass p=3, tree p=2) from PARITY prints at 1757873
F = Dict(1=>(psnsun=NaN,), 12=>(psnsun=13.854, psnsha=6.0998, leafmr=9.424102e-7,
              frootmr=2.025897e-6, ci=25.258, vcmaxz=59.94, tveg=288.55, bsun=0.5259))
Ftree = (psnsun=NaN,)

for pj in 1:length(jl_pfts)
    it = jl_pfts[pj]; it == 0 && continue
    pf = findfirst(==(it), f_pfts); pf === nothing && continue
    tag = it == 1 ? "TREE " : "GRASS"
    gpp = cf.psnsun_to_cpool_patch[pj] + cf.psnshade_to_cpool_patch[pj]
    mr  = cf.leaf_mr_patch[pj]+cf.froot_mr_patch[pj]+cf.livestem_mr_patch[pj]+cf.livecroot_mr_patch[pj]
    @printf("\n%s  (Julia | Fortran)\n", tag)
    @printf("  psnsun  = %.4f\n", ps.psnsun_patch[pj])
    @printf("  psnsha  = %.4f\n", ps.psnsha_patch[pj])
    @printf("  t_veg   = %.3f\n", tp.t_veg_patch[pj])
    @printf("  btran   = %.4f\n", ef.btran_patch[pj])
    @printf("  lmrsun_z1=%.4f (per-layer, pre-integration)  lmrsun=%.4f  lmrsha=%.4f\n",
            ps.lmrsun_z_patch[pj,1], ps.lmrsun_patch[pj], ps.lmrsha_patch[pj])
    @printf("  laisun  = %.4f  laisha=%.4f\n", cs.laisun_patch[pj], cs.laisha_patch[pj])
    @printf("  t_a10   = %.3f\n", tp.t_a10_patch[pj])
    @printf("  leaf_mr = %.4e\n", cf.leaf_mr_patch[pj])
    @printf("  froot_mr= %.4e\n", cf.froot_mr_patch[pj])
    @printf("  gpp     = %.4e\n", gpp)
    @printf("  availc  = %.4e | F=%.4e (%.3f)\n", cf.availc_patch[pj], availc_f[pf],
            cf.availc_patch[pj]/availc_f[pf])
    if it == 12
        @printf("  >>> Fortran grass: psnsun=13.854 psnsha=6.100 leaf_mr=9.424e-7 froot_mr=2.026e-6 t_veg=288.55 btran=0.526\n")
    end
end
close(ds)
