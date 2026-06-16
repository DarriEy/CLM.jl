# Probe: multi-step grass+tree availc decomposition vs Fortran dump availc.
# vegwp already verified <0.5% — this isolates the GPP/mr residual.
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const DD = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup"

function run_step(nstep)
    inst, bounds = run_one_parity_step!(nstep; use_cn=true, dumpdir=DD,
        use_hydrstress=true, use_luna=true,
        step_date=DateTime(2002,1,1)+Hour(nstep-1753153),
        forcing_file=replace(FFORCING, "clmforc.2003.nc"=>"clmforc.2002.nc"))
    ds = NCDataset(joinpath(DD, "pdump_after_hydrologydrainage_n$(nstep).nc"), "r")
    f_pfts = Int.(ds["pfts1d_itypveg"][:]); jl_pfts = Int.(inst.patch.itype)
    availc_f = haskey(ds,"availc") ? Float64.(ds["availc"][:]) : fill(NaN,length(f_pfts))
    ps = inst.photosyns; cf = inst.bgc_vegetation.cnveg_carbonflux_inst
    ef = inst.energyflux
    for pj in 1:length(jl_pfts)
        it = jl_pfts[pj]; it == 0 && continue
        pf = findfirst(==(it), f_pfts); pf === nothing && continue
        tag = it == 1 ? "TREE " : "GRASS"
        gpp_bd = cf.gpp_before_downreg_patch[pj]
        gpp_cp = cf.psnsun_to_cpool_patch[pj] + cf.psnshade_to_cpool_patch[pj]
        mr = cf.leaf_mr_patch[pj]+cf.froot_mr_patch[pj]+cf.livestem_mr_patch[pj]+cf.livecroot_mr_patch[pj]
        avc_j = cf.availc_patch[pj]; avc_f = availc_f[pf]
        @printf("n%d %s avc J=%.4e F=%.4e (%.3f)  btran=%.4f  gpp_bd=%.4e gpp_cp=%.4e mr=%.4e  psnsun=%.3f\n",
                nstep, tag, avc_j, avc_f, avc_j/avc_f, ef.btran_patch[pj],
                gpp_bd, gpp_cp, mr, ps.psnsun_patch[pj])
    end
    close(ds)
end
for n in (1757874,1757876,1757878,1757880,1757882)
    try; run_step(n); catch e; println("n$n FAIL: ", e); end
end
