# Grass-GPP polish: compare Julia psn/ci/vcmax at nstep 1757873 against the
# Fortran PARITYPHS converged print. The injectable IC is the 1757872 end-state
# (after_hydrologydrainage), aliased as before_step_n1757873 in /tmp/psn_scratch.
# One step from an exact IC → drift contamination ~0.024%, negligible vs the ~1.2%
# grass residual we are chasing.
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
const SC = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup"  # clean before_step dumps (psn_probe run)
const NSTEP = 1757873
inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=SC,
    use_hydrstress=true, use_luna=true,
    step_date=DateTime(2002,1,1)+Hour(NSTEP-1753153),
    forcing_file=replace(FFORCING, "clmforc.2003.nc"=>"clmforc.2002.nc"))
ps = inst.photosyns; ef = inst.energyflux; tv = inst.temperature
# Fortran converged (last iter) at nstep 1757873
F = Dict(
 2 => (psnsun=8.28403, psnsha=0.794426, vcmax=53.1953, ci=24.8585, bsun=0.672235, tveg=299.446),
 3 => (psnsun=13.7486, psnsha=6.63364, vcmax=91.1934, ci=25.0443, bsun=0.421422, tveg=292.908))
@printf("%-4s %-8s | %12s %12s %7s\n","p","field","Julia","Fortran","ratio")
for p in 2:3
  f = F[p]
  J = (psnsun=Float64(ps.psnsun_patch[p]), psnsha=Float64(ps.psnsha_patch[p]),
       vcmax=Float64(ps.vcmax_z_patch[p,1]), ci=Float64(ps.cisun_z_patch[p,1]),
       bsun=Float64(ef.bsun_patch[p]), tveg=Float64(tv.t_veg_patch[p]))
  for k in (:psnsun,:psnsha,:vcmax,:ci,:bsun,:tveg)
    jv=getfield(J,k); fv=getfield(f,k)
    @printf("p%-3d %-8s | %12.5e %12.5e %7.4f\n", p, k, jv, fv, jv/fv)
  end
  println("-"^48)
end
# Deeper: PHS sunlit vcmax (vcmax_z_phs_patch[p,sun=1,1]) vs Fortran vcmax_z_sun,
# plus absorbed-PAR sunlit (parsun_z) — to localize the residual to vcmax(T) vs Aj.
sa = inst.solarabs
println("\n== PHS sunlit vcmax + absorbed PAR ==")
Fv = Dict(2=>53.1953, 3=>91.1934)  # Fortran vcmax_z_sun
for p in 2:3
  jvc = isempty(ps.vcmax_z_phs_patch) ? NaN : Float64(ps.vcmax_z_phs_patch[p,1,1])
  jpar = isempty(sa.parsun_z_patch) ? NaN : Float64(sa.parsun_z_patch[p,1])
  @printf("p%d  vcmax_z_phs_sun J=%.5e F=%.5e ratio=%.4f | parsun_z J=%.5e\n",
          p, jvc, Fv[p], jvc/Fv[p], jpar)
end

# Aj-path inputs vs dump: parsun (absorbed PAR sunlit), jmx25t, SABV
println("\n== Aj-path inputs: Julia vs Fortran dump (before_step_n1757873) ==")
let b=NCDataset("/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/pdump_before_step_n1757873.nc")
  fpar=Float64.(replace(b["parsun"][1,:],missing=>NaN)); fsabv=Float64.(replace(b["SABV_P"][:],missing=>NaN))
  fjmx=haskey(b,"jmx25t") ? Float64.(replace(b["jmx25t"][:],missing=>NaN)) : fill(NaN,3)
  for p in 2:3
    jpar=Float64(sa.parsun_z_patch[p,1]); jsabv=Float64(sa.sabv_patch[p])
    @printf("p%d  parsun J=%.5e F=%.5e ratio=%.4f | SABV J=%.5e F=%.5e ratio=%.4f | jmx25t F=%.4f\n",
            p, jpar, fpar[p], jpar/fpar[p], jsabv, fsabv[p], jsabv/fsabv[p], fjmx[p])
  end
  close(b)
end

# Farquhar co-limitation breakdown (PHS sunlit layer 1) vs Fortran PARITYAJ @1757873
println("\n== Farquhar ac/aj/ap/ag (PHS sunlit) vs Fortran PARITYAJ ==")
FAJ = Dict(  # par_z, je, jmax_z, ac, aj, ap, ag
  2=>(ac=8.68752, aj=15.0944, ap=26.6508, ag=8.28403, je=91.1209, jmax=98.1546, par=227.753),
  3=>(ac=14.2379, aj=34.4410, ap=45.6879, ag=13.7486, je=185.017, jmax=250.146, par=175.292))
for p in 2:3
  f=FAJ[p]
  jac=Float64(ps.ac_phs_patch[p,1,1]); jaj=Float64(ps.aj_phs_patch[p,1,1])
  jap=Float64(ps.ap_phs_patch[p,1,1]); jag=Float64(ps.ag_phs_patch[p,1,1])
  @printf("p%d ac J=%.4f F=%.4f (%.4f) | aj J=%.4f F=%.4f (%.4f) | ap J=%.4f F=%.4f (%.4f) | ag J=%.4f F=%.4f (%.4f)\n",
          p, jac,f.ac,jac/f.ac, jaj,f.aj,jaj/f.aj, jap,f.ap,jap/f.ap, jag,f.ag,jag/f.ag)
end

# ac-denominator inputs: O2(oair), pbot, kc/ko/cp (Julia) vs Fortran dump
println("\n== ac M-M inputs: Julia vs Fortran (dump po2_240=16493.0 pbot=78983.5 cp_F~2.57) ==")
let a=inst.atm2lnd
  oair = Float64(a.forc_po2_grc[1])
  pbot = hasproperty(a,:forc_pbot_downscaled_col) && !isempty(a.forc_pbot_downscaled_col) ? Float64(a.forc_pbot_downscaled_col[1]) : NaN
  @printf("Julia oair(po2)=%.4f  pbot=%.4f\n", oair, pbot)
  for p in 3:3
    kc=Float64(ps.kc_patch[p]); ko=Float64(ps.ko_patch[p]); cp=Float64(ps.cp_patch[p])
    X = kc*(1+oair/ko)
    @printf("p%d grass: kc=%.5f ko=%.5f cp=%.5f | X=kc(1+oair/ko)=%.4f (Fortran-derived 118.89)\n",p,kc,ko,cp,X)
  end
end
