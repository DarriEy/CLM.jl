# Validate LUNA's end-of-day vcmax/jmax acclimation update against Fortran.
#
# LUNA's per-step climate accumulators are in no standard dump, so we instrument
# Fortran LunaMod (PARITYLUNA inputs + PARITYLUNO outputs at the EOD n1757880,
# reachable from the 2202-07-16-57600 restart). Here we inject those exact inputs
# + the pre-update vcmx25_z (before_step_n1757880) into a use_luna Julia inst and
# call update_photosynthesis_capacity! standalone, isolating the update function.
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

# Fortran ground truth at nstep 1757880 (from LunaMod PARITYLUNA/PARITYLUNO).
# pre-update vcmx25_z/jmx25_z + the persisted pnlc_z/enzs_z from
# pdump_before_step_n1757880.nc. pnlc IS the LUNA optimizer's initial guess
# (PNlcold): the nitrogen-allocation loop is a start-dependent discrete hill-climb
# (increase_flag latches on the first iteration), so injecting Fortran's REAL
# per-patch pnlc_z — not a uniform placeholder — is required for a valid vcmx25_z
# comparison. With it, Julia matches Fortran to <1e-4 (see PR: the earlier ~5%
# vcmax "divergence" was this harness placeholder, NOT a LunaMod port error).
const F = Dict(
 2 => (tvd10=296.8046, tvn10=280.0893, gppday=10424.38, rh10=0.198959, rb10=42.6012,
       par240d=164.3697, par240x=248.9905, t10=285.5351, lnca=2.97826, tlai=0.0476,
       vcpre=49.5631, jmpre=92.9850, vcpost=49.8445, jmpost=92.9847,
       pnlc=0.0499999999999986, enzs=1.00900022224095),
 3 => (tvd10=291.0896, tvn10=279.7898, gppday=119618.3, rh10=0.390046, rb10=51.8405,
       par240d=125.3865, par240x=180.0593, t10=285.1617, lnca=8.35148, tlai=0.2354,
       vcpre=136.8619, jmpre=336.0765, vcpost=136.1970, jmpost=334.2716,
       pnlc=0.0899999999999997, enzs=1.00494913755897))
const CO2_240=28.96185; const O2_240=16493.26; const PBOT240=78915.11; const DAYLF=0.944373

(inst, bounds, filt, tm) = build_bow_inst(; dtime=3600, start_date=DateTime(2002,7,15,23),
                                            use_cn=true, use_luna=true)
ps=inst.photosyns; temp=inst.temperature; cs=inst.canopystate; alb=inst.surfalb
sa=inst.solarabs; wdb=inst.water.waterdiagnosticbulk_inst; fv=inst.frictionvel
pch=inst.patch; grc=inst.gridcell; a2l=inst.atm2lnd
np=bounds.endp

for p in 2:3
    f=F[p]
    ps.vcmx25_z_patch[p,1]=f.vcpre;  ps.jmx25_z_patch[p,1]=f.jmpre
    ps.vcmx25_z_last_valid_patch[p,1]=f.vcpre; ps.jmx25_z_last_valid_patch[p,1]=f.jmpre
    ps.lnca_patch[p]=f.lnca; ps.fpsn24_patch[p]=f.gppday
    ps.pnlc_z_patch[p,1]=f.pnlc; ps.enzs_z_patch[p,1]=f.enzs
    temp.t_veg10_day_patch[p]=f.tvd10; temp.t_veg10_night_patch[p]=f.tvn10
    temp.t_a10_patch[p]=f.t10; temp.t_veg_day_patch[p]=290.0  # non-SPVAL → pass first-day guard
    wdb.rh10_af_patch[p]=f.rh10; fv.rb10_patch[p]=f.rb10
    sa.par240d_z_patch[p,1]=f.par240d; sa.par240x_z_patch[p,1]=f.par240x
    alb.nrad_patch[p]=1; alb.tlai_z_patch[p,1]=f.tlai; cs.tlai_patch[p]=f.tlai
    a2l.forc_pco2_240_patch[p]=CO2_240; a2l.forc_po2_240_patch[p]=O2_240
    a2l.forc_pbot240_downscaled_patch[p]=PBOT240
end

# gridcell daylength (used by the update for hourpd + the 10-day leaf-temp mean).
for g in 1:bounds.endg; grc.dayl[g] = 57263.81; grc.max_dayl[g] = 57263.81 / 0.971589; end
mask=falses(np); mask[2]=true; mask[3]=true
dayl_factor=fill(DAYLF, np); o3coefjmax=fill(1.0, np)
pftcon=CLM.pftcon

CLM.update_photosynthesis_capacity!(ps, temp, cs, alb, sa, wdb, fv, pch, grc,
    mask, 1:np, dayl_factor,
    a2l.forc_pbot240_downscaled_patch, a2l.forc_pco2_240_patch, a2l.forc_po2_240_patch,
    pftcon.c3psn, pftcon.slatop, pftcon.leafcn, pftcon.rhol, pftcon.taul,
    o3coefjmax, CLM.luna_params_inst, 3600.0, CLM.NLEVCAN)

println("== LUNA EOD update (nstep 1757880): Julia vs Fortran ==")
@printf("%-6s %-10s %12s %12s %12s %10s\n","patch","field","pre","Julia post","Fortran post","ratio")
for p in 2:3
    f=F[p]
    jvc=Float64(ps.vcmx25_z_patch[p,1]); jjm=Float64(ps.jmx25_z_patch[p,1])
    @printf("p%-5d %-10s %12.5f %12.5f %12.5f %10.5f\n",p,"vcmx25_z",f.vcpre,jvc,f.vcpost,jvc/f.vcpost)
    @printf("%-6s %-10s %12.5f %12.5f %12.5f %10.5f\n","","jmx25_z",f.jmpre,jjm,f.jmpost,jjm/f.jmpost)
end
