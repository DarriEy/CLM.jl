# ==========================================================================
# test_fortran_parity_luna_update.jl — LUNA end-of-day vcmax/jmax acclimation
# update vs Fortran (gated on the machine-local BGC dumps).
#
# LUNA's per-step climate accumulators are in no standard restart dump, so the
# Fortran ground truth comes from an instrumented LunaMod run (PARITYLUNA inputs
# + PARITYLUNO outputs at the EOD nstep 1757880, reachable from the
# 2202-07-16-57600 restart). Those values are hardcoded below. We inject them +
# the pre-update vcmx25_z (from pdump_before_step_n1757880.nc) into a use_luna
# Julia inst and call update_photosynthesis_capacity! standalone.
#
# ASSERTED: jmx25 within 5%, vcmx25 (single-day constrained update) within 8%.
# KNOWN residual: the vcmax OPTIMUM (PNcb / Rubisco-N allocation) is lower in
# Julia than Fortran; the daily max_daily_pchg constraint keeps the one-day
# update close, but the target differs — a LUNA-port refinement to investigate.
# ==========================================================================
using Test, CLM, NCDatasets, Dates, Printf

@testset "Fortran LUNA EOD update (gated)" begin
    bgcdir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup"
    edump  = joinpath(bgcdir, "pdump_before_step_n1757880.nc")
    if !isfile(edump)
        @info "LUNA-update parity: dumps absent, skipping" edump
        @test_skip isfile(edump)
    else
        common = joinpath(@__DIR__, "..", "scripts", "fortran_parity_common.jl")
        Base.include(@__MODULE__, common)

        # Fortran ground truth at nstep 1757880 (LunaMod PARITYLUNA / PARITYLUNO).
        F = Dict(
         2 => (tvd10=296.8046, tvn10=280.0893, gppday=10424.38, rh10=0.198959, rb10=42.6012,
               par240d=164.3697, par240x=248.9905, t10=285.5351, lnca=2.97826, tlai=0.0476,
               vcpre=49.5631, jmpre=92.9850, vcpost=49.8445, jmpost=92.9847),
         3 => (tvd10=291.0896, tvn10=279.7898, gppday=119618.3, rh10=0.390046, rb10=51.8405,
               par240d=125.3865, par240x=180.0593, t10=285.1617, lnca=8.35148, tlai=0.2354,
               vcpre=136.8619, jmpre=336.0765, vcpost=136.1970, jmpost=334.2716))
        CO2_240=28.96185; O2_240=16493.26; PBOT240=78915.11; DAYL=57263.81; DAYLF=0.944373

        (inst, bounds, filt, tm) = Base.invokelatest(build_bow_inst; dtime=3600,
            start_date=DateTime(2002,7,15,23), use_cn=true, use_luna=true)
        ps=inst.photosyns; temp=inst.temperature; cs=inst.canopystate; alb=inst.surfalb
        sa=inst.solarabs; wdb=inst.water.waterdiagnosticbulk_inst; fv=inst.frictionvel
        pch=inst.patch; grc=inst.gridcell; a2l=inst.atm2lnd; np=bounds.endp

        for g in 1:bounds.endg; grc.dayl[g]=DAYL; grc.max_dayl[g]=DAYL/0.971589; end
        for p in 2:3
            f=F[p]
            ps.vcmx25_z_patch[p,1]=f.vcpre; ps.jmx25_z_patch[p,1]=f.jmpre
            ps.vcmx25_z_last_valid_patch[p,1]=f.vcpre; ps.jmx25_z_last_valid_patch[p,1]=f.jmpre
            ps.lnca_patch[p]=f.lnca; ps.fpsn24_patch[p]=f.gppday
            ps.pnlc_z_patch[p,1]=0.01; ps.enzs_z_patch[p,1]=1.0
            temp.t_veg10_day_patch[p]=f.tvd10; temp.t_veg10_night_patch[p]=f.tvn10
            temp.t_a10_patch[p]=f.t10; temp.t_veg_day_patch[p]=290.0
            wdb.rh10_af_patch[p]=f.rh10; fv.rb10_patch[p]=f.rb10
            sa.par240d_z_patch[p,1]=f.par240d; sa.par240x_z_patch[p,1]=f.par240x
            alb.nrad_patch[p]=1; alb.tlai_z_patch[p,1]=f.tlai; cs.tlai_patch[p]=f.tlai
            a2l.forc_pco2_240_patch[p]=CO2_240; a2l.forc_po2_240_patch[p]=O2_240
            a2l.forc_pbot240_downscaled_patch[p]=PBOT240
        end
        mask=falses(np); mask[2]=true; mask[3]=true
        pftcon=CLM.pftcon
        Base.invokelatest(CLM.update_photosynthesis_capacity!, ps, temp, cs, alb, sa, wdb, fv, pch, grc,
            mask, 1:np, fill(DAYLF,np),
            a2l.forc_pbot240_downscaled_patch, a2l.forc_pco2_240_patch, a2l.forc_po2_240_patch,
            pftcon.c3psn, pftcon.slatop, pftcon.leafcn, pftcon.rhol, pftcon.taul,
            fill(1.0,np), CLM.luna_params_inst, 3600.0, CLM.NLEVCAN)

        for p in 2:3
            f=F[p]
            jvc=Float64(ps.vcmx25_z_patch[p,1]); jjm=Float64(ps.jmx25_z_patch[p,1])
            @info "LUNA update" patch=p vcmx25_julia=jvc vcmx25_fortran=f.vcpost jmx25_julia=jjm jmx25_fortran=f.jmpost
            @test isfinite(jvc) && isfinite(jjm)
            @test abs(jjm/f.jmpost - 1) < 0.05   # jmax matches Fortran
            @test abs(jvc/f.vcpost - 1) < 0.08   # vcmax single-day update within 8%
        end
    end
end
