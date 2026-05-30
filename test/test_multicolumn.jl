# ==========================================================================
# test_multicolumn.jl — batched (N-gridcell) run reproduces single-gridcell run
#
# clm_initialize!(; ncopies=N) tiles the single-point surface data into N
# identical independent gridcells (the data-parallel layout the GPU port batches
# over). N tiled gridcells must reproduce, per gridcell, a single-gridcell run.
# ==========================================================================

using Test
using CLM

@testset "Multi-column harness (ncopies)" begin
    fsurdat = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
    paramfile = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"
    if !isfile(fsurdat) || !isfile(paramfile)
        @warn "Skipping multi-column test: input files not found"
        @test true
        return
    end

    function setup_forcing!(a2l, T, ng)
        for g in 1:ng
            a2l.forc_t_not_downscaled_grc[g] = T
            a2l.forc_pbot_not_downscaled_grc[g] = 85000.0
            a2l.forc_th_not_downscaled_grc[g] = T * (100000.0 / 85000.0)^(CLM.RAIR / CLM.CPAIR)
            a2l.forc_rho_not_downscaled_grc[g] = 85000.0 / (CLM.RAIR * T)
            a2l.forc_lwrad_not_downscaled_grc[g] = 300.0
            a2l.forc_vp_grc[g] = 800.0
            a2l.forc_hgt_grc[g] = 30.0
            a2l.forc_hgt_u_grc[g] = 30.0; a2l.forc_hgt_t_grc[g] = 30.0; a2l.forc_hgt_q_grc[g] = 30.0
            a2l.forc_topo_grc[g] = 0.0
            a2l.forc_wind_grc[g] = 3.0; a2l.forc_u_grc[g] = 3.0; a2l.forc_v_grc[g] = 0.0
            for b in 1:CLM.NUMRAD
                a2l.forc_solad_not_downscaled_grc[g, b] = 200.0
                a2l.forc_solai_grc[g, b] = 80.0
            end
            a2l.forc_solar_not_downscaled_grc[g] = 560.0
            a2l.forc_rain_not_downscaled_grc[g] = 0.0001
            a2l.forc_snow_not_downscaled_grc[g] = 0.0
        end
    end

    function run(ncopies, ksteps)
        (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=fsurdat, paramfile=paramfile, ncopies=ncopies)
        ng = bounds.endg
        config = CLM.CLMDriverConfig()
        fia = CLM.clump_filter_inactive_and_active
        (declin, _) = CLM.compute_orbital(120.0)
        nsw = 120.0 + 1800.0 / CLM.SECSPDAY
        setup_forcing!(inst.atm2lnd, 285.0, ng)
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        for n in 1:ksteps
            CLM.clm_drv!(config, inst, filt, fia, bounds, true, nsw, declin, declin,
                CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
                nstep=n, is_first_step=(n == 1), is_beg_curr_day=(n == 1),
                dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
        end
        return inst, bounds
    end

    # Group elements by gridcell (columns/patches are ordered landunit-major).
    function matches(arrN, ref, owner_grc, ng)
        for g in 1:ng
            idx = findall(==(g), owner_grc)
            length(idx) == length(ref) || return false
            for i in eachindex(ref)
                isequal(arrN[idx[i]], ref[i]) || return false
            end
        end
        return true
    end

    inst1, b1 = run(1, 5)
    ref_tg = copy(inst1.temperature.t_grnd_col)
    ref_lh = copy(inst1.energyflux.eflx_lh_tot_patch)

    N = 3
    instN, bN = run(N, 5)

    @test bN.endg == N
    @test bN.endc == N * b1.endc
    @test bN.endp == N * b1.endp
    @test matches(instN.temperature.t_grnd_col, ref_tg, instN.column.gridcell[1:bN.endc], N)
    @test matches(instN.energyflux.eflx_lh_tot_patch, ref_lh, instN.patch.gridcell[1:bN.endp], N)
    println("  multi-column: $N tiled gridcells reproduce single-gridcell run exactly")
end
