# ==========================================================================
# test_cn_active_layer_roll.jl — the 1 January active-layer roll on the CN path.
#
# Fortran ActiveLayerMod.F90:114 rolls the annual maxima on the first step of
# 1 January (NH; SH 1 July, :126):
#     if (mon == 1 .and. day == 1 .and. sec/dtime == 1)
#        altmax_lastyear_indx(c) = altmax_indx(c); altmax_indx(c) = 0 ...
# SoilBiogeochemVerticalProfileMod.F90:145-214 then integrates the N-deposition,
# N-fixation, leaf/stem and root profiles over 1:altmax_lastyear_indx. With the
# cold-start index 0 every profile is put in the top decomposition layer.
#
# clm_drv! used to call alt_calc! with a constant `sec = 0` (only the FATES path
# passed the real seconds-of-day, PR #357), so a CN run kept the cold-start index
# for ever and every year from year 2 on used top-layer-only profiles.
#
# A. alt_calc! roll → soil_bgc_vertical_profile! spreads the profiles (always runs).
# B. The native restart round-trips the active-layer state (always runs).
# C. A live CN clm_drv! step at 1 January 01:00 rolls the index and the profile
#    the same step uses it (data-gated: Bow site files).
# ==========================================================================
using Test
using CLM
using Dates
using NCDatasets
include(joinpath(@__DIR__, "testdata.jl"))

@testset "CN active-layer roll reaches the vertical profiles" begin
    @testset "A. roll → vertical profile" begin
        saved_zsoi = copy(CLM.zsoi[])
        saved_vp = (CLM.varpar.nlevgrnd, CLM.varpar.nlevsno)
        nlevgrnd = 10; nlevsno = 5; nlev = 5
        try
            CLM.varpar.nlevgrnd = nlevgrnd
            CLM.varpar.nlevsno  = nlevsno
            CLM.zsoi[] = collect(range(0.05, 3.0; length = nlevgrnd))
            al = CLM.ActiveLayerData{Float64}()
            CLM.active_layer_init!(al, 1)
            # cold start: nothing rolled yet; last year thawed down to layer 8
            al.alt_col .= 0.0; al.altmax_col .= 2.0; al.altmax_lastyear_col .= 0.0
            al.alt_indx_col .= 0; al.altmax_indx_col .= 8; al.altmax_lastyear_indx_col .= 0
            temp = CLM.TemperatureData()
            temp.t_soisno_col = fill(CLM.TFRZ + 2.0, 1, nlevsno + nlevgrnd)   # thawed
            col = CLM.ColumnData(); col.gridcell = [1]
            col.nbedrock = [nlev]; col.is_fates = [false]
            grc = CLM.GridcellData(); grc.lat = [deg2rad(51.36)]
            dtime = 3600

            bgc = CLM.SoilBiogeochemStateData()
            bgc.leaf_prof_patch = zeros(1, nlev); bgc.froot_prof_patch = zeros(1, nlev)
            bgc.croot_prof_patch = zeros(1, nlev); bgc.stem_prof_patch = zeros(1, nlev)
            bgc.nfixation_prof_col = zeros(1, nlev); bgc.ndep_prof_col = zeros(1, nlev)
            ss = CLM.SoilStateData(); ss.crootfr_patch = fill(1.0 / nlev, 1, nlev)
            pch = CLM.PatchData(); pch.itype = [1]; pch.column = [1]; pch.wtcol = [1.0]
            dz = fill(0.1, nlev); zs = [(j - 0.5) * 0.1 for j in 1:nlev]
            prof!() = CLM.soil_bgc_vertical_profile!(bgc, al, ss, col, pch;
                mask_bgc_soilc = trues(1), mask_bgc_vegp = trues(1),
                nlevdecomp = nlev, nlevdecomp_full = nlev, dzsoi_decomp = dz, zsoi = zs)

            # The old driver constant: no roll, so the profiles stay in layer 1.
            CLM.alt_calc!(al, trues(1), temp, col, grc; mon = 1, day = 1, sec = 0, dtime = dtime)
            @test al.altmax_lastyear_indx_col[1] == 0
            prof!()
            @test bgc.nfixation_prof_col[1, 1] ≈ 1 / dz[1]
            @test all(iszero, bgc.nfixation_prof_col[1, 2:end])
            @test all(iszero, bgc.froot_prof_patch[1, 2:end])

            # The first step of the year (sec == dtime) rolls last year's maximum.
            CLM.alt_calc!(al, trues(1), temp, col, grc; mon = 1, day = 1, sec = dtime, dtime = dtime)
            @test al.altmax_lastyear_indx_col[1] == 8
            @test al.altmax_lastyear_col[1] == 2.0
            prof!()
            # roots spread over 1:min(8, nlevdecomp), each profile integrates to 1
            @test all(>(0), bgc.nfixation_prof_col[1, :])
            @test all(>(0), bgc.froot_prof_patch[1, :])
            @test sum(bgc.nfixation_prof_col[1, :] .* dz) ≈ 1.0
            @test sum(bgc.froot_prof_patch[1, :] .* dz) ≈ 1.0
            @test bgc.nfixation_prof_col[1, 1] < 1 / dz[1]
        finally
            CLM.zsoi[] = saved_zsoi
            CLM.varpar.nlevgrnd, CLM.varpar.nlevsno = saved_vp
        end
    end

    @testset "B. native restart keeps altmax_lastyear" begin
        names = [rv.name for rv in CLM._restart_registry_biogeophys()]
        for nm in ("altmax", "altmax_lastyear", "altmax_indx", "altmax_lastyear_indx")
            @test nm in names
        end
    end

    @testset "C. live CN driver: 1 January 01:00 rolls and the profile uses it" begin
        bow = joinpath(bow_domain_dir(), "optimization/CLM/dds_run_1/final_evaluation/settings/CLM/parameters")
        fsurdat = joinpath(bow, "surfdata_clm.nc")
        fparam  = joinpath(bow, "clm5_params.nc")
        forcing = joinpath(bow_domain_dir(), "data/forcing/CLM_input/clmforc.2003.nc")
        if !(isfile(fsurdat) && isfile(fparam) && isfile(forcing))
            testdata_missing("CN active-layer roll (live driver)", fsurdat, fparam, forcing)
            return
        end
        CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
        CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
        # Clock at 2002-12-31 23:00: the first step ends at 2003-01-01 00:00
        # (sec = 0, no roll), the second at 01:00 (sec == dtime, roll).
        (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat = fsurdat, paramfile = fparam,
            start_date = DateTime(2002, 12, 31, 23), dtime = 3600, use_cn = true,
            use_luna = true, use_bedrock = true, use_aquifer_layer = false, h2osfcflag = 0)
        config  = CLM.CLMDriverConfig(use_cn = true, use_aquifer_layer = false, use_luna = true)
        filt_ia = CLM.clump_filter_inactive_and_active
        al = inst.active_layer
        c = findfirst(filt.bgc_soilc)
        @test c !== nothing
        @test al.altmax_lastyear_indx_col[c] == 0          # cold start
        bgcs = inst.soilbiogeochem_state
        dzd = CLM.dzsoi_decomp[]
        nlev = CLM.varpar.nlevdecomp

        fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, forcing)
        ng, nc = bounds.endg, bounds.endc
        dtime = 3600.0
        for step in 1:2
            cur = tm.current_date
            calday = CLM.get_curr_calday(tm)
            (declin, _)   = CLM.compute_orbital(calday)
            (declinm1, _) = CLM.compute_orbital(calday - dtime / CLM.SECSPDAY)
            step == 1 && CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
            CLM.advance_timestep!(tm)
            CLM.read_forcing_step!(fr, inst.atm2lnd, DateTime(2003, 1, 1) + (cur - DateTime(2002, 12, 31, 23)), ng, nc)
            CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
            (yr, mon, dy, tod) = CLM.get_curr_date(tm)
            if step == 1
                @test (mon, dy, tod) == (1, 1, 0)
                # last year's maximum, as a thawed Bow column reaches by summer
                al.altmax_col[c] = 3.0
                al.altmax_indx_col[c] = 12
            else
                @test (mon, dy, tod) == (1, 1, 3600)
            end
            CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, calday + dtime / CLM.SECSPDAY,
                declin, declin, CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
                nstep = tm.nstep, is_first_step = false,
                is_beg_curr_day = CLM.is_beg_curr_day(tm), is_end_curr_day = CLM.is_end_curr_day(tm),
                is_beg_curr_year = CLM.is_beg_curr_year(tm), dtime = dtime, mon = mon, day = dy,
                secs = tod, jday = 1, year = yr, photosyns = inst.photosyns)
            if step == 1
                # 00:00 is not the first step of the year (sec/dtime == 0): no roll,
                # and the profiles still sit in the top layer
                @test al.altmax_lastyear_indx_col[c] == 0
                @test bgcs.nfixation_prof_col[c, 1] ≈ 1 / dzd[1]
                @test all(iszero, bgcs.nfixation_prof_col[c, 2:nlev])
            end
        end
        CLM.forcing_reader_close!(fr)

        @test al.altmax_lastyear_indx_col[c] >= 12
        @test al.altmax_lastyear_col[c] >= 3.0
        jmax = min(al.altmax_lastyear_indx_col[c], nlev)
        @test all(>(0), bgcs.nfixation_prof_col[c, 1:jmax])
        @test sum(bgcs.nfixation_prof_col[c, 1:nlev] .* dzd[1:nlev]) ≈ 1.0
        @test bgcs.nfixation_prof_col[c, 1] < 1 / dzd[1]

        # and a native restart written now carries the rolled index forward
        f = tempname() * ".nc"
        try
            CLM.write_restart(inst, f; bounds = bounds, use_cn = true)
            saved = al.altmax_lastyear_indx_col[c]
            al.altmax_lastyear_indx_col .= 0
            CLM.read_restart!(inst, f; bounds = bounds, use_cn = true)
            @test al.altmax_lastyear_indx_col[c] == saved
        finally
            isfile(f) && rm(f)
        end
    end
end
