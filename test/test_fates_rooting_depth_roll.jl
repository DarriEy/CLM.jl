# ==========================================================================
# test_fates_rooting_depth_roll.jl — last year's active layer must reach FATES.
#
# Fortran: ActiveLayerMod.F90:114 rolls the annual maxima on the first step of
# 1 January (NH; SH 1 July, :126),
#     if (mon == 1 .and. day == 1 .and. sec/dtime == 1)
#        altmax_lastyear_indx(c) = altmax_indx(c); altmax_indx(c) = 0 ...
# and clmfates_interfaceMod.F90:1131-1132 hands that index to FATES every day:
#     bc_in(s)%max_rooting_depth_index_col = min(nlevsoil, altmax_lastyear_indx_col(c))
#
# The port called alt_calc! with `sec = 0`, so `sec/dtime == 1` could never be true.
# `altmax_lastyear_indx_col` kept its cold-start 0 (ActiveLayerMod.F90:305) forever,
# and FATES ran every year as if it were year 1: the drought-phenology soil-moisture
# memory, btran and the root-litter profiles were all confined to the top 2 soil
# layers. On the 4-year fates_parity_1pt reference this was "Defect 2": year-1
# recruitment matched Fortran, then years 2-4 ran 1.4-2.2x high, because the frozen
# top layers switched the hard drought-deciduous grasses (PFTs 13, 14) off ~2 months
# early in autumn and back on ~2 months early in spring
# (docs/FATES_MULTIYEAR_ENERGY.md).
#
# A. alt_calc! roll semantics (always runs).
# B. The live FATES driver passes the real seconds-of-day: a 1 January 01:00 step
#    rolls the index, and the daily FATES bc pack then sees it (data-gated; the
#    Bow site files live outside the repo).
# ==========================================================================
using Test
using CLM
using Dates
include(joinpath(@__DIR__, "testdata.jl"))

@testset "FATES rooting depth: annual active-layer roll" begin
    @testset "A. alt_calc! rolls only on the first step of 1 January" begin
        saved_zsoi = copy(CLM.zsoi[])
        saved_vp = (CLM.varpar.nlevgrnd, CLM.varpar.nlevsno)
        nlevgrnd = 25
        nlevsno  = 12
        try
            CLM.varpar.nlevgrnd = nlevgrnd
            CLM.varpar.nlevsno  = nlevsno
            CLM.zsoi[] = collect(range(0.01, 40.0; length = nlevgrnd))
            al = CLM.ActiveLayerData{Float64}()
            CLM.active_layer_init!(al, 1)
            al.alt_col .= 0.0; al.altmax_col .= 3.0; al.altmax_lastyear_col .= 0.0
            al.alt_indx_col .= 0; al.altmax_indx_col .= 12; al.altmax_lastyear_indx_col .= 0
            temp = CLM.TemperatureData()
            temp.t_soisno_col = fill(CLM.TFRZ - 5.0, 1, nlevsno + nlevgrnd)   # all frozen
            col = CLM.ColumnData(); col.gridcell = [1]
            grc = CLM.GridcellData(); grc.lat = [deg2rad(51.36)]
            mask = trues(1)
            dtime = 3600

            # the constant the driver used to pass: never rolls
            CLM.alt_calc!(al, mask, temp, col, grc; mon = 1, day = 1, sec = 0, dtime = dtime)
            @test al.altmax_lastyear_indx_col[1] == 0
            @test al.altmax_indx_col[1] == 12

            # the first step of the year rolls last year's maximum
            CLM.alt_calc!(al, mask, temp, col, grc; mon = 1, day = 1, sec = dtime, dtime = dtime)
            @test al.altmax_lastyear_indx_col[1] == 12
            @test al.altmax_lastyear_col[1] == 3.0
            @test al.altmax_indx_col[1] == 0          # frozen column: new-year max restarts at 0

            # any other step of 1 January leaves it alone
            al.altmax_indx_col .= 7
            CLM.alt_calc!(al, mask, temp, col, grc; mon = 1, day = 1, sec = 2dtime, dtime = dtime)
            @test al.altmax_lastyear_indx_col[1] == 12

            # NH columns do not roll on 1 July
            CLM.alt_calc!(al, mask, temp, col, grc; mon = 7, day = 1, sec = dtime, dtime = dtime)
            @test al.altmax_lastyear_indx_col[1] == 12
        finally
            CLM.zsoi[] = saved_zsoi
            CLM.varpar.nlevgrnd, CLM.varpar.nlevsno = saved_vp
        end
    end

    @testset "B. live FATES driver: 1 January roll reaches max_rooting_depth_index_col" begin
        sitedir = symfluence_path("domain_Bow_at_Banff_lumped_era5", "settings", "CLM", "parameters")
        fsurdat = joinpath(sitedir, "surfdata_clm.nc")
        fparam  = joinpath(sitedir, "clm5_params.nc")
        forcing = symfluence_path("domain_Bow_at_Banff_lumped_era5", "data", "forcing",
                                  "CLM_input", "clmforc.2003.nc")
        if !(isfile(fsurdat) && isfile(fparam) && isfile(forcing))
            testdata_missing("FATES rooting-depth roll (live driver)", fsurdat, fparam, forcing)
            return
        end
        ENV["SYMFLUENCE_DATA"] = symfluence_data_root()
        m = Module(:_FatesRootRollHarness)
        Core.eval(m, :(include(p::AbstractString) = Base.include($m, p)))
        Base.include(m, joinpath(@__DIR__, "..", "scripts", "fates_fortran_parity.jl"))
        inst, bounds, filt, filt_ia, config, _tm = Base.invokelatest(getglobal(m, :build))
        @test config.use_fates
        al = inst.active_layer
        nlevsoil = inst.fates.bc_in[1].nlevsoil
        # cold start: FATES roots live in the top layers until the first roll
        @test al.altmax_lastyear_indx_col[1] == 0
        # pretend last year's active layer thawed down to layer 15 (the Bow column
        # thaws through in summer; any index in 3..nlevsoil exercises the roll)
        al.altmax_col[1] = 2.0
        al.altmax_indx_col[1] = 15

        fr = CLM.ForcingReader()
        CLM.forcing_reader_init!(fr, forcing)
        dtime = 3600.0
        ng, nc = bounds.endg, bounds.endc
        # one step whose seconds-of-day satisfy `sec/dtime == 1` on 1 January
        cur = DateTime(2003, 1, 1, 1)
        tod = 3600
        calday = 1 + tod / CLM.SECSPDAY
        (declin, _) = CLM.compute_orbital(calday)
        (declinm1, _) = CLM.compute_orbital(calday - dtime / CLM.SECSPDAY)
        CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
        CLM.read_forcing_step!(fr, inst.atm2lnd, cur, ng, nc)
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, calday + dtime / CLM.SECSPDAY,
            declin, declin, CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep = 2, is_first_step = false, is_beg_curr_day = false, is_end_curr_day = false,
            is_beg_curr_year = false, dtime = dtime, mon = 1, day = 1, secs = tod,
            jday = 1, year = 2003, photosyns = inst.photosyns)
        CLM.forcing_reader_close!(fr)

        @test al.altmax_lastyear_indx_col[1] == 15
        nlevdecomp = length(inst.fates.bc_in[1].w_scalar_sisl)
        CLM.fates_pack_bcin_daily!(inst; s = 1, c = 1, nlevsoil = nlevsoil, nlevdecomp = nlevdecomp)
        @test inst.fates.bc_in[1].max_rooting_depth_index_col == min(nlevsoil, 15)
        @test inst.fates.bc_in[1].max_rooting_depth_index_col > 2
    end
end
