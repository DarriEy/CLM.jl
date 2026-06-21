@testset "Crop GDD accumulation" begin

    TFRZ = CLM.TFRZ            # 273.15 K
    DTIME = 86400              # one full day per step → dtime/SECSPDAY == 1

    # Build a manager with a single runaccum GDD field for `basetemp`.
    function make_gdd_mgr(basetemp::Int, npts::Int; active=fill(true, npts))
        mgr = CLM.AccumManager()
        CLM.init_accum_field!(mgr;
            name = "GDD$(basetemp)",
            units = "K",
            desc  = "growing degree days base $(basetemp)C",
            accum_type = "runaccum",
            accum_period = typemax(Int),   # runaccum: never auto-reset by period
            numlev = 1,
            subgrid_type = "pft",
            init_value = 0.0,
            active = active,
            npts = npts)
        return mgr
    end

    # A single NH grid cell (lat 45N). One patch (begp=endp=1).
    latdeg   = [45.0]
    gridcell = [1]
    active   = [true]
    itype    = [0]

    # ------------------------------------------------------------------
    @testset "basic accumulation (base 0, NH summer)" begin
        mgr  = make_gdd_mgr(0, 1)
        temp = CLM.TemperatureData()
        temp.t_ref2m_patch = [TFRZ + 10.0]   # 10 C above freezing
        gddx = [0.0]

        # July (in NH accumulation season). Each daily step adds 10 GDD
        # (max(0, min(26, 10)) * 86400/86400).
        for nstep in 1:3
            CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr, gddx, 0;
                begp=1, endp=1, month=7, day=15+nstep, secs=DTIME,
                dtime=DTIME, nstep=nstep, jday=196+nstep,
                latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        end
        @test gddx[1] ≈ 30.0 atol=1e-12   # 3 days * 10 GDD
    end

    # ------------------------------------------------------------------
    @testset "freezing cutoff (no negative GDD)" begin
        mgr  = make_gdd_mgr(0, 1)
        temp = CLM.TemperatureData()
        temp.t_ref2m_patch = [TFRZ - 5.0]    # 5 C below freezing → GDD clamped to 0
        gddx = [0.0]

        for nstep in 1:3
            CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr, gddx, 0;
                begp=1, endp=1, month=7, day=15+nstep, secs=DTIME,
                dtime=DTIME, nstep=nstep, jday=196+nstep,
                latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        end
        @test gddx[1] == 0.0
    end

    # ------------------------------------------------------------------
    @testset "max_accum cap (26 for base 0, 30 for base 8)" begin
        # Very hot day: 50 C above freezing.
        temp = CLM.TemperatureData()
        temp.t_ref2m_patch = [TFRZ + 50.0]

        mgr0 = make_gdd_mgr(0, 1)
        g0 = [0.0]
        CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr0, g0, 0;
            begp=1, endp=1, month=7, day=16, secs=DTIME,
            dtime=DTIME, nstep=1, jday=197,
            latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        @test g0[1] ≈ 26.0 atol=1e-12      # capped at max_accum=26

        # base 8: subtract 8, then cap at 30. 50-8 = 42 → capped to 30.
        mgr8 = make_gdd_mgr(8, 1)
        g8 = [0.0]
        CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr8, g8, 8;
            begp=1, endp=1, month=7, day=16, secs=DTIME,
            dtime=DTIME, nstep=1, jday=197,
            latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        @test g8[1] ≈ 30.0 atol=1e-12      # capped at max_accum=30
    end

    # ------------------------------------------------------------------
    @testset "base 8 subtracts base temperature" begin
        mgr  = make_gdd_mgr(8, 1)
        temp = CLM.TemperatureData()
        temp.t_ref2m_patch = [TFRZ + 20.0]   # 20 C; base 8 → 12 GDD/day
        gddx = [0.0]
        CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr, gddx, 8;
            begp=1, endp=1, month=7, day=16, secs=DTIME,
            dtime=DTIME, nstep=1, jday=197,
            latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        @test gddx[1] ≈ 12.0 atol=1e-12
    end

    # ------------------------------------------------------------------
    @testset "Jan-1 reset" begin
        mgr  = make_gdd_mgr(0, 1)
        temp = CLM.TemperatureData()
        temp.t_ref2m_patch = [TFRZ + 10.0]
        gddx = [0.0]

        # Accumulate some GDD in July.
        for nstep in 1:3
            CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr, gddx, 0;
                begp=1, endp=1, month=7, day=15+nstep, secs=DTIME,
                dtime=DTIME, nstep=nstep, jday=196+nstep,
                latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        end
        @test gddx[1] ≈ 30.0 atol=1e-12

        # Jan 1, first timestep of day (secs == dtime): accumulator resets.
        # runaccum cannot reset AND accumulate in the same call, so this
        # step zeroes the field (increment applied next step).
        CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr, gddx, 0;
            begp=1, endp=1, month=1, day=1, secs=DTIME,
            dtime=DTIME, nstep=4, jday=1,
            latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        @test gddx[1] == 0.0

        # Jan is NH out-of-season → next step adds nothing and stays 0.
        CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr, gddx, 0;
            begp=1, endp=1, month=1, day=2, secs=DTIME,
            dtime=DTIME, nstep=5, jday=2,
            latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        @test gddx[1] == 0.0
    end

    # ------------------------------------------------------------------
    @testset "out-of-season contributes zero (NH winter)" begin
        mgr  = make_gdd_mgr(0, 1)
        temp = CLM.TemperatureData()
        temp.t_ref2m_patch = [TFRZ + 10.0]   # warm, but it's December in the NH
        gddx = [0.0]
        # December (month 12) in NH is out of season (season is Apr–Sep).
        CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr, gddx, 0;
            begp=1, endp=1, month=12, day=15, secs=DTIME,
            dtime=DTIME, nstep=1, jday=349,
            latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        @test gddx[1] == 0.0
    end

    # ------------------------------------------------------------------
    @testset "Southern hemisphere season is flipped" begin
        mgr  = make_gdd_mgr(0, 1)
        temp = CLM.TemperatureData()
        temp.t_ref2m_patch = [TFRZ + 10.0]
        gddx = [0.0]
        sh_lat = [-30.0]
        # January in the SH IS in season (Oct–Mar) → accumulates.
        CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr, gddx, 0;
            begp=1, endp=1, month=1, day=2, secs=DTIME,
            dtime=DTIME, nstep=1, jday=2,
            latdeg=sh_lat, gridcell=gridcell, active=active, itype=itype)
        @test gddx[1] ≈ 10.0 atol=1e-12
    end

    # ------------------------------------------------------------------
    @testset "read-in gdd20 season overrides latitude (crop patch)" begin
        npcropmin = 15
        # Crop patch (itype >= npcropmin) with a valid read-in season DOY 100–200.
        mgr  = make_gdd_mgr(0, 1)
        temp = CLM.TemperatureData()
        temp.t_ref2m_patch = [TFRZ + 10.0]
        crop_itype = [npcropmin]
        season_start = [100.0]
        season_end   = [200.0]

        # NH December (month 12) → latitude fallback would say OUT of season,
        # but jday=150 is inside the read-in season 100–200 → accumulates.
        gddx = [0.0]
        CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr, gddx, 0;
            begp=1, endp=1, month=12, day=15, secs=DTIME,
            dtime=DTIME, nstep=1, jday=150,
            latdeg=latdeg, gridcell=gridcell, active=active, itype=crop_itype,
            npcropmin=npcropmin,
            gdd20_season_start=season_start, gdd20_season_end=season_end)
        @test gddx[1] ≈ 10.0 atol=1e-12

        # jday=250 is OUTSIDE the read-in season → no accumulation even in
        # NH summer (month 7) where the latitude fallback would say in-season.
        mgr2  = make_gdd_mgr(0, 1)
        gddx2 = [0.0]
        CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr2, gddx2, 0;
            begp=1, endp=1, month=7, day=15, secs=DTIME,
            dtime=DTIME, nstep=1, jday=250,
            latdeg=latdeg, gridcell=gridcell, active=active, itype=crop_itype,
            npcropmin=npcropmin,
            gdd20_season_start=season_start, gdd20_season_end=season_end)
        @test gddx2[1] == 0.0
    end

    # ------------------------------------------------------------------
    @testset "inactive patch is skipped" begin
        mgr  = make_gdd_mgr(0, 1; active=[false])
        temp = CLM.TemperatureData()
        temp.t_ref2m_patch = [TFRZ + 10.0]
        gddx = [0.0]
        CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr, gddx, 0;
            begp=1, endp=1, month=7, day=16, secs=DTIME,
            dtime=DTIME, nstep=1, jday=197,
            latdeg=latdeg, gridcell=gridcell, active=[false], itype=itype)
        @test gddx[1] == 0.0   # inactive: update_accum_field ignores it
    end

    # ------------------------------------------------------------------
    @testset "partial timestep scales by dtime/SECSPDAY" begin
        mgr  = make_gdd_mgr(0, 1)
        temp = CLM.TemperatureData()
        temp.t_ref2m_patch = [TFRZ + 10.0]
        gddx = [0.0]
        half_day = 43200   # 12 h
        # secs != dtime so not the Jan-1 reset trigger even if month/day matched.
        CLM.temperature_update_acc_vars_crop_gdds!(temp, mgr, gddx, 0;
            begp=1, endp=1, month=7, day=16, secs=half_day,
            dtime=half_day, nstep=1, jday=197,
            latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        @test gddx[1] ≈ 10.0 * (half_day / CLM.SECSPDAY) atol=1e-12  # = 5.0
    end

    # ==================================================================
    # Wired path: temperature_init_acc_buffer! + temperature_update_acc_vars!
    # ==================================================================

    # Build a TemperatureData sized for `np` patches with the GDD fields
    # initialized (the driver-call path reads/writes temp.gddN_patch in place).
    function make_temp(np::Int; t2m=TFRZ + 10.0)
        temp = CLM.TemperatureData{Float64}()
        # arrays touched by the wired update (kernel + crop-GDD)
        temp.t_veg_patch    = fill(t2m, np)
        temp.t_ref2m_patch  = fill(t2m, np)
        temp.t_veg24_patch  = fill(NaN, np)
        temp.t_veg240_patch = fill(NaN, np)
        temp.t_a10_patch    = fill(NaN, np)
        temp.gdd0_patch     = zeros(np)
        temp.gdd8_patch     = zeros(np)
        temp.gdd10_patch    = zeros(np)
        temp.gdd020_patch   = zeros(np)
        temp.gdd820_patch   = zeros(np)
        temp.gdd1020_patch  = zeros(np)
        return temp
    end

    # ------------------------------------------------------------------
    @testset "init_acc_buffer registers GDD fields only when use_crop" begin
        mgr = CLM.AccumManager()
        CLM.temperature_init_acc_buffer!(mgr, 1; use_crop=false)
        @test isempty(mgr.fields)

        mgr2 = CLM.AccumManager()
        CLM.temperature_init_acc_buffer!(mgr2, 1; use_crop=true)
        names = [f.name for f in mgr2.fields]
        @test names == ["GDD0", "GDD8", "GDD10", "GDD020", "GDD820", "GDD1020"]
        # runaccum vs runmean types and 20-yr period
        @test mgr2.fields[CLM.find_field(mgr2, "GDD0")].acctype == CLM.ACCTYPE_RUNACCUM
        gdd020 = mgr2.fields[CLM.find_field(mgr2, "GDD020")]
        @test gdd020.acctype == CLM.ACCTYPE_RUNMEAN
        @test gdd020.period == 20 * 365 * round(Int, CLM.SECSPDAY) ÷ 1800
    end

    # ------------------------------------------------------------------
    @testset "update_acc_vars drives GDD0/8/10 through the manager" begin
        np = 1
        temp = make_temp(np; t2m=TFRZ + 20.0)   # 20 C
        lun  = CLM.LandunitData{Float64}()
        pch  = CLM.PatchData{Float64}()
        mgr  = CLM.AccumManager()
        CLM.temperature_init_acc_buffer!(mgr, np; use_crop=true)

        bc = 1:1
        bp = 1:1
        # Three NH-summer days at 86400 s/step. GDD0 += 20/day (cap 26),
        # GDD8 += 12/day, GDD10 += 10/day.
        for nstep in 1:3
            CLM.temperature_update_acc_vars!(temp, bc, bp, lun, pch;
                secs=DTIME, dtime=DTIME, nstep=nstep,
                mgr=mgr, use_crop=true,
                month=7, day=15+nstep, jday=196+nstep,
                latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        end
        @test temp.gdd0_patch[1]  ≈ 60.0 atol=1e-12   # 3 * 20
        @test temp.gdd8_patch[1]  ≈ 36.0 atol=1e-12   # 3 * 12
        @test temp.gdd10_patch[1] ≈ 30.0 atol=1e-12   # 3 * 10
        # The always-on veg running mean also ran.
        @test temp.t_veg24_patch[1] ≈ TFRZ + 20.0 atol=1e-9
    end

    # ------------------------------------------------------------------
    @testset "no manager / use_crop=false skips GDD (still runs running means)" begin
        np = 1
        temp = make_temp(np; t2m=TFRZ + 20.0)
        lun  = CLM.LandunitData{Float64}()
        pch  = CLM.PatchData{Float64}()
        bc = 1:1; bp = 1:1
        CLM.temperature_update_acc_vars!(temp, bc, bp, lun, pch;
            secs=DTIME, dtime=DTIME, nstep=1)   # no mgr, use_crop=false
        @test temp.gdd0_patch[1] == 0.0          # untouched
        @test temp.t_veg24_patch[1] ≈ TFRZ + 20.0 atol=1e-9   # running mean ran
    end

    # ------------------------------------------------------------------
    @testset "end-of-year folds annual GDD into 20-yr runmean" begin
        np = 1
        temp = make_temp(np; t2m=TFRZ + 20.0)
        lun  = CLM.LandunitData{Float64}()
        pch  = CLM.PatchData{Float64}()
        mgr  = CLM.AccumManager()
        CLM.temperature_init_acc_buffer!(mgr, np; use_crop=true)
        CLM.varctl.flush_gdd20 = false
        bc = 1:1; bp = 1:1

        # Accumulate one summer day so gdd0=20, then trigger end-of-year.
        CLM.temperature_update_acc_vars!(temp, bc, bp, lun, pch;
            secs=DTIME, dtime=DTIME, nstep=1, mgr=mgr, use_crop=true,
            month=7, day=16, jday=197,
            latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        @test temp.gdd0_patch[1] ≈ 20.0 atol=1e-12

        # End of year: first sample of the 20-yr runmean == the annual value.
        CLM.temperature_update_acc_vars!(temp, bc, bp, lun, pch;
            secs=DTIME, dtime=DTIME, nstep=2, mgr=mgr, use_crop=true,
            month=12, day=31, jday=365, is_end_curr_year=true,
            latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        @test temp.gdd020_patch[1]  ≈ 20.0 atol=1e-12
        @test temp.gdd820_patch[1]  ≈ 12.0 atol=1e-12   # base 8 → 12/day
        @test temp.gdd1020_patch[1] ≈ 10.0 atol=1e-12   # base 10 → 10/day
    end

    # Run a full "year": accumulate `gdd0_target` worth of GDD0 in season
    # (one in-season day at the right t_ref2m), then end-of-year at Dec 31.
    # Returns nothing; mutates temp/mgr. `nstep_acc`/`nstep_eoy` are the step
    # indices to use for the in-season and year-end calls.
    function run_year!(temp, mgr, lun, pch, gdd0_target, nstep_acc, nstep_eoy)
        # GDD0 (base 0) accumulates min(26, t2m-TFRZ) per full day, so set
        # t_ref2m so that the daily increment equals gdd0_target.
        temp.t_ref2m_patch[1] = TFRZ + gdd0_target
        CLM.temperature_update_acc_vars!(temp, 1:1, 1:1, lun, pch;
            secs=DTIME, dtime=DTIME, nstep=nstep_acc, mgr=mgr, use_crop=true,
            month=7, day=16, jday=197,
            latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        # Year-end fold (NH Dec 31 is out of season → GDD0 unchanged this step).
        CLM.temperature_update_acc_vars!(temp, 1:1, 1:1, lun, pch;
            secs=DTIME, dtime=DTIME, nstep=nstep_eoy, mgr=mgr, use_crop=true,
            month=12, day=31, jday=365, is_end_curr_year=true,
            latdeg=latdeg, gridcell=gridcell, active=active, itype=itype)
        return nothing
    end

    # ------------------------------------------------------------------
    @testset "GDD20 runmean averages successive years" begin
        np = 1
        temp = make_temp(np)
        lun  = CLM.LandunitData{Float64}()
        pch  = CLM.PatchData{Float64}()
        mgr  = CLM.AccumManager()
        CLM.temperature_init_acc_buffer!(mgr, np; use_crop=true)
        CLM.varctl.flush_gdd20 = false

        # Year 1: gdd0 = 20 at year end (single in-season day).
        run_year!(temp, mgr, lun, pch, 20.0, 1, 2)
        @test temp.gdd0_patch[1]   ≈ 20.0 atol=1e-12
        @test temp.gdd020_patch[1] ≈ 20.0 atol=1e-12

        # Year 2: a second in-season day adds 10 more → annual gdd0 = 30
        # (runaccum persists across the out-of-season Dec 31; reset only on Jan 1).
        # runmean of {20, 30} after 2 samples = 25.
        run_year!(temp, mgr, lun, pch, 10.0, 3, 4)
        @test temp.gdd0_patch[1]   ≈ 30.0 atol=1e-12
        @test temp.gdd020_patch[1] ≈ 25.0 atol=1e-12
    end

    # ------------------------------------------------------------------
    @testset "flush_gdd20 resets the 20-yr runmean once" begin
        np = 1
        temp = make_temp(np)
        lun  = CLM.LandunitData{Float64}()
        pch  = CLM.PatchData{Float64}()
        mgr  = CLM.AccumManager()
        CLM.temperature_init_acc_buffer!(mgr, np; use_crop=true)
        CLM.varctl.flush_gdd20 = false

        # Year 1: seed a prior 20-yr mean from annual gdd0 = 25.
        run_year!(temp, mgr, lun, pch, 25.0, 1, 2)
        @test temp.gdd020_patch[1] ≈ 25.0 atol=1e-12

        # Year 2 with a flush requested: the year-end resets the runmean, then
        # folds in only the new annual gdd0 (25+5 = 30 because runaccum persists)
        # → runmean == 30, not (25+30)/2.
        CLM.varctl.flush_gdd20 = true
        run_year!(temp, mgr, lun, pch, 5.0, 3, 4)
        @test temp.gdd0_patch[1]   ≈ 30.0 atol=1e-12
        @test temp.gdd020_patch[1] ≈ 30.0 atol=1e-12
        @test CLM.varctl.flush_gdd20 == false   # consumed
    end

end
