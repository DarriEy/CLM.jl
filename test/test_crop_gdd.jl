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

end
