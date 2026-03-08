# ==========================================================================
# Integration test for clm_run! end-to-end simulation pipeline
# ==========================================================================

include("generate_forcing.jl")

@testset "CLM End-to-End Simulation" begin

    # Input data paths
    fsurdat = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
    paramfile = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"

    # Skip if input files not available
    if !isfile(fsurdat) || !isfile(paramfile)
        @warn "Skipping integration test: input files not found"
        @test true  # placeholder pass
        return
    end

    # ---- Test 1: Generate synthetic forcing file ----
    @testset "Synthetic forcing generation" begin
        forcing_path = tempname() * "_forcing.nc"
        generate_forcing(forcing_path;
            start_date=DateTime(2000, 1, 1),
            end_date=DateTime(2000, 1, 2),
            dtime=1800)
        @test isfile(forcing_path)

        # Verify forcing file structure
        ds = NCDataset(forcing_path, "r")
        @test haskey(ds, "time")
        @test haskey(ds, "TBOT")
        @test haskey(ds, "PSRF")
        @test haskey(ds, "WIND")
        @test haskey(ds, "FLDS")
        @test haskey(ds, "FSDS")
        @test haskey(ds, "PRECTmms")
        @test haskey(ds, "QBOT")
        ntimes = length(ds["time"])
        @test ntimes == 48  # 24 hours / 0.5 hour = 48 steps
        @test all(ds["TBOT"][:] .> 250.0)  # reasonable temperatures
        @test all(ds["TBOT"][:] .< 290.0)
        @test all(ds["PSRF"][:] .> 80000.0)  # reasonable pressure
        @test all(ds["FSDS"][:] .>= 0.0)    # non-negative shortwave
        close(ds)
        rm(forcing_path)
    end

    # ---- Test 2: Forcing reader ----
    @testset "Forcing reader" begin
        forcing_path = tempname() * "_forcing.nc"
        generate_forcing(forcing_path;
            start_date=DateTime(2000, 1, 1),
            end_date=DateTime(2000, 1, 2),
            dtime=1800)

        fr = CLM.ForcingReader()
        CLM.forcing_reader_init!(fr, forcing_path)
        @test fr.ntimes == 48
        @test length(fr.times) == 48

        # Read a step
        a2l = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a2l, 1, 1, 1)
        CLM.read_forcing_step!(fr, a2l, DateTime(2000, 1, 1, 12, 0), 1, 1)

        @test a2l.forc_t_not_downscaled_grc[1] > 250.0
        @test a2l.forc_t_not_downscaled_grc[1] < 290.0
        @test a2l.forc_pbot_not_downscaled_grc[1] > 80000.0
        @test a2l.forc_wind_grc[1] > 0.0
        @test a2l.forc_th_not_downscaled_grc[1] > 0.0  # potential temp computed
        @test a2l.forc_rho_not_downscaled_grc[1] > 0.0  # density computed

        CLM.forcing_reader_close!(fr)
        rm(forcing_path)
    end

    # ---- Test 3: Orbital parameters ----
    @testset "Orbital parameters" begin
        # Summer solstice: calday ~172 (June 21)
        (declin_summer, eccf_summer) = CLM.compute_orbital(172.0)
        @test declin_summer > 0.0    # Northern hemisphere summer → positive declination
        @test declin_summer < 0.5    # Less than ~28 degrees
        @test eccf_summer > 0.9      # Earth-sun distance factor ~1

        # Winter solstice: calday ~356 (Dec 22)
        (declin_winter, eccf_winter) = CLM.compute_orbital(356.0)
        @test declin_winter < 0.0    # Northern hemisphere winter → negative declination
        @test eccf_winter > 0.9

        # Equinox: calday ~80 (March 21)
        (declin_eq, _) = CLM.compute_orbital(80.5)
        @test abs(declin_eq) < 0.02  # ~0 at equinox
    end

    # ---- Test 4: Downscale forcings (trivial case) ----
    @testset "Downscale forcings" begin
        # Single gridcell: downscaling should basically copy values
        (inst, bounds, filt, tm) = CLM.clm_initialize!(;
            fsurdat=fsurdat, paramfile=paramfile)

        ng = bounds.endg; nc = bounds.endc
        a2l = inst.atm2lnd

        # Set gridcell forcings
        for g in 1:ng
            a2l.forc_t_not_downscaled_grc[g] = 270.0
            a2l.forc_pbot_not_downscaled_grc[g] = 85000.0
            a2l.forc_th_not_downscaled_grc[g] = 275.0
            a2l.forc_rho_not_downscaled_grc[g] = 1.1
            a2l.forc_lwrad_not_downscaled_grc[g] = 250.0
            a2l.forc_vp_grc[g] = 300.0
            a2l.forc_hgt_grc[g] = 30.0
            a2l.forc_topo_grc[g] = 0.0
            for b in 1:CLM.NUMRAD
                a2l.forc_solad_not_downscaled_grc[g, b] = 100.0
            end
            a2l.forc_solar_not_downscaled_grc[g] = 200.0
            a2l.forc_rain_not_downscaled_grc[g] = 0.001
            a2l.forc_snow_not_downscaled_grc[g] = 0.0
        end

        CLM.downscale_forcings!(bounds, a2l, inst.column, inst.landunit, inst.topo)

        # For flat terrain (topo=0), downscaled should match gridcell
        @test a2l.forc_t_downscaled_col[1] ≈ 270.0 atol=1.0
        @test a2l.forc_pbot_downscaled_col[1] ≈ 85000.0 atol=100.0
        @test a2l.forc_lwrad_downscaled_col[1] ≈ 250.0 atol=1.0
        # Rain/snow should be partitioned
        @test a2l.forc_rain_downscaled_col[1] >= 0.0
        @test a2l.forc_snow_downscaled_col[1] >= 0.0
        total_precip = a2l.forc_rain_downscaled_col[1] + a2l.forc_snow_downscaled_col[1]
        @test total_precip ≈ 0.001 atol=1e-6
    end

    # ---- Test 5: History writer ----
    @testset "History writer" begin
        hist_path = tempname() * "_history.nc"
        fields = CLM.default_hist_fields()
        hw = CLM.HistoryWriter()
        CLM.history_writer_init!(hw, hist_path, fields, 1, 5, 10)

        @test isfile(hist_path)
        CLM.history_writer_close!(hw)

        ds = NCDataset(hist_path, "r")
        @test haskey(ds, "time")
        @test haskey(ds, "T_GRND")
        @test haskey(ds, "QFLX_EVAP_TOT")
        @test haskey(ds, "FSA")
        @test haskey(ds, "H2OSNO")
        @test haskey(ds, "SNOWDP")
        close(ds)
        rm(hist_path)
    end

    # ---- Test 5b: Vegetation masks for ZWT/BTRAN history diagnostics ----
    @testset "History vegetation masks (ZWT/BTRAN)" begin
        (inst, bounds, _, _) = CLM.clm_initialize!(;
            fsurdat=fsurdat, paramfile=paramfile)

        np = bounds.endp
        nc = bounds.endc
        @test np >= 1
        @test nc >= 1

        # Force all patches to bare-ground, then make one vegetated patch.
        inst.patch.itype .= CLM.noveg
        inst.patch.itype[1] = 1

        # Use running daily-min source path for BTRAN (preferred over completed-day min).
        inst.energyflux.btran_min_patch .= 0.77
        inst.energyflux.btran_min_inst_patch .= 0.66
        inst.energyflux.btran_patch .= 0.55
        btran_hist = CLM.history_btran_daily_min_patch(inst)
        @test btran_hist[1] ≈ 0.66
        for p in 2:np
            @test isnan(btran_hist[p])
        end

        # ZWT should be retained only for columns that host at least one vegetated patch.
        for c in 1:nc
            inst.soilhydrology.zwt_col[c] = 1.0 + c
        end
        zwt_hist = CLM.history_zwt_veg_col(inst)
        veg_col = inst.patch.column[1]
        @test zwt_hist[veg_col] ≈ inst.soilhydrology.zwt_col[veg_col]
        for c in 1:nc
            if c != veg_col
                @test isnan(zwt_hist[c])
            end
        end
    end

    # ---- Test 6: Full end-to-end simulation (1 day, 48 timesteps) ----
    @testset "Full simulation (1 day)" begin
        # Generate forcing
        forcing_path = tempname() * "_forcing.nc"
        hist_path = tempname() * "_history.nc"

        generate_forcing(forcing_path;
            start_date=DateTime(2000, 1, 1),
            end_date=DateTime(2000, 1, 2),
            dtime=1800)

        # Run simulation
        inst = CLM.clm_run!(;
            fsurdat=fsurdat,
            paramfile=paramfile,
            fforcing=forcing_path,
            fhistory=hist_path,
            start_date=DateTime(2000, 1, 1),
            end_date=DateTime(2000, 1, 2),
            dtime=1800,
            use_cn=false,
            verbose=false)

        @test inst isa CLM.CLMInstances

        # Verify history output
        @test isfile(hist_path)
        ds = NCDataset(hist_path, "r")
        @test haskey(ds, "time")
        ntimes = length(ds["time"])
        @test ntimes == 1  # daily averaging: 48 half-hour steps → 1 daily record

        # Check T_GRND: should be physically reasonable (wide range for cold-start transient)
        # 0.0 is a valid default for inactive columns
        tgrnd = ds["T_GRND"][:, :]
        @test all(x -> x == -9999.0 || x == 0.0 || (100.0 < x < 400.0), tgrnd)

        # Check at least some non-fill values exist
        @test any(x -> x != -9999.0 && x != 0.0, tgrnd)

        # Check H2OSNO: non-negative
        h2osno = ds["H2OSNO"][:, :]
        @test all(x -> x == -9999.0 || x >= 0.0, h2osno)

        close(ds)

        # Cleanup
        rm(forcing_path)
        rm(hist_path)
    end

    # ---- Test 7: No-aquifer mode smoke test ----
    @testset "Full simulation (no aquifer, smoke)" begin
        forcing_path = tempname() * "_forcing.nc"
        hist_path = tempname() * "_history.nc"

        generate_forcing(forcing_path;
            start_date=DateTime(2000, 1, 1),
            end_date=DateTime(2000, 1, 1, 1),
            dtime=1800)

        inst = CLM.clm_run!(;
            fsurdat=fsurdat,
            paramfile=paramfile,
            fforcing=forcing_path,
            fhistory=hist_path,
            start_date=DateTime(2000, 1, 1),
            end_date=DateTime(2000, 1, 1, 1),
            dtime=1800,
            use_cn=false,
            use_aquifer_layer=false,
            verbose=false)

        @test inst isa CLM.CLMInstances
        @test isfile(hist_path)

        ds = NCDataset(hist_path, "r")
        @test length(ds["time"]) == 1  # partial day flushed on close → 1 record
        tgrnd = ds["T_GRND"][:, :]
        @test all(isfinite, tgrnd[tgrnd .!= -9999.0])
        close(ds)

        rm(forcing_path)
        rm(hist_path)
    end
end
