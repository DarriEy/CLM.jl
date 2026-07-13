# ==========================================================================
# Tests for the atmospheric N-deposition stream (infrastructure/ndep_stream.jl)
# and the N-input wiring it feeds (n_deposition!, cn_lag_npp_update!).
#
# The stream is the ONLY path by which nitrogen enters the CLM ecosystem.
# Before this was wired, forc_ndep_grc was initialised to 0.0 and never
# assigned, so ndep_to_sminn was identically zero.
# ==========================================================================

using Test
using NCDatasets
using CLM

# Build a small synthetic ndep file. `vals` is (nlon, nlat, ntime) in `units`.
function _write_ndep_file(path, lon, lat, tdays, vals, units)
    isfile(path) && rm(path)
    ds = NCDataset(path, "c")
    defDim(ds, "lon", length(lon)); defDim(ds, "lat", length(lat))
    defDim(ds, "time", length(tdays))
    defVar(ds, "lon", collect(Float64, lon), ("lon",))
    defVar(ds, "lat", collect(Float64, lat), ("lat",))
    defVar(ds, "time", collect(Float64, tdays), ("time",);
           attrib = ["units" => "days since 2000-01-01 00:00:00",
                     "calendar" => "gregorian"])
    defVar(ds, "NDEP_month", Array{Float64}(vals), ("lon", "lat", "time");
           attrib = ["units" => units])
    close(ds)
    return path
end

@testset "ndep stream" begin
    tmp = mktempdir()
    # 12 monthly mid-month day-of-year samples on a coarse global grid.
    lon = 0.0:90.0:270.0            # 4
    lat = -60.0:60.0:60.0           # 3
    tt  = [15.5, 45.0, 74.5, 105.0, 135.5, 166.0,
           196.5, 227.5, 258.0, 288.5, 319.0, 349.5] .- 1.0

    @testset "units handling" begin
        # g(N)/m2/s -> no conversion
        f = _write_ndep_file(joinpath(tmp, "s.nc"), lon, lat, tt,
                             fill(2.0e-8, 4, 3, 12), "g(N)/m2/s")
        s = CLM.NDepStreamData()
        CLM.ndep_stream_init!(s, f)
        @test s.active
        @test s.divide_by_secs_per_yr == false

        # g(N)/m2/yr -> divide by secspday*dayspyr
        f2 = _write_ndep_file(joinpath(tmp, "y.nc"), lon, lat, tt,
                              fill(0.5, 4, 3, 12), "g(N)/m2/yr")
        s2 = CLM.NDepStreamData()
        CLM.ndep_stream_init!(s2, f2)
        @test s2.divide_by_secs_per_yr == true

        # Anything else must ERROR, exactly as CTSM's check_units does. A silently
        # mis-scaled N input is the bug class this module exists to remove.
        f3 = _write_ndep_file(joinpath(tmp, "bad.nc"), lon, lat, tt,
                              fill(1.0, 4, 3, 12), "kg/m2/day")
        @test_throws ErrorException CLM.ndep_stream_init!(CLM.NDepStreamData(), f3)

        # Missing file / missing variable are errors too.
        @test_throws ErrorException CLM.ndep_stream_init!(CLM.NDepStreamData(),
                                                          joinpath(tmp, "nope.nc"))
    end

    @testset "uniform + constant stream is reproduced EXACTLY" begin
        # This is the reduction the Fortran parity reference relies on: for a
        # spatially uniform, time-constant stream, bilinear weights sum to 1 and the
        # linear time weights sum to 1, so the interpolated value is the stream value
        # exactly -- in BOTH codes. That is what lets the parity test isolate the
        # N-cycle physics from forcing-interpolation fidelity.
        val = 1.0e-8
        f = _write_ndep_file(joinpath(tmp, "u.nc"), lon, lat, tt,
                             fill(val, 4, 3, 12), "g(N)/m2/s")
        s = CLM.NDepStreamData(); CLM.ndep_stream_init!(s, f)

        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, 3)
        grc.latdeg .= [51.36, -12.0, 33.3]
        grc.londeg .= [244.0, 10.0, 355.0]

        a2l = CLM.Atm2LndData(); CLM.atm2lnd_init!(a2l, 3, 3, 3)

        for calday in (1.0, 45.7, 200.3, 288.5, 364.9)
            CLM.ndep_interp!(a2l, s, grc; calday = calday, dayspyr = 365.0,
                             bounds_grc = 1:3)
            for g in 1:3
                @test a2l.forc_ndep_grc[g] ≈ val rtol=1e-14
            end
        end
    end

    @testset "g(N)/m2/yr converts to g(N)/m2/s" begin
        yr = 0.5   # gN/m2/yr
        f = _write_ndep_file(joinpath(tmp, "y2.nc"), lon, lat, tt,
                             fill(yr, 4, 3, 12), "g(N)/m2/yr")
        s = CLM.NDepStreamData(); CLM.ndep_stream_init!(s, f)

        grc = CLM.GridcellData(); CLM.gridcell_init!(grc, 1)
        grc.latdeg .= [51.36]; grc.londeg .= [244.0]
        a2l = CLM.Atm2LndData(); CLM.atm2lnd_init!(a2l, 1, 1, 1)

        CLM.ndep_interp!(a2l, s, grc; calday = 100.0, dayspyr = 365.0, bounds_grc = 1:1)
        @test a2l.forc_ndep_grc[1] ≈ yr / (86400.0 * 365.0) rtol=1e-12
    end

    @testset "bilinear is exact on a field linear in lon/lat" begin
        # Bilinear interpolation reproduces a bilinear field exactly. Use a field
        # linear in latitude only (longitude is periodic, so a globally linear lon
        # field is not well defined).
        vals = zeros(4, 3, 12)
        for i in 1:4, j in 1:3, k in 1:12
            vals[i, j, k] = 1.0e-8 + 1.0e-10 * lat[j]
        end
        f = _write_ndep_file(joinpath(tmp, "lin.nc"), lon, lat, tt, vals, "g(N)/m2/s")
        s = CLM.NDepStreamData(); CLM.ndep_stream_init!(s, f)

        grc = CLM.GridcellData(); CLM.gridcell_init!(grc, 1)
        a2l = CLM.Atm2LndData(); CLM.atm2lnd_init!(a2l, 1, 1, 1)
        for tlat in (-30.0, 0.0, 17.5, 45.0)
            grc.latdeg .= [tlat]; grc.londeg .= [123.0]
            CLM.ndep_interp!(a2l, s, grc; calday = 150.0, dayspyr = 365.0, bounds_grc = 1:1)
            @test a2l.forc_ndep_grc[1] ≈ 1.0e-8 + 1.0e-10 * tlat rtol=1e-12
        end
    end

    @testset "inactive stream leaves forc_ndep untouched (default path)" begin
        s = CLM.NDepStreamData()            # never initialised
        grc = CLM.GridcellData(); CLM.gridcell_init!(grc, 1)
        grc.latdeg .= [51.0]; grc.londeg .= [244.0]
        a2l = CLM.Atm2LndData(); CLM.atm2lnd_init!(a2l, 1, 1, 1)
        a2l.forc_ndep_grc[1] = 7.0
        CLM.ndep_interp!(a2l, s, grc; calday = 1.0, dayspyr = 365.0, bounds_grc = 1:1)
        @test a2l.forc_ndep_grc[1] == 7.0   # untouched => default run is unchanged
    end
end

@testset "n_deposition! gathers forc_ndep to the column" begin
    nc = 3
    nf = CLM.SoilBiogeochemNitrogenFluxData()
    CLM.soil_bgc_nitrogen_flux_init!(nf, nc, 1, 7, 7)

    forc_ndep = [1.0e-8, 5.0e-9]        # 2 gridcells
    col_gridcell = [1, 2, 2]

    CLM.n_deposition!(nf; forc_ndep = forc_ndep, col_gridcell = col_gridcell,
                      bounds = 1:nc)

    # CNNDeposition is a pure gridcell->column gather: ALL atmospheric N deposition
    # goes to the soil mineral N pool.
    @test nf.ndep_to_sminn_col[1] == 1.0e-8
    @test nf.ndep_to_sminn_col[2] == 5.0e-9
    @test nf.ndep_to_sminn_col[3] == 5.0e-9
end

@testset "cn_lag_npp_update! (feeds CNNFixation)" begin
    nc = 2
    cf = CLM.CNVegCarbonFluxData()
    CLM.cnveg_carbon_flux_init!(cf, 2, nc, 1; nlevdecomp_full = 1, ndecomp_pools = 7)
    mask = trues(nc)

    cf.npp_col .= [1.0e-6, 2.0e-6]
    # Fresh state is SPVAL => Fortran's "first timestep" branch seeds lag = npp.
    @test all(cf.lag_npp_col .== CLM.SPVAL)
    CLM.cn_lag_npp_update!(cf, mask; dt = 3600.0, nfix_timeconst = 10.0)
    @test cf.lag_npp_col ≈ [1.0e-6, 2.0e-6]

    # Subsequent steps relax exponentially toward npp with tau = nfix_timeconst days.
    cf.npp_col .= [0.0, 0.0]
    decay = exp(-3600.0 / (10.0 * CLM.SECSPDAY))
    CLM.cn_lag_npp_update!(cf, mask; dt = 3600.0, nfix_timeconst = 10.0)
    @test cf.lag_npp_col ≈ [1.0e-6 * decay, 2.0e-6 * decay] rtol=1e-12

    # Outside the (0, 500) window Fortran does not maintain lag_npp at all.
    before = copy(cf.lag_npp_col)
    CLM.cn_lag_npp_update!(cf, mask; dt = 3600.0, nfix_timeconst = 0.0)
    @test cf.lag_npp_col == before
end
