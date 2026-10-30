# ==========================================================================
# SPITFIRE fire modes 2-5 (FATESFireFactoryMod: 2 lightning_from_data,
# 3 successful_ignitions, 4 anthro_ignitions, 5 anthro_suppression).
#
#   1. hdm stream year window (popd_streams stream_year_first/last_popdens,
#      model_year_align_popdens; CDEPS taxmode 'extend' + tintalgo 'nearest').
#   2. Per-mode ignition logic of area_burnt_intensity! (SFMainMod.F90:499-535).
#   3. Data-gated, real CTSM streams + committed Fortran extracts
#      (scripts/validation/fates_fortran_parity/fates_spitfire_m{2,3,4}real_thr5_fortran.csv,
#      Bow FATES case, fates_fire_threshold 5 kW/m): host LNFM/HDM, and the per-mode
#      FATES_IGNITIONS ratio on the first burn day.
# ==========================================================================
using Test, NCDatasets
include("testdata.jl")

# Save/restore the SPITFIRE globals this file touches.
function _with_sf_globals(f)
    saved = (CLM.SFParams[], CLM.EDParams[], CLM.hlm_masterproc[], CLM.hlm_spitfire_mode[],
             CLM.hlm_sf_nofire_def[], CLM.hlm_sf_scalar_lightning_def[],
             CLM.hlm_sf_successful_ignitions_def[], CLM.hlm_sf_anthro_ignitions_def[])
    try
        f()
    finally
        (CLM.SFParams[], CLM.EDParams[], CLM.hlm_masterproc[], CLM.hlm_spitfire_mode[],
         CLM.hlm_sf_nofire_def[], CLM.hlm_sf_scalar_lightning_def[],
         CLM.hlm_sf_successful_ignitions_def[], CLM.hlm_sf_anthro_ignitions_def[]) = saved
    end
end

# A one-patch site that burns (FI far above threshold), default FATES fire parameters
# (fates_params_default.cdl), CTSM's mode enumeration (clmfates_interfaceMod.F90:464-467).
function _sf_mode_site(fdi::Real)
    sfp = CLM.sf_params_type()
    sfp.SF_val_fdi_alpha      = 0.00037
    sfp.SF_val_fuel_energy    = 18000.0
    sfp.SF_val_max_durat      = 240.0
    sfp.SF_val_durat_slope    = -11.06
    sfp.SF_val_fire_threshold = 5.0
    CLM.SFParams[] = sfp
    edp = CLM.ed_params_type()
    edp.ED_val_nignitions = 15.0
    edp.cg_strikes        = 0.2
    CLM.EDParams[] = edp
    CLM.hlm_masterproc[]                  = CLM.ifalse
    CLM.hlm_sf_nofire_def[]               = 0
    CLM.hlm_sf_scalar_lightning_def[]     = 1
    CLM.hlm_sf_successful_ignitions_def[] = 3
    CLM.hlm_sf_anthro_ignitions_def[]     = 4

    p = CLM.fates_patch_type()
    p.area = CLM.area; p.total_tree_area = 0.3 * CLM.area
    p.nocomp_pft_label = 1; p.patchno = 1
    p.ros_front = 2.0; p.ros_back = 0.5; p.tfc_ros = 0.5
    site = CLM.ed_site_type()
    site.oldest_patch = p; site.youngest_patch = p
    nw = CLM.nesterov_index()
    nw.fire_weather_index  = -log(1 - fdi) / sfp.SF_val_fdi_alpha   # FDI = 1 - exp(-alpha NI)
    nw.effective_windspeed = 30.0
    site.fireWeather = nw
    return site, p
end

# NF and NF_successful for one mode, with lightning24 [#/km2/day] and pop_density.
function _sf_mode_nf(mode, fdi, lightning24, pop_density)
    site, p = _sf_mode_site(fdi)
    CLM.hlm_spitfire_mode[] = mode
    bc = CLM.bc_in_type()
    bc.lightning24 = [lightning24]; bc.pop_density = [pop_density]
    CLM.area_burnt_intensity!(site, bc)
    @test p.fire == 1
    return site.NF, site.fdi, site.NF_successful
end

@testset "SPITFIRE modes: hdm stream year window (CDEPS extend + nearest)" begin
    dir = mktempdir()
    fh = joinpath(dir, "hdm_multiyear.nc")
    years = 1850:1855
    NCDataset(fh, "c") do ds
        ds.dim["lon"] = 2; ds.dim["lat"] = 2; ds.dim["time"] = length(years)
        defVar(ds, "lon", [0.0, 180.0], ("lon",))
        defVar(ds, "lat", [-45.0, 45.0], ("lat",))
        # HYDE layout: one mid-year record per year, noleap days since 1850-01-01
        defVar(ds, "time", [182.5 + 365.0 * k for k in 0:length(years)-1], ("time",),
               attrib = ["units" => "days since 1850-01-01", "calendar" => "noleap"])
        defVar(ds, "hdm", Float64[y for _ in 1:2, _ in 1:2, y in years], ("lon", "lat", "time"),
               attrib = ["units" => "counts/km^2"])
    end
    grc = CLM.GridcellData(); grc.londeg = [10.0]; grc.latdeg = [0.0]
    hdm_at(d; year, calday) = (CLM.fates_fire_interp!(d, grc; calday = calday, dayspyr = 365.0,
                                                      bounds_grc = 1:1, year = year); d.forc_hdm[1])

    # No window: nearest by day of year over all records -> record 1 (all share a doy).
    d0 = CLM.fates_fire_data_init!(CLM.FATESFireData(), 1; fhdm = fh)
    @test d0.stream.hdm.year_first == 0
    @test hdm_at(d0; year = 2002, calday = 100.0) ≈ 1850.0
    @test d0.tod_offset == 0                      # streams at the driver's (jday, secs)

    # The CTSM setting for a present-day run: first = last = 2000, align 1. Model
    # year maps past the window (dYear = 2000 + 2002 - 1) -> the GVD record.
    d1 = CLM.fates_fire_data_init!(CLM.FATESFireData(), 1; fhdm = fh,
        stream_year_first_popdens = 1852, stream_year_last_popdens = 1852,
        model_year_align_popdens = 1)
    for (y, cd) in ((2002, 1.0), (2002, 200.0), (1, 1.0), (3000, 365.9))
        @test hdm_at(d1; year = y, calday = cd) ≈ 1852.0
    end
    # align defaults to stream_year_first_popdens
    d1b = CLM.fates_fire_data_init!(CLM.FATESFireData(), 1; fhdm = fh,
        stream_year_first_popdens = 1853, stream_year_last_popdens = 1853)
    @test d1b.stream.hdm.year_align == 1853

    # Multi-year window 1851-1854 aligned to model year 2001.
    d2 = CLM.fates_fire_data_init!(CLM.FATESFireData(), 1; fhdm = fh,
        stream_year_first_popdens = 1851, stream_year_last_popdens = 1854,
        model_year_align_popdens = 2001)
    @test hdm_at(d2; year = 2000, calday = 300.0) ≈ 1851.0   # before the LVD -> LVD
    @test hdm_at(d2; year = 2001, calday = 100.0) ≈ 1851.0   # before the LVD -> LVD
    @test hdm_at(d2; year = 2001, calday = 200.0) ≈ 1851.0   # nearer 1851-07-02
    @test hdm_at(d2; year = 2002, calday = 1.0)   ≈ 1851.0   # exact tie -> lower bound
    @test hdm_at(d2; year = 2002, calday = 2.0)   == 1852.0   # nearer the upper bound
    @test hdm_at(d2; year = 2003, calday = 250.0) ≈ 1853.0
    @test hdm_at(d2; year = 2004, calday = 250.0) ≈ 1854.0   # at/after the GVD -> GVD
    @test hdm_at(d2; year = 2040, calday = 10.0)  == 1854.0
    # 1855 is outside the window and is never used
    @test hdm_at(d2; year = 2005, calday = 300.0) ≈ 1854.0

    @test_throws ErrorException CLM.fates_fire_data_init!(CLM.FATESFireData(), 1; fhdm = fh,
        stream_year_first_popdens = 1900, stream_year_last_popdens = 1900)
end

@testset "SPITFIRE modes: ignition logic per mode (SFMainMod area_burnt_intensity)" begin
    _with_sf_globals() do
        fdi = 0.8; L = 0.03; H = 25.0
        cg = 0.2
        anthro = 0.0035 * 6.8 * H^0.43 / 30.0          # Li et al. 2012, SFMainMod.F90:532
        nf1, f1, s1 = _sf_mode_nf(1, fdi, L, H)
        @test nf1 ≈ 15.0 / 365.0 * cg                  # scalar: ED_val_nignitions/yr
        @test f1 ≈ fdi
        nf2, f2, s2 = _sf_mode_nf(2, fdi, L, H)        # lightning_from_data
        @test nf2 ≈ L * cg
        @test f2 ≈ fdi
        @test s2 ≈ nf2 * f2                            # one patch covering the site
        nf3, f3, s3 = _sf_mode_nf(3, fdi, L, H)        # successful_ignitions
        @test f3 == 1.0                                # FDI forced to 1
        @test nf3 ≈ L                                  # cloud_to_ground = 1
        nf4, f4, s4 = _sf_mode_nf(4, fdi, L, H)        # anthro_ignitions
        @test nf4 ≈ L * cg + anthro
        @test f4 ≈ fdi
        nf5, f5, s5 = _sf_mode_nf(5, fdi, L, H)        # anthro_suppression: no FATES
        @test (nf5, f5, s5) == (nf2, f2, s2)           # branch, behaves as mode 2
        @test s3 / s2 ≈ 1 / (cg * fdi)
        @test s4 / s2 ≈ 1 + anthro / (L * cg)
    end
end

# ---------------------------------------------------------------------------
# Data-gated: real CTSM fire streams + committed Fortran extracts.
# ---------------------------------------------------------------------------
const _SF_VDIR = joinpath(@__DIR__, "..", "scripts", "validation", "fates_fortran_parity")
_sf_extract(m) = joinpath(_SF_VDIR, "fates_spitfire_m$(m)real_thr5_fortran.csv")
function _sf_read_csv(path)
    l = readlines(path); h = split(l[1], ',')
    rows = [parse.(Float64, split(x, ',')) for x in l[2:end] if !isempty(x)]
    return Dict(String(h[i]) => [r[i] for r in rows] for i in eachindex(h))
end
function _sf_real_streams()
    for root in (cesm_inputdata_root(), symfluence_path("installs", "cesm-inputdata"))
        fl = joinpath(root, "atm", "datm7", "NASA_LIS",
                      "clmforc.Li_2012_climo1995-2011.T62.lnfm_Total_c140423.nc")
        fh = joinpath(root, "lnd", "clm2", "firedata",
                      "clmforc.Li_2017_HYDEv3.2_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2016_c180202.nc")
        isfile(fl) && isfile(fh) && return (fl, fh)
    end
    return nothing
end

@testset "SPITFIRE modes: real streams vs Fortran (data-gated)" begin
    streams = _sf_real_streams()
    if streams === nothing || !isfile(_sf_extract(2))
        @test_skip "real LIS lightning / HYDE popdens streams not staged"
    else
        flnfm, fhdm = streams
        F = Dict(m => _sf_read_csv(_sf_extract(m)) for m in 2:4)
        # Fortran site (fates_parity_1pt h0 lon/lat)
        grc = CLM.GridcellData(); grc.londeg = [-116.0118]; grc.latdeg = [51.36051]
        kw = (stream_year_first_popdens = 2000, stream_year_last_popdens = 2000,
              model_year_align_popdens = 1)
        d = redirect_stderr(devnull) do      # the files' unit strings differ from CTSM's
            CLM.fates_fire_data_init!(CLM.FATESFireData(), 1; flnfm, fhdm, kw...)
        end
        d.tod_offset = 3600                  # step-start (jday, secs) as in the harness

        # Step the host fire data hourly from 2002-07-01 00:00 (noleap), as the
        # driver does: FireInterp at the step's (start) time + tod_offset, then
        # UpdateAccVars. Record day means of forc_lnfm and lnfm24 at each day start.
        ndays = 401
        lnfm_day = zeros(ndays); lnfm24_daystart = zeros(ndays + 1)
        jd0 = 182                                       # 1 July, noleap
        for k in 0:ndays*24-1
            (k % 24 == 0) && (lnfm24_daystart[k ÷ 24 + 1] = d.lnfm24[1])
            doy = mod(jd0 - 1 + k ÷ 24, 365) + 1
            yr  = 2002 + (jd0 - 1 + k ÷ 24) ÷ 365
            CLM.fates_fire_interp!(d, grc; calday = doy + (3600 * (k % 24) + d.tod_offset) / 86400,
                                   dayspyr = 365.0, bounds_grc = 1:1, year = yr)
            CLM.fates_fire_update_acc_vars!(d, 1:1; dtime = 3600)
            lnfm_day[k ÷ 24 + 1] += d.forc_lnfm[1] / 24
        end
        lnfm24_daystart[end] = d.lnfm24[1]

        # Host streams vs Fortran LNFM / HDM history (identical in every mode > 1)
        FL = F[2]["LNFM"][1:ndays]
        rel = abs.(lnfm_day .- FL) ./ max.(abs.(FL), 1e-3 * maximum(FL))
        @test sum(rel) / ndays < 2e-4
        @test maximum(rel) < 5e-3
        @test all(F[m]["LNFM"] == F[2]["LNFM"] && F[m]["HDM"] == F[2]["HDM"] for m in 3:4)
        @test d.forc_hdm[1] ≈ F[2]["HDM"][1] rtol = 1e-3   # 1.399 /km2 (year 2000)
        # Without the year window the 1850 record (0.02 /km2) would be used.
        d0 = redirect_stderr(devnull) do
            CLM.fates_fire_data_init!(CLM.FATESFireData(), 1; fhdm)
        end
        CLM.fates_fire_interp!(d0, grc; calday = 1.0, dayspyr = 365.0, bounds_grc = 1:1, year = 2002)
        @test d0.forc_hdm[1] < 0.05

        # Per mode, first burn day (402): the stand is still the same in every mode
        # and the burn decision (FI > threshold) does not depend on NF or FDI, so
        # FATES_IGNITIONS(mode)/FATES_IGNITIONS(mode 2) = NF*FDI ratio. Drive
        # area_burnt_intensity! with the host lightning24 = lnfm24*24 at the day
        # start (the daily pack, clmfates_interfaceMod.F90:1099-1106) and FDI from
        # the Fortran mode-2 history.
        dburn = findfirst(>(0), F[2]["FATES_BURNFRAC"])
        @test dburn == 402
        @test all(findfirst(>(0), F[m]["FATES_BURNFRAC"]) == dburn for m in 3:4)
        L24 = lnfm24_daystart[dburn] * 24
        fdi = F[2]["FATES_FDI"][dburn]
        _with_sf_globals() do
            ign = Dict{Int,Float64}()
            for m in 2:5
                _, f, s = _sf_mode_nf(m, fdi, L24, m >= 4 ? d.forc_hdm[1] : 0.0)
                ign[m] = s
                m == 3 ? (@test f == 1.0) : (@test f ≈ fdi)
            end
            r3F = F[3]["FATES_IGNITIONS"][dburn] / F[2]["FATES_IGNITIONS"][dburn]
            r4F = F[4]["FATES_IGNITIONS"][dburn] / F[2]["FATES_IGNITIONS"][dburn]
            @test ign[3] / ign[2] ≈ r3F rtol = 1e-5      # 6.1977 = 1/(cg FDI)
            @test ign[4] / ign[2] ≈ r4F rtol = 2e-3      # 1.1598: anthro vs lightning
            @test ign[5] == ign[2]                         # Fortran m5 == m2 bitwise
        end
    end
end
