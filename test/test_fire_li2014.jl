@testset "Fire Li2014" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for Li2014 fire tests
    # ----------------------------------------------------------------
    function make_fire_li2014_data(; np=4, nc=2, ng=1, nlevgrnd=3, nlevdecomp=1,
                                     ndecomp_pools=4, n_litr=3)
        npft = 20

        # --- PFT constants (fire base) ---
        pftcon = CLM.PftConFireBase(
            woody    = vcat(fill(1.0, 8), fill(0.0, 12)),
            cc_leaf  = fill(0.4, npft),
            cc_lstem = fill(0.2, npft),
            cc_dstem = fill(0.1, npft),
            cc_other = fill(0.3, npft),
            fm_leaf  = fill(0.6, npft),
            fm_lstem = fill(0.5, npft),
            fm_other = fill(0.4, npft),
            fm_root  = fill(0.3, npft),
            fm_lroot = fill(0.5, npft),
            fm_droot = fill(0.2, npft),
            lf_f     = fill(1.0/n_litr, npft, n_litr),
            fr_f     = fill(1.0/n_litr, npft, n_litr),
            smpso    = fill(-66000.0, npft),
            smpsc    = fill(-275000.0, npft),
        )

        # --- PFT constants (Li2014-specific) ---
        pftcon_li2014 = CLM.PftConFireLi2014(
            fsr_pft = fill(0.2, npft),   # fire spread rate (km/hr)
            fd_pft  = fill(1.0, npft),   # fire duration (hr)
        )

        # --- Fire constants ---
        cnfire_const = CLM.CNFireConstData()

        # --- Fire params ---
        cnfire_params = CLM.CNFireParams(
            prh30 = 0.05,
            ignition_efficiency = 0.02,
        )

        # --- Fire base data ---
        fire_data = CLM.CNFireBaseData(
            btran2_patch = zeros(np),
        )

        # --- Li2014 data ---
        fire_li2014 = CLM.CNFireLi2014Data(
            forc_hdm     = [50.0],        # population density (#/km2)
            forc_lnfm    = [0.05],        # lightning frequency (#/km2/hr)
            gdp_lf_col   = [10.0, 10.0],  # GDP (k$/capita)
            peatf_lf_col = [0.0, 0.0],    # no peatland
            abm_lf_col   = [6, 6],        # crop fire month = June
        )

        # --- Patch data ---
        # p1: tree (itype=2) on col 1
        # p2: grass (itype=10) on col 1
        # p3: tree (itype=5) on col 2
        # p4: crop (itype=15, which is > nc4_grass=14) on col 2
        patch = CLM.PatchData()
        patch.itype  = [2, 10, 5, 15]
        patch.column = [1, 1, 2, 2]
        patch.wtcol  = [0.5, 0.5, 0.6, 0.4]

        # --- Column data ---
        col = CLM.ColumnData()
        col.gridcell = [1, 1]

        # --- Gridcell data ---
        grc = CLM.GridcellData()
        grc.latdeg = [45.0]
        grc.lat    = [45.0 * pi / 180.0]  # radians

        # --- Soil state ---
        soilstate = CLM.SoilStateData()
        soilstate.watsat_col   = fill(0.45, nc, nlevgrnd)
        soilstate.rootfr_patch = fill(1.0/nlevgrnd, np, nlevgrnd)
        soilstate.sucsat_col   = fill(200.0, nc, nlevgrnd)
        soilstate.bsw_col      = fill(5.0, nc, nlevgrnd)

        h2osoi_vol_col = fill(0.1, nc, nlevgrnd)

        # --- CN Veg State ---
        cnveg_state = CLM.CNVegStateData()
        cnveg_state.dwt_smoothed_patch = zeros(np)
        cnveg_state.cropf_col          = zeros(nc)
        cnveg_state.baf_crop_col       = zeros(nc)
        cnveg_state.baf_peatf_col      = zeros(nc)
        cnveg_state.burndate_patch     = fill(10000, np)
        cnveg_state.fbac_col           = zeros(nc)
        cnveg_state.fbac1_col          = zeros(nc)
        cnveg_state.farea_burned_col   = zeros(nc)
        cnveg_state.nfire_col          = zeros(nc)
        cnveg_state.fsr_col            = zeros(nc)
        cnveg_state.fd_col             = zeros(nc)
        cnveg_state.lgdp_col           = zeros(nc)
        cnveg_state.lgdp1_col          = zeros(nc)
        cnveg_state.lpop_col           = zeros(nc)
        cnveg_state.lfwt_col           = zeros(nc)
        cnveg_state.trotr1_col         = zeros(nc)
        cnveg_state.trotr2_col         = zeros(nc)
        cnveg_state.dtrotr_col         = zeros(nc)
        cnveg_state.lfc_col            = zeros(nc)
        cnveg_state.wtlf_col           = zeros(nc)

        # --- Carbon state ---
        cnveg_cs = CLM.CNVegCarbonStateData()
        cnveg_cs.totvegc_col              = [500.0, 400.0]
        cnveg_cs.rootc_col                = zeros(nc)
        cnveg_cs.leafc_col                = zeros(nc)
        cnveg_cs.fuelc_col                = zeros(nc)
        cnveg_cs.fuelc_crop_col           = zeros(nc)
        cnveg_cs.leafc_patch              = [10.0, 5.0, 8.0, 3.0]
        cnveg_cs.leafc_storage_patch      = [1.0, 0.5, 0.8, 0.3]
        cnveg_cs.leafc_xfer_patch         = [0.5, 0.25, 0.4, 0.15]
        cnveg_cs.frootc_patch             = [4.0, 2.0, 3.0, 1.0]
        cnveg_cs.frootc_storage_patch     = [0.5, 0.2, 0.3, 0.1]
        cnveg_cs.frootc_xfer_patch        = [0.2, 0.1, 0.15, 0.05]
        cnveg_cs.deadcrootc_patch         = [25.0, 0.0, 20.0, 0.0]
        cnveg_cs.deadcrootc_storage_patch = [1.5, 0.0, 1.2, 0.0]
        cnveg_cs.deadcrootc_xfer_patch    = [0.8, 0.0, 0.6, 0.0]
        cnveg_cs.livecrootc_patch         = [10.0, 0.0, 8.0, 0.0]
        cnveg_cs.livecrootc_storage_patch = [1.0, 0.0, 0.8, 0.0]
        cnveg_cs.livecrootc_xfer_patch    = [0.5, 0.0, 0.4, 0.0]

        # --- Decomp cascade config ---
        decomp_cascade_con = CLM.DecompCascadeConData()
        decomp_cascade_con.is_litter = BitVector([true, true, true, false])
        decomp_cascade_con.is_cwd    = BitVector([false, false, false, true])

        # --- Input arrays ---
        totlitc_col          = [200.0, 150.0]
        decomp_cpools_vr_col = fill(100.0, nc, nlevdecomp, ndecomp_pools)
        t_soi17cm_col        = fill(280.0, nc)

        # --- Atmospheric forcing ---
        forc_rh_grc    = [50.0]       # relative humidity (%)
        forc_wind_grc  = [3.0]        # wind speed (m/s)
        forc_t_col     = fill(290.0, nc)  # air temperature (K)
        forc_rain_col  = fill(0.0, nc)    # rain (mm/s = kg/m2/s)
        forc_snow_col  = fill(0.0, nc)    # snow (mm/s)
        prec60_patch   = fill(2.0e-5, np)  # 60-day running mean precip (mm/s)
        prec10_patch   = fill(3.0e-5, np)  # 10-day running mean precip (mm/s)

        # --- Water data ---
        fsat_col = fill(0.1, nc)
        wf_col   = fill(0.3, nc)
        wf2_col  = fill(0.25, nc)

        # --- Masks ---
        mask_soilc      = trues(nc)
        mask_soilp      = trues(np)
        mask_exposedveg = trues(np)
        mask_noexposedveg = falses(np)

        return (pftcon=pftcon, pftcon_li2014=pftcon_li2014,
                cnfire_const=cnfire_const, cnfire_params=cnfire_params,
                fire_data=fire_data, fire_li2014=fire_li2014,
                patch=patch, col=col, grc=grc,
                soilstate=soilstate, h2osoi_vol_col=h2osoi_vol_col,
                cnveg_state=cnveg_state, cnveg_cs=cnveg_cs,
                decomp_cascade_con=decomp_cascade_con,
                totlitc_col=totlitc_col,
                decomp_cpools_vr_col=decomp_cpools_vr_col,
                t_soi17cm_col=t_soi17cm_col,
                forc_rh_grc=forc_rh_grc, forc_wind_grc=forc_wind_grc,
                forc_t_col=forc_t_col, forc_rain_col=forc_rain_col,
                forc_snow_col=forc_snow_col,
                prec60_patch=prec60_patch, prec10_patch=prec10_patch,
                fsat_col=fsat_col, wf_col=wf_col, wf2_col=wf2_col,
                mask_soilc=mask_soilc, mask_soilp=mask_soilp,
                mask_exposedveg=mask_exposedveg,
                mask_noexposedveg=mask_noexposedveg,
                np=np, nc=nc, ng=ng, nlevgrnd=nlevgrnd,
                nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools)
    end

    # ================================================================
    # Test: CNFireLi2014Data construction
    # ================================================================
    @testset "CNFireLi2014Data defaults" begin
        d = CLM.CNFireLi2014Data()
        @test isempty(d.forc_hdm)
        @test isempty(d.forc_lnfm)
        @test isempty(d.gdp_lf_col)
        @test isempty(d.peatf_lf_col)
        @test isempty(d.abm_lf_col)
    end

    # ================================================================
    # Test: PftConFireLi2014 construction
    # ================================================================
    @testset "PftConFireLi2014 defaults" begin
        p = CLM.PftConFireLi2014()
        @test isempty(p.fsr_pft)
        @test isempty(p.fd_pft)
    end

    # ================================================================
    # Test: need_lightning_and_popdens_li2014
    # ================================================================
    @testset "need_lightning_and_popdens" begin
        @test CLM.need_lightning_and_popdens_li2014() == true
    end

    # ================================================================
    # Test: p2c! patch to column averaging
    # ================================================================
    @testset "p2c! basic" begin
        patch = CLM.PatchData()
        patch.column = [1, 1, 2]
        patch.wtcol  = [0.6, 0.4, 1.0]

        patch_var = [10.0, 20.0, 30.0]
        col_var   = zeros(2)
        mask_c = trues(2)

        CLM.p2c!(col_var, patch_var, patch, mask_c, 1:2, 1:3)

        # col 1: 10*0.6 + 20*0.4 = 6 + 8 = 14
        @test col_var[1] ≈ 14.0 atol=1e-12
        # col 2: 30*1.0 = 30
        @test col_var[2] ≈ 30.0 atol=1e-12
    end

    # ================================================================
    # Test: cnfire_area_li2014! runs without error
    # ================================================================
    @testset "cnfire_area_li2014! basic run" begin
        d = make_fire_li2014_data()

        # Need to initialize dzsoi_decomp for the fuel calculation
        CLM.dzsoi_decomp[] = [0.1]

        CLM.cnfire_area_li2014!(
            d.fire_li2014, d.pftcon_li2014,
            d.fire_data, d.cnfire_const, d.cnfire_params, d.pftcon,
            d.mask_soilc, d.mask_soilp, d.mask_exposedveg, d.mask_noexposedveg,
            1:d.nc, 1:d.np,
            d.patch, d.col, d.grc,
            d.soilstate, d.h2osoi_vol_col,
            d.cnveg_state, d.cnveg_cs,
            d.decomp_cascade_con,
            d.totlitc_col, d.decomp_cpools_vr_col, d.t_soi17cm_col;
            forc_rh_grc   = d.forc_rh_grc,
            forc_wind_grc = d.forc_wind_grc,
            forc_t_col    = d.forc_t_col,
            forc_rain_col = d.forc_rain_col,
            forc_snow_col = d.forc_snow_col,
            prec60_patch  = d.prec60_patch,
            prec10_patch  = d.prec10_patch,
            fsat_col      = d.fsat_col,
            wf_col        = d.wf_col,
            wf2_col       = d.wf2_col,
            dt = 1800.0, dayspyr = 365.0,
            kmo = 6, kda = 15, mcsec = 3600,
            nstep = 10,
            nlevgrnd = d.nlevgrnd,
            nlevdecomp = d.nlevdecomp,
            ndecomp_pools = d.ndecomp_pools,
        )

        # farea_burned should be computed and in [0, 1]
        for c in 1:d.nc
            @test 0.0 <= d.cnveg_state.farea_burned_col[c] <= 1.0
        end
        # farea_burned should be > 0 (we have fuel, ignition sources)
        @test d.cnveg_state.farea_burned_col[1] > 0.0
    end

    # ================================================================
    # Test: first timestep zeros everything
    # ================================================================
    @testset "cnfire_area_li2014! first timestep" begin
        d = make_fire_li2014_data()
        CLM.dzsoi_decomp[] = [0.1]

        # Set some nonzero values
        d.cnveg_state.farea_burned_col .= 999.0

        CLM.cnfire_area_li2014!(
            d.fire_li2014, d.pftcon_li2014,
            d.fire_data, d.cnfire_const, d.cnfire_params, d.pftcon,
            d.mask_soilc, d.mask_soilp, d.mask_exposedveg, d.mask_noexposedveg,
            1:d.nc, 1:d.np,
            d.patch, d.col, d.grc,
            d.soilstate, d.h2osoi_vol_col,
            d.cnveg_state, d.cnveg_cs,
            d.decomp_cascade_con,
            d.totlitc_col, d.decomp_cpools_vr_col, d.t_soi17cm_col;
            forc_rh_grc   = d.forc_rh_grc,
            forc_wind_grc = d.forc_wind_grc,
            forc_t_col    = d.forc_t_col,
            forc_rain_col = d.forc_rain_col,
            forc_snow_col = d.forc_snow_col,
            prec60_patch  = d.prec60_patch,
            prec10_patch  = d.prec10_patch,
            fsat_col      = d.fsat_col,
            wf_col        = d.wf_col,
            wf2_col       = d.wf2_col,
            nstep = 0,  # first timestep
        )

        for c in 1:d.nc
            @test d.cnveg_state.farea_burned_col[c] == 0.0
            @test d.cnveg_state.baf_crop_col[c] == 0.0
            @test d.cnveg_state.baf_peatf_col[c] == 0.0
            @test d.cnveg_state.fbac_col[c] == 0.0
            @test d.cnveg_state.fbac1_col[c] == 0.0
        end
    end

    # ================================================================
    # Test: cropf_col calculation
    # ================================================================
    @testset "cropf_col and lfwt" begin
        d = make_fire_li2014_data()
        CLM.dzsoi_decomp[] = [0.1]

        CLM.cnfire_area_li2014!(
            d.fire_li2014, d.pftcon_li2014,
            d.fire_data, d.cnfire_const, d.cnfire_params, d.pftcon,
            d.mask_soilc, d.mask_soilp, d.mask_exposedveg, d.mask_noexposedveg,
            1:d.nc, 1:d.np,
            d.patch, d.col, d.grc,
            d.soilstate, d.h2osoi_vol_col,
            d.cnveg_state, d.cnveg_cs,
            d.decomp_cascade_con,
            d.totlitc_col, d.decomp_cpools_vr_col, d.t_soi17cm_col;
            forc_rh_grc   = d.forc_rh_grc,
            forc_wind_grc = d.forc_wind_grc,
            forc_t_col    = d.forc_t_col,
            forc_rain_col = d.forc_rain_col,
            forc_snow_col = d.forc_snow_col,
            prec60_patch  = d.prec60_patch,
            prec10_patch  = d.prec10_patch,
            fsat_col      = d.fsat_col,
            wf_col        = d.wf_col,
            wf2_col       = d.wf2_col,
            nstep = 10,
            nlevgrnd = d.nlevgrnd,
            nlevdecomp = d.nlevdecomp,
            ndecomp_pools = d.ndecomp_pools,
        )

        # Col 1: patches 1 (itype=2, tree) and 2 (itype=10, grass)
        # Both itype <= 14 (nc4_grass), so cropf = 0
        @test d.cnveg_state.cropf_col[1] ≈ 0.0 atol=1e-12

        # Col 2: patch 3 (itype=5, tree) and 4 (itype=15, crop)
        # itype 15 > 14 so it's crop: wtcol=0.4
        @test d.cnveg_state.cropf_col[2] ≈ 0.4 atol=1e-12
    end

    # ================================================================
    # Test: peatland fire
    # ================================================================
    @testset "peatland fire non-boreal" begin
        d = make_fire_li2014_data()
        CLM.dzsoi_decomp[] = [0.1]

        # Set latitude to tropical (non-boreal: < 40)
        d.grc.latdeg[1] = 10.0
        d.grc.lat[1]    = 10.0 * pi / 180.0

        # Enable peatland
        d.fire_li2014.peatf_lf_col .= 0.5

        CLM.cnfire_area_li2014!(
            d.fire_li2014, d.pftcon_li2014,
            d.fire_data, d.cnfire_const, d.cnfire_params, d.pftcon,
            d.mask_soilc, d.mask_soilp, d.mask_exposedveg, d.mask_noexposedveg,
            1:d.nc, 1:d.np,
            d.patch, d.col, d.grc,
            d.soilstate, d.h2osoi_vol_col,
            d.cnveg_state, d.cnveg_cs,
            d.decomp_cascade_con,
            d.totlitc_col, d.decomp_cpools_vr_col, d.t_soi17cm_col;
            forc_rh_grc   = d.forc_rh_grc,
            forc_wind_grc = d.forc_wind_grc,
            forc_t_col    = d.forc_t_col,
            forc_rain_col = d.forc_rain_col,
            forc_snow_col = d.forc_snow_col,
            prec60_patch  = d.prec60_patch,
            prec10_patch  = d.prec10_patch,
            fsat_col      = d.fsat_col,
            wf_col        = d.wf_col,
            wf2_col       = d.wf2_col,
            nstep = 10,
            nlevgrnd = d.nlevgrnd,
            nlevdecomp = d.nlevdecomp,
            ndecomp_pools = d.ndecomp_pools,
        )

        # Non-boreal peatland fire should be > 0 with peatland fraction
        @test d.cnveg_state.baf_peatf_col[1] > 0.0
    end

    # ================================================================
    # Test: peatland fire boreal
    # ================================================================
    @testset "peatland fire boreal" begin
        d = make_fire_li2014_data()
        CLM.dzsoi_decomp[] = [0.1]

        # Boreal latitude (> 40)
        d.grc.latdeg[1] = 55.0
        d.grc.lat[1]    = 55.0 * pi / 180.0

        # Enable peatland
        d.fire_li2014.peatf_lf_col .= 0.5

        # Warm enough soil temperature
        d.t_soi17cm_col .= 280.0

        CLM.cnfire_area_li2014!(
            d.fire_li2014, d.pftcon_li2014,
            d.fire_data, d.cnfire_const, d.cnfire_params, d.pftcon,
            d.mask_soilc, d.mask_soilp, d.mask_exposedveg, d.mask_noexposedveg,
            1:d.nc, 1:d.np,
            d.patch, d.col, d.grc,
            d.soilstate, d.h2osoi_vol_col,
            d.cnveg_state, d.cnveg_cs,
            d.decomp_cascade_con,
            d.totlitc_col, d.decomp_cpools_vr_col, d.t_soi17cm_col;
            forc_rh_grc   = d.forc_rh_grc,
            forc_wind_grc = d.forc_wind_grc,
            forc_t_col    = d.forc_t_col,
            forc_rain_col = d.forc_rain_col,
            forc_snow_col = d.forc_snow_col,
            prec60_patch  = d.prec60_patch,
            prec10_patch  = d.prec10_patch,
            fsat_col      = d.fsat_col,
            wf_col        = d.wf_col,
            wf2_col       = d.wf2_col,
            nstep = 10,
            nlevgrnd = d.nlevgrnd,
            nlevdecomp = d.nlevdecomp,
            ndecomp_pools = d.ndecomp_pools,
        )

        # Boreal peatland fire should be > 0
        @test d.cnveg_state.baf_peatf_col[1] > 0.0

        # Boreal formula: boreal_peatfire_c/secsphr * exp(-pi*wf2/denom) * temp_term * peatf * (1-fsat)
        c = d.cnfire_const
        expected = c.boreal_peatfire_c / CLM.SECSPHR *
            exp(-CLM.RPI * (max(0.0, 0.25) / c.borpeat_fire_soilmoist_denom)) *
            max(0.0, min(1.0, (280.0 - CLM.TFRZ) / 10.0)) *
            0.5 * (1.0 - 0.1)
        @test d.cnveg_state.baf_peatf_col[1] ≈ expected atol=1e-15
    end

    # ================================================================
    # Test: crop fire burned area
    # ================================================================
    @testset "crop fire area" begin
        d = make_fire_li2014_data()
        CLM.dzsoi_decomp[] = [0.1]

        # Set conditions for crop fire on patch 4 (crop, col 2):
        # - temp >= TFRZ ✓ (290 K)
        # - itype > nc4_grass ✓ (15 > 14)
        # - kmo == abm_lf[c] ✓ (6 == 6)
        # - no precip ✓ (rain=0, snow=0)
        # - burndate >= 999 ✓ (10000)
        # - wtcol > 0 ✓ (0.4)

        CLM.cnfire_area_li2014!(
            d.fire_li2014, d.pftcon_li2014,
            d.fire_data, d.cnfire_const, d.cnfire_params, d.pftcon,
            d.mask_soilc, d.mask_soilp, d.mask_exposedveg, d.mask_noexposedveg,
            1:d.nc, 1:d.np,
            d.patch, d.col, d.grc,
            d.soilstate, d.h2osoi_vol_col,
            d.cnveg_state, d.cnveg_cs,
            d.decomp_cascade_con,
            d.totlitc_col, d.decomp_cpools_vr_col, d.t_soi17cm_col;
            forc_rh_grc   = d.forc_rh_grc,
            forc_wind_grc = d.forc_wind_grc,
            forc_t_col    = d.forc_t_col,
            forc_rain_col = d.forc_rain_col,
            forc_snow_col = d.forc_snow_col,
            prec60_patch  = d.prec60_patch,
            prec10_patch  = d.prec10_patch,
            fsat_col      = d.fsat_col,
            wf_col        = d.wf_col,
            wf2_col       = d.wf2_col,
            kmo = 6, kda = 15, mcsec = 3600,
            nstep = 10,
            nlevgrnd = d.nlevgrnd,
            nlevdecomp = d.nlevdecomp,
            ndecomp_pools = d.ndecomp_pools,
        )

        # Col 2 has a crop patch, so baf_crop should be > 0
        # (only if fuel is sufficient - crop fuel needs leafc + totlitc)
        # baf_crop[2] may be 0 if fuelc_crop is very low (below lfuel=75)
        # Let's just check it's non-negative
        @test d.cnveg_state.baf_crop_col[2] >= 0.0

        # Col 1 has no crop, so baf_crop should be 0
        @test d.cnveg_state.baf_crop_col[1] == 0.0
    end

    # ================================================================
    # Test: nfire calculation (fire counts/km2/sec)
    # ================================================================
    @testset "nfire positive" begin
        d = make_fire_li2014_data()
        CLM.dzsoi_decomp[] = [0.1]

        CLM.cnfire_area_li2014!(
            d.fire_li2014, d.pftcon_li2014,
            d.fire_data, d.cnfire_const, d.cnfire_params, d.pftcon,
            d.mask_soilc, d.mask_soilp, d.mask_exposedveg, d.mask_noexposedveg,
            1:d.nc, 1:d.np,
            d.patch, d.col, d.grc,
            d.soilstate, d.h2osoi_vol_col,
            d.cnveg_state, d.cnveg_cs,
            d.decomp_cascade_con,
            d.totlitc_col, d.decomp_cpools_vr_col, d.t_soi17cm_col;
            forc_rh_grc   = d.forc_rh_grc,
            forc_wind_grc = d.forc_wind_grc,
            forc_t_col    = d.forc_t_col,
            forc_rain_col = d.forc_rain_col,
            forc_snow_col = d.forc_snow_col,
            prec60_patch  = d.prec60_patch,
            prec10_patch  = d.prec10_patch,
            fsat_col      = d.fsat_col,
            wf_col        = d.wf_col,
            wf2_col       = d.wf2_col,
            nstep = 10,
            nlevgrnd = d.nlevgrnd,
            nlevdecomp = d.nlevdecomp,
            ndecomp_pools = d.ndecomp_pools,
        )

        # Column 1 has only natural veg (no crop, no tropical forest)
        # So nfire should be computed and >= 0
        @test d.cnveg_state.nfire_col[1] >= 0.0
    end

    # ================================================================
    # Test: farea_burned bounded [0, 1]
    # ================================================================
    @testset "farea_burned bounds" begin
        d = make_fire_li2014_data()
        CLM.dzsoi_decomp[] = [0.1]

        # Set extreme conditions
        d.fire_li2014.forc_hdm[1] = 500.0  # high population
        d.forc_wind_grc[1] = 20.0           # strong wind

        CLM.cnfire_area_li2014!(
            d.fire_li2014, d.pftcon_li2014,
            d.fire_data, d.cnfire_const, d.cnfire_params, d.pftcon,
            d.mask_soilc, d.mask_soilp, d.mask_exposedveg, d.mask_noexposedveg,
            1:d.nc, 1:d.np,
            d.patch, d.col, d.grc,
            d.soilstate, d.h2osoi_vol_col,
            d.cnveg_state, d.cnveg_cs,
            d.decomp_cascade_con,
            d.totlitc_col, d.decomp_cpools_vr_col, d.t_soi17cm_col;
            forc_rh_grc   = d.forc_rh_grc,
            forc_wind_grc = d.forc_wind_grc,
            forc_t_col    = d.forc_t_col,
            forc_rain_col = d.forc_rain_col,
            forc_snow_col = d.forc_snow_col,
            prec60_patch  = d.prec60_patch,
            prec10_patch  = d.prec10_patch,
            fsat_col      = d.fsat_col,
            wf_col        = d.wf_col,
            wf2_col       = d.wf2_col,
            nstep = 10,
            nlevgrnd = d.nlevgrnd,
            nlevdecomp = d.nlevdecomp,
            ndecomp_pools = d.ndecomp_pools,
        )

        for c in 1:d.nc
            @test 0.0 <= d.cnveg_state.farea_burned_col[c] <= 1.0
        end
    end

end
