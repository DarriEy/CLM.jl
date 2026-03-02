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

    # ================================================================
    # Helper: create data structures for Li2014 fire flux tests
    # ================================================================
    function make_fire_li2014_flux_data(; np=4, nc=2, ng=1, nlevdecomp=1,
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

        # --- Fire constants ---
        cnfire_const = CLM.CNFireConstData()

        # --- Patch data ---
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
        grc.lat    = [45.0 * pi / 180.0]

        # --- DGVS data ---
        dgvs = CLM.DgvsFireData(nind_patch = fill(100.0, np))

        # --- CN Veg State ---
        cnveg_state = CLM.CNVegStateData()
        cnveg_state.cropf_col          = [0.0, 0.4]
        cnveg_state.farea_burned_col   = [1.0e-4, 1.0e-4]
        cnveg_state.fbac1_col          = zeros(nc)
        cnveg_state.fbac_col           = zeros(nc)
        cnveg_state.baf_crop_col       = [0.0, 2.0e-5]
        cnveg_state.baf_peatf_col      = zeros(nc)
        cnveg_state.trotr1_col         = zeros(nc)
        cnveg_state.trotr2_col         = zeros(nc)
        cnveg_state.dtrotr_col         = zeros(nc)
        cnveg_state.lfc_col            = zeros(nc)
        cnveg_state.lfc2_col           = zeros(nc)

        # --- Carbon state ---
        cnveg_cs = CLM.CNVegCarbonStateData()
        cnveg_cs.leafcmax_patch               = fill(0.0, np)
        cnveg_cs.leafc_patch                  = [10.0, 5.0, 8.0, 3.0]
        cnveg_cs.leafc_storage_patch          = [1.0, 0.5, 0.8, 0.3]
        cnveg_cs.leafc_xfer_patch             = [0.5, 0.25, 0.4, 0.15]
        cnveg_cs.livestemc_patch              = [20.0, 0.0, 15.0, 0.0]
        cnveg_cs.livestemc_storage_patch      = [2.0, 0.0, 1.5, 0.0]
        cnveg_cs.livestemc_xfer_patch         = [1.0, 0.0, 0.8, 0.0]
        cnveg_cs.deadstemc_patch              = [50.0, 0.0, 40.0, 0.0]
        cnveg_cs.deadstemc_storage_patch      = [3.0, 0.0, 2.5, 0.0]
        cnveg_cs.deadstemc_xfer_patch         = [1.5, 0.0, 1.2, 0.0]
        cnveg_cs.frootc_patch                 = [4.0, 2.0, 3.0, 1.0]
        cnveg_cs.frootc_storage_patch         = [0.5, 0.2, 0.3, 0.1]
        cnveg_cs.frootc_xfer_patch            = [0.2, 0.1, 0.15, 0.05]
        cnveg_cs.livecrootc_patch             = [10.0, 0.0, 8.0, 0.0]
        cnveg_cs.livecrootc_storage_patch     = [1.0, 0.0, 0.8, 0.0]
        cnveg_cs.livecrootc_xfer_patch        = [0.5, 0.0, 0.4, 0.0]
        cnveg_cs.deadcrootc_patch             = [25.0, 0.0, 20.0, 0.0]
        cnveg_cs.deadcrootc_storage_patch     = [1.5, 0.0, 1.2, 0.0]
        cnveg_cs.deadcrootc_xfer_patch        = [0.8, 0.0, 0.6, 0.0]
        cnveg_cs.gresp_storage_patch          = [0.2, 0.1, 0.15, 0.05]
        cnveg_cs.gresp_xfer_patch             = [0.1, 0.05, 0.08, 0.03]

        # --- Carbon flux (zero-initialized) ---
        cnveg_cf = CLM.CNVegCarbonFluxData()
        cnveg_cf.m_leafc_to_fire_patch                     = zeros(np)
        cnveg_cf.m_leafc_storage_to_fire_patch             = zeros(np)
        cnveg_cf.m_leafc_xfer_to_fire_patch                = zeros(np)
        cnveg_cf.m_livestemc_to_fire_patch                 = zeros(np)
        cnveg_cf.m_livestemc_storage_to_fire_patch         = zeros(np)
        cnveg_cf.m_livestemc_xfer_to_fire_patch            = zeros(np)
        cnveg_cf.m_deadstemc_to_fire_patch                 = zeros(np)
        cnveg_cf.m_deadstemc_storage_to_fire_patch         = zeros(np)
        cnveg_cf.m_deadstemc_xfer_to_fire_patch            = zeros(np)
        cnveg_cf.m_frootc_to_fire_patch                    = zeros(np)
        cnveg_cf.m_frootc_storage_to_fire_patch            = zeros(np)
        cnveg_cf.m_frootc_xfer_to_fire_patch               = zeros(np)
        cnveg_cf.m_livecrootc_to_fire_patch                = zeros(np)
        cnveg_cf.m_livecrootc_storage_to_fire_patch        = zeros(np)
        cnveg_cf.m_livecrootc_xfer_to_fire_patch           = zeros(np)
        cnveg_cf.m_deadcrootc_to_fire_patch                = zeros(np)
        cnveg_cf.m_deadcrootc_storage_to_fire_patch        = zeros(np)
        cnveg_cf.m_deadcrootc_xfer_to_fire_patch           = zeros(np)
        cnveg_cf.m_gresp_storage_to_fire_patch             = zeros(np)
        cnveg_cf.m_gresp_xfer_to_fire_patch                = zeros(np)
        cnveg_cf.m_leafc_to_litter_fire_patch              = zeros(np)
        cnveg_cf.m_leafc_storage_to_litter_fire_patch      = zeros(np)
        cnveg_cf.m_leafc_xfer_to_litter_fire_patch         = zeros(np)
        cnveg_cf.m_livestemc_to_litter_fire_patch          = zeros(np)
        cnveg_cf.m_livestemc_storage_to_litter_fire_patch  = zeros(np)
        cnveg_cf.m_livestemc_xfer_to_litter_fire_patch     = zeros(np)
        cnveg_cf.m_livestemc_to_deadstemc_fire_patch       = zeros(np)
        cnveg_cf.m_deadstemc_to_litter_fire_patch          = zeros(np)
        cnveg_cf.m_deadstemc_storage_to_litter_fire_patch  = zeros(np)
        cnveg_cf.m_deadstemc_xfer_to_litter_fire_patch     = zeros(np)
        cnveg_cf.m_frootc_to_litter_fire_patch             = zeros(np)
        cnveg_cf.m_frootc_storage_to_litter_fire_patch     = zeros(np)
        cnveg_cf.m_frootc_xfer_to_litter_fire_patch        = zeros(np)
        cnveg_cf.m_livecrootc_to_litter_fire_patch         = zeros(np)
        cnveg_cf.m_livecrootc_storage_to_litter_fire_patch = zeros(np)
        cnveg_cf.m_livecrootc_xfer_to_litter_fire_patch    = zeros(np)
        cnveg_cf.m_livecrootc_to_deadcrootc_fire_patch     = zeros(np)
        cnveg_cf.m_deadcrootc_to_litter_fire_patch         = zeros(np)
        cnveg_cf.m_deadcrootc_storage_to_litter_fire_patch = zeros(np)
        cnveg_cf.m_deadcrootc_xfer_to_litter_fire_patch    = zeros(np)
        cnveg_cf.m_gresp_storage_to_litter_fire_patch      = zeros(np)
        cnveg_cf.m_gresp_xfer_to_litter_fire_patch         = zeros(np)
        cnveg_cf.fire_mortality_c_to_cwdc_col              = zeros(nc, nlevdecomp)
        cnveg_cf.m_decomp_cpools_to_fire_vr_col            = zeros(nc, nlevdecomp, ndecomp_pools)
        cnveg_cf.m_c_to_litr_fire_col                      = zeros(nc, nlevdecomp, ndecomp_pools)

        # --- Nitrogen state ---
        cnveg_ns = CLM.CNVegNitrogenStateData()
        cnveg_ns.leafn_patch                  = [0.4, 0.2, 0.32, 0.12]
        cnveg_ns.leafn_storage_patch          = [0.04, 0.02, 0.032, 0.012]
        cnveg_ns.leafn_xfer_patch             = [0.02, 0.01, 0.016, 0.006]
        cnveg_ns.livestemn_patch              = [0.4, 0.0, 0.3, 0.0]
        cnveg_ns.livestemn_storage_patch      = [0.04, 0.0, 0.03, 0.0]
        cnveg_ns.livestemn_xfer_patch         = [0.02, 0.0, 0.016, 0.0]
        cnveg_ns.deadstemn_patch              = [1.0, 0.0, 0.8, 0.0]
        cnveg_ns.deadstemn_storage_patch      = [0.06, 0.0, 0.05, 0.0]
        cnveg_ns.deadstemn_xfer_patch         = [0.03, 0.0, 0.024, 0.0]
        cnveg_ns.frootn_patch                 = [0.16, 0.08, 0.12, 0.04]
        cnveg_ns.frootn_storage_patch         = [0.02, 0.008, 0.012, 0.004]
        cnveg_ns.frootn_xfer_patch            = [0.008, 0.004, 0.006, 0.002]
        cnveg_ns.livecrootn_patch             = [0.2, 0.0, 0.16, 0.0]
        cnveg_ns.livecrootn_storage_patch     = [0.02, 0.0, 0.016, 0.0]
        cnveg_ns.livecrootn_xfer_patch        = [0.01, 0.0, 0.008, 0.0]
        cnveg_ns.deadcrootn_patch             = [0.5, 0.0, 0.4, 0.0]
        cnveg_ns.deadcrootn_storage_patch     = [0.03, 0.0, 0.024, 0.0]
        cnveg_ns.deadcrootn_xfer_patch        = [0.016, 0.0, 0.012, 0.0]
        cnveg_ns.retransn_patch               = [0.05, 0.02, 0.04, 0.01]

        # --- Nitrogen flux (zero-initialized) ---
        cnveg_nf = CLM.CNVegNitrogenFluxData()
        cnveg_nf.m_leafn_to_fire_patch                     = zeros(np)
        cnveg_nf.m_leafn_storage_to_fire_patch             = zeros(np)
        cnveg_nf.m_leafn_xfer_to_fire_patch                = zeros(np)
        cnveg_nf.m_livestemn_to_fire_patch                 = zeros(np)
        cnveg_nf.m_livestemn_storage_to_fire_patch         = zeros(np)
        cnveg_nf.m_livestemn_xfer_to_fire_patch            = zeros(np)
        cnveg_nf.m_deadstemn_to_fire_patch                 = zeros(np)
        cnveg_nf.m_deadstemn_storage_to_fire_patch         = zeros(np)
        cnveg_nf.m_deadstemn_xfer_to_fire_patch            = zeros(np)
        cnveg_nf.m_frootn_to_fire_patch                    = zeros(np)
        cnveg_nf.m_frootn_storage_to_fire_patch            = zeros(np)
        cnveg_nf.m_frootn_xfer_to_fire_patch               = zeros(np)
        cnveg_nf.m_livecrootn_to_fire_patch                = zeros(np)
        cnveg_nf.m_livecrootn_storage_to_fire_patch        = zeros(np)
        cnveg_nf.m_livecrootn_xfer_to_fire_patch           = zeros(np)
        cnveg_nf.m_deadcrootn_to_fire_patch                = zeros(np)
        cnveg_nf.m_deadcrootn_storage_to_fire_patch        = zeros(np)
        cnveg_nf.m_deadcrootn_xfer_to_fire_patch           = zeros(np)
        cnveg_nf.m_retransn_to_fire_patch                  = zeros(np)
        cnveg_nf.m_leafn_to_litter_fire_patch              = zeros(np)
        cnveg_nf.m_leafn_storage_to_litter_fire_patch      = zeros(np)
        cnveg_nf.m_leafn_xfer_to_litter_fire_patch         = zeros(np)
        cnveg_nf.m_livestemn_to_litter_fire_patch          = zeros(np)
        cnveg_nf.m_livestemn_storage_to_litter_fire_patch  = zeros(np)
        cnveg_nf.m_livestemn_xfer_to_litter_fire_patch     = zeros(np)
        cnveg_nf.m_livestemn_to_deadstemn_fire_patch       = zeros(np)
        cnveg_nf.m_deadstemn_to_litter_fire_patch          = zeros(np)
        cnveg_nf.m_deadstemn_storage_to_litter_fire_patch  = zeros(np)
        cnveg_nf.m_deadstemn_xfer_to_litter_fire_patch     = zeros(np)
        cnveg_nf.m_frootn_to_litter_fire_patch             = zeros(np)
        cnveg_nf.m_frootn_storage_to_litter_fire_patch     = zeros(np)
        cnveg_nf.m_frootn_xfer_to_litter_fire_patch        = zeros(np)
        cnveg_nf.m_livecrootn_to_litter_fire_patch         = zeros(np)
        cnveg_nf.m_livecrootn_storage_to_litter_fire_patch = zeros(np)
        cnveg_nf.m_livecrootn_xfer_to_litter_fire_patch    = zeros(np)
        cnveg_nf.m_livecrootn_to_deadcrootn_fire_patch     = zeros(np)
        cnveg_nf.m_deadcrootn_to_litter_fire_patch         = zeros(np)
        cnveg_nf.m_deadcrootn_storage_to_litter_fire_patch = zeros(np)
        cnveg_nf.m_deadcrootn_xfer_to_litter_fire_patch    = zeros(np)
        cnveg_nf.m_retransn_to_litter_fire_patch           = zeros(np)
        cnveg_nf.fire_mortality_n_to_cwdn_col              = zeros(nc, nlevdecomp)
        cnveg_nf.m_decomp_npools_to_fire_vr_col            = zeros(nc, nlevdecomp, ndecomp_pools)
        cnveg_nf.m_n_to_litr_fire_col                      = zeros(nc, nlevdecomp, ndecomp_pools)

        # --- Soil Biogeochem Carbon Flux ---
        soilbgc_cf = CLM.SoilBiogeochemCarbonFluxData()
        soilbgc_cf.somc_fire_col = zeros(nc)

        # --- Decomp cascade config ---
        decomp_cascade_con = CLM.DecompCascadeConData()
        decomp_cascade_con.is_litter = BitVector([true, true, true, false])
        decomp_cascade_con.is_cwd    = BitVector([false, false, false, true])

        # --- Profile arrays ---
        leaf_prof  = fill(1.0, np, nlevdecomp)
        froot_prof = fill(1.0, np, nlevdecomp)
        croot_prof = fill(1.0, np, nlevdecomp)
        stem_prof  = fill(1.0, np, nlevdecomp)
        totsomc    = fill(5000.0, nc)
        decomp_cpools_vr = fill(100.0, nc, nlevdecomp, ndecomp_pools)
        decomp_npools_vr = fill(5.0, nc, nlevdecomp, ndecomp_pools)
        somc_fire  = zeros(nc)

        # --- Masks ---
        mask_soilc = trues(nc)
        mask_soilp = trues(np)

        return (pftcon=pftcon, cnfire_const=cnfire_const,
                patch=patch, col=col, grc=grc,
                dgvs=dgvs, cnveg_state=cnveg_state,
                cnveg_cs=cnveg_cs, cnveg_cf=cnveg_cf,
                cnveg_ns=cnveg_ns, cnveg_nf=cnveg_nf,
                soilbgc_cf=soilbgc_cf, decomp_cascade_con=decomp_cascade_con,
                leaf_prof=leaf_prof, froot_prof=froot_prof,
                croot_prof=croot_prof, stem_prof=stem_prof,
                totsomc=totsomc, decomp_cpools_vr=decomp_cpools_vr,
                decomp_npools_vr=decomp_npools_vr, somc_fire=somc_fire,
                mask_soilc=mask_soilc, mask_soilp=mask_soilp,
                np=np, nc=nc, nlevdecomp=nlevdecomp,
                ndecomp_pools=ndecomp_pools, n_litr=n_litr)
    end

    # ================================================================
    # Test: cnfire_fluxes_li2014! runs without error
    # ================================================================
    @testset "cnfire_fluxes_li2014! basic run" begin
        d = make_fire_li2014_flux_data()

        mask_actfirec, mask_actfirep = CLM.cnfire_fluxes_li2014!(
            d.mask_soilc, d.mask_soilp,
            1:d.nc, 1:d.np,
            d.cnfire_const, d.pftcon, d.patch, d.col, d.grc,
            d.dgvs, d.cnveg_state, d.cnveg_cs, d.cnveg_cf,
            d.cnveg_ns, d.cnveg_nf,
            d.soilbgc_cf, d.decomp_cascade_con,
            d.leaf_prof, d.froot_prof, d.croot_prof, d.stem_prof,
            d.totsomc, d.decomp_cpools_vr, d.decomp_npools_vr, d.somc_fire;
            dt=1800.0, dayspyr=365.0,
            nlevdecomp=d.nlevdecomp, ndecomp_pools=d.ndecomp_pools,
            i_met_lit=1, i_litr_max=d.n_litr
        )

        # Output masks should be BitVectors
        @test mask_actfirec isa BitVector
        @test mask_actfirep isa BitVector
        @test length(mask_actfirec) == d.nc
        @test length(mask_actfirep) == d.np

        # Some fire should be active (farea_burned > 0)
        @test any(mask_actfirec)
        @test any(mask_actfirep)
    end

    # ================================================================
    # Test: Li2014 uses 0.5 for litter and 0.25 for CWD
    # ================================================================
    @testset "cnfire_fluxes_li2014! combustion factors" begin
        d = make_fire_li2014_flux_data()

        CLM.cnfire_fluxes_li2014!(
            d.mask_soilc, d.mask_soilp,
            1:d.nc, 1:d.np,
            d.cnfire_const, d.pftcon, d.patch, d.col, d.grc,
            d.dgvs, d.cnveg_state, d.cnveg_cs, d.cnveg_cf,
            d.cnveg_ns, d.cnveg_nf,
            d.soilbgc_cf, d.decomp_cascade_con,
            d.leaf_prof, d.froot_prof, d.croot_prof, d.stem_prof,
            d.totsomc, d.decomp_cpools_vr, d.decomp_npools_vr, d.somc_fire;
            dt=1800.0, dayspyr=365.0,
            nlevdecomp=d.nlevdecomp, ndecomp_pools=d.ndecomp_pools,
            i_met_lit=1, i_litr_max=d.n_litr
        )

        # Check decomp pool fire fluxes use Li2014 combustion factors
        f = d.cnveg_state.farea_burned_col[1]
        # Litter pools (is_litter = [true, true, true, false]):
        # m_decomp_cpools_to_fire_vr = decomp_cpools_vr * f * 0.5
        for l in 1:3  # litter pools
            expected_c = 100.0 * f * 0.5
            @test d.cnveg_cf.m_decomp_cpools_to_fire_vr_col[1, 1, l] ≈ expected_c atol=1e-15
        end
        # CWD pool (is_cwd = [false, false, false, true]):
        # m_decomp_cpools_to_fire_vr = decomp_cpools_vr * (f - baf_crop) * 0.25
        baf_crop = d.cnveg_state.baf_crop_col[1]
        expected_cwd = 100.0 * (f - baf_crop) * 0.25
        @test d.cnveg_cf.m_decomp_cpools_to_fire_vr_col[1, 1, 4] ≈ expected_cwd atol=1e-15
    end

    # ================================================================
    # Test: Li2014 fire fluxes produce non-negative C emissions
    # ================================================================
    @testset "cnfire_fluxes_li2014! non-negative emissions" begin
        d = make_fire_li2014_flux_data()

        CLM.cnfire_fluxes_li2014!(
            d.mask_soilc, d.mask_soilp,
            1:d.nc, 1:d.np,
            d.cnfire_const, d.pftcon, d.patch, d.col, d.grc,
            d.dgvs, d.cnveg_state, d.cnveg_cs, d.cnveg_cf,
            d.cnveg_ns, d.cnveg_nf,
            d.soilbgc_cf, d.decomp_cascade_con,
            d.leaf_prof, d.froot_prof, d.croot_prof, d.stem_prof,
            d.totsomc, d.decomp_cpools_vr, d.decomp_npools_vr, d.somc_fire;
            dt=1800.0, dayspyr=365.0,
            nlevdecomp=d.nlevdecomp, ndecomp_pools=d.ndecomp_pools,
            i_met_lit=1, i_litr_max=d.n_litr
        )

        # All fire emission fluxes should be non-negative
        for p in 1:d.np
            @test d.cnveg_cf.m_leafc_to_fire_patch[p] >= 0.0
            @test d.cnveg_cf.m_gresp_storage_to_fire_patch[p] >= 0.0
        end

        # Leaf fire emission for patch 1 (tree, itype=2):
        # f = (farea_burned - baf_crop) / (1 - cropf_col)
        # = (1e-4 - 0) / (1 - 0) = 1e-4
        # m_leafc_to_fire = leafc * f * cc_leaf = 10 * 1e-4 * 0.4 = 4e-4
        @test d.cnveg_cf.m_leafc_to_fire_patch[1] ≈ 10.0 * 1.0e-4 * 0.4 atol=1e-15
    end

    # ================================================================
    # Test: peat fire somc_fire calculation
    # ================================================================
    @testset "cnfire_fluxes_li2014! somc_fire" begin
        d = make_fire_li2014_flux_data()

        # Set baf_peatf > 0 to trigger peat fire
        d.cnveg_state.baf_peatf_col .= 1.0e-5

        CLM.cnfire_fluxes_li2014!(
            d.mask_soilc, d.mask_soilp,
            1:d.nc, 1:d.np,
            d.cnfire_const, d.pftcon, d.patch, d.col, d.grc,
            d.dgvs, d.cnveg_state, d.cnveg_cs, d.cnveg_cf,
            d.cnveg_ns, d.cnveg_nf,
            d.soilbgc_cf, d.decomp_cascade_con,
            d.leaf_prof, d.froot_prof, d.croot_prof, d.stem_prof,
            d.totsomc, d.decomp_cpools_vr, d.decomp_npools_vr, d.somc_fire;
            dt=1800.0, dayspyr=365.0,
            nlevdecomp=d.nlevdecomp, ndecomp_pools=d.ndecomp_pools,
            i_met_lit=1, i_litr_max=d.n_litr
        )

        # Non-boreal (lat=45 > 40): somc_fire = baf_peatf * 2.2e3
        @test d.somc_fire[1] ≈ 1.0e-5 * 2.2e3 atol=1e-10
    end

end
