@testset "Fire methods (Li2016/2021/2024/NoFire + factory)" begin

    # ----------------------------------------------------------------
    # Helper: minimal AREA fixture (mirrors test_fire_li2014's make_fire_li2014_data)
    # plus the extra inputs the non-default methods need (rswf, rh30, prec30, croplive).
    # ----------------------------------------------------------------
    function make_fire_methods_area_data(; np=4, nc=2, ng=1, nlevgrnd=3,
                                           nlevdecomp=1, ndecomp_pools=4, n_litr=3)
        npft = 20

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
            rswf_min = fill(0.1, npft),   # Li2021/2024
            rswf_max = fill(0.9, npft),   # Li2021/2024
        )

        pftcon_li2014 = CLM.PftConFireLi2014(
            fsr_pft = fill(0.2, npft),
            fd_pft  = fill(1.0, npft),
        )

        cnfire_const  = CLM.CNFireConstData()
        cnfire_params = CLM.CNFireParams(prh30 = 0.05, ignition_efficiency = 0.02)
        fire_data     = CLM.CNFireBaseData(btran2_patch = zeros(np))

        fire_li2014 = CLM.CNFireLi2014Data(
            forc_hdm     = [50.0],
            forc_lnfm    = [0.05],
            gdp_lf_col   = [10.0, 10.0],
            peatf_lf_col = [0.0, 0.0],
            abm_lf_col   = [6, 6],
        )

        patch = CLM.PatchData()
        patch.itype   = [2, 10, 5, 15]
        patch.column  = [1, 1, 2, 2]
        patch.wtcol   = [0.5, 0.5, 0.6, 0.4]
        patch.wtgcell = [0.5, 0.5, 0.6, 0.4]

        col = CLM.ColumnData()
        col.gridcell = [1, 1]
        col.wtgcell  = [1.0, 1.0]

        grc = CLM.GridcellData()
        grc.latdeg = [45.0]
        grc.lat    = [45.0 * pi / 180.0]

        soilstate = CLM.SoilStateData()
        soilstate.watsat_col   = fill(0.45, nc, nlevgrnd)
        soilstate.rootfr_patch = fill(1.0/nlevgrnd, np, nlevgrnd)
        soilstate.sucsat_col   = fill(200.0, nc, nlevgrnd)
        soilstate.bsw_col      = fill(5.0, nc, nlevgrnd)

        h2osoi_vol_col = fill(0.1, nc, nlevgrnd)

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

        cnveg_cs = CLM.CNVegCarbonStateData()
        cnveg_cs.totvegc_col              = [500.0, 400.0]
        cnveg_cs.rootc_col                = zeros(nc)
        cnveg_cs.leafc_col                = zeros(nc)
        cnveg_cs.deadstemc_col            = zeros(nc)
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
        cnveg_cs.deadstemc_patch          = [50.0, 0.0, 40.0, 0.0]

        decomp_cascade_con = CLM.DecompCascadeConData()
        decomp_cascade_con.is_litter = BitVector([true, true, true, false])
        decomp_cascade_con.is_cwd    = BitVector([false, false, false, true])
        decomp_cascade_con.spinup_factor = ones(ndecomp_pools)

        totlitc_col          = [200.0, 150.0]
        decomp_cpools_vr_col = fill(100.0, nc, nlevdecomp, ndecomp_pools)
        t_soi17cm_col        = fill(280.0, nc)

        forc_rh_grc    = [50.0]
        forc_wind_grc  = [3.0]
        forc_t_col     = fill(290.0, nc)
        forc_rain_col  = fill(0.0, nc)
        forc_snow_col  = fill(0.0, nc)
        prec60_patch   = fill(2.0e-5, np)
        prec10_patch   = fill(3.0e-5, np)
        prec30_patch   = fill(2.5e-5, np)
        rh30_patch     = fill(60.0, np)

        fsat_col = fill(0.1, nc)
        wf_col   = fill(0.3, nc)
        wf2_col  = fill(0.25, nc)

        croplive_patch = falses(np)

        mask_soilc        = trues(nc)
        mask_soilp        = trues(np)
        mask_exposedveg   = trues(np)
        mask_noexposedveg = falses(np)

        return (; pftcon, pftcon_li2014, cnfire_const, cnfire_params,
                  fire_data, fire_li2014, patch, col, grc, soilstate,
                  h2osoi_vol_col, cnveg_state, cnveg_cs, decomp_cascade_con,
                  totlitc_col, decomp_cpools_vr_col, t_soi17cm_col,
                  forc_rh_grc, forc_wind_grc, forc_t_col, forc_rain_col,
                  forc_snow_col, prec60_patch, prec10_patch, prec30_patch,
                  rh30_patch, fsat_col, wf_col, wf2_col, croplive_patch,
                  mask_soilc, mask_soilp, mask_exposedveg, mask_noexposedveg,
                  np, nc, ng, nlevgrnd, nlevdecomp, ndecomp_pools)
    end

    # Run one area method via the factory and return farea_burned (copy).
    function run_area(method::Symbol; kw...)
        d = make_fire_methods_area_data()
        CLM.dzsoi_decomp[] = [0.1]
        CLM.cnfire_area!(
            method,
            d.fire_li2014, d.pftcon_li2014, d.fire_data, d.cnfire_const,
            d.cnfire_params, d.pftcon,
            d.mask_soilc, d.mask_soilp, d.mask_exposedveg, d.mask_noexposedveg,
            1:d.nc, 1:d.np, d.patch, d.col, d.grc,
            d.soilstate, d.h2osoi_vol_col, d.cnveg_state, d.cnveg_cs,
            d.decomp_cascade_con, d.totlitc_col, d.decomp_cpools_vr_col,
            d.t_soi17cm_col;
            forc_rh_grc=d.forc_rh_grc, forc_wind_grc=d.forc_wind_grc,
            forc_t_col=d.forc_t_col, forc_rain_col=d.forc_rain_col,
            forc_snow_col=d.forc_snow_col,
            prec60_patch=d.prec60_patch, prec10_patch=d.prec10_patch,
            prec30_patch=d.prec30_patch, rh30_patch=d.rh30_patch,
            fsat_col=d.fsat_col, wf_col=d.wf_col, wf2_col=d.wf2_col,
            croplive_patch=d.croplive_patch,
            dt=1800.0, dayspyr=365.0, kmo=6, kda=15, mcsec=3600, nstep=10,
            nlevgrnd=d.nlevgrnd, nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools, kw...)
        return (d=d, farea=copy(d.cnveg_state.farea_burned_col),
                nfire=copy(d.cnveg_state.nfire_col),
                baf_peatf=copy(d.cnveg_state.baf_peatf_col))
    end

    # ================================================================
    # Factory: method-symbol mapping + lightning/popdens flags
    # ================================================================
    @testset "cnfire_method_symbol mapping" begin
        @test CLM.cnfire_method_symbol("nofire")        === :nofire
        @test CLM.cnfire_method_symbol("li2014qianfrc") === :li2014
        @test CLM.cnfire_method_symbol("li2016crufrc")  === :li2016
        @test CLM.cnfire_method_symbol("li2021gswpfrc") === :li2021
        @test CLM.cnfire_method_symbol("li2024gswpfrc") === :li2024
        @test CLM.cnfire_method_symbol("li2024crujra")  === :li2024
        @test CLM.cnfire_method_symbol(" LI2014QIANFRC ") === :li2014  # case/space
        @test_throws ErrorException CLM.cnfire_method_symbol("bogus")
    end

    @testset "need_lightning_and_popdens" begin
        @test CLM.need_lightning_and_popdens(:nofire) == false
        @test CLM.need_lightning_and_popdens(:li2014) == true
        @test CLM.need_lightning_and_popdens(:li2016) == true
        @test CLM.need_lightning_and_popdens(:li2021) == true
        @test CLM.need_lightning_and_popdens(:li2024) == true
    end

    # ================================================================
    # Each method: finite + physical burned area in [0,1]
    # ================================================================
    @testset "finite + physical burned area" begin
        for m in (:li2014, :li2016, :li2021, :li2024, :nofire)
            r = run_area(m)
            for c in 1:r.d.nc
                @test isfinite(r.farea[c])
                @test 0.0 <= r.farea[c] <= 1.0
                @test isfinite(r.nfire[c]) && r.nfire[c] >= 0.0
                @test isfinite(r.baf_peatf[c]) && r.baf_peatf[c] >= 0.0
            end
        end
    end

    # ================================================================
    # NoFire => exactly zero burned area / fbac everywhere
    # ================================================================
    @testset "nofire zeros everything" begin
        r = run_area(:nofire)
        d = r.d
        for c in 1:d.nc
            @test d.cnveg_state.farea_burned_col[c] == 0.0
            @test d.cnveg_state.baf_crop_col[c]     == 0.0
            @test d.cnveg_state.baf_peatf_col[c]    == 0.0
            @test d.cnveg_state.fbac_col[c]         == 0.0
            @test d.cnveg_state.fbac1_col[c]        == 0.0
            @test d.cnveg_state.cropf_col[c]        == 0.0
            @test d.cnveg_state.lfc_col[c]          == 0.0
        end
    end

    # ================================================================
    # Default cnfire_method reproduces Li2014 byte-for-byte
    # ================================================================
    @testset "default == li2014 (byte-identical)" begin
        r_factory = run_area(:li2014)

        # Direct li2014 call with the same fixture
        d = make_fire_methods_area_data()
        CLM.dzsoi_decomp[] = [0.1]
        CLM.cnfire_area_li2014!(
            d.fire_li2014, d.pftcon_li2014, d.fire_data, d.cnfire_const,
            d.cnfire_params, d.pftcon,
            d.mask_soilc, d.mask_soilp, d.mask_exposedveg, d.mask_noexposedveg,
            1:d.nc, 1:d.np, d.patch, d.col, d.grc,
            d.soilstate, d.h2osoi_vol_col, d.cnveg_state, d.cnveg_cs,
            d.decomp_cascade_con, d.totlitc_col, d.decomp_cpools_vr_col,
            d.t_soi17cm_col;
            forc_rh_grc=d.forc_rh_grc, forc_wind_grc=d.forc_wind_grc,
            forc_t_col=d.forc_t_col, forc_rain_col=d.forc_rain_col,
            forc_snow_col=d.forc_snow_col,
            prec60_patch=d.prec60_patch, prec10_patch=d.prec10_patch,
            fsat_col=d.fsat_col, wf_col=d.wf_col, wf2_col=d.wf2_col,
            dt=1800.0, dayspyr=365.0, kmo=6, kda=15, mcsec=3600, nstep=10,
            nlevgrnd=d.nlevgrnd, nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools)

        for c in 1:d.nc
            @test r_factory.farea[c] === d.cnveg_state.farea_burned_col[c]
            @test r_factory.nfire[c] === d.cnveg_state.nfire_col[c]
        end
    end

    # ================================================================
    # Li2016/2021/2024 differ from Li2014 where the formulas differ
    # ================================================================
    @testset "non-default methods differ from li2014" begin
        f2014 = run_area(:li2014).farea
        f2016 = run_area(:li2016).farea
        f2021 = run_area(:li2021).farea
        f2024 = run_area(:li2024).farea

        # Column 1 (natural veg, fire-spread active): the revised fire-occurrence
        # formulation (afuel/arh30 + GDP/pop factors) makes Li2016 differ from Li2014.
        @test f2016[1] != f2014[1]
        # Li2021 differs from Li2016 (rswf btran rescale + no bt clamp + cli form).
        @test f2021[1] != f2016[1]
        # Li2024 differs from Li2021 (lfwt/ig/fd rework + prec30).
        @test f2024[1] != f2021[1]
    end

    # ----------------------------------------------------------------
    # FLUX fixture (mirrors test_fire_li2014's make_fire_li2014_flux_data)
    # ----------------------------------------------------------------
    function make_fire_methods_flux_data(; np=4, nc=2, nlevdecomp=1,
                                          ndecomp_pools=4, n_litr=3)
        npft = 20
        pftcon = CLM.PftConFireBase(
            woody    = vcat(fill(1.0, 8), fill(0.0, 12)),
            cc_leaf  = fill(0.4, npft), cc_lstem = fill(0.2, npft),
            cc_dstem = fill(0.1, npft), cc_other = fill(0.3, npft),
            fm_leaf  = fill(0.6, npft), fm_lstem = fill(0.5, npft),
            fm_other = fill(0.4, npft), fm_root  = fill(0.3, npft),
            fm_lroot = fill(0.5, npft), fm_droot = fill(0.2, npft),
            lf_f     = fill(1.0/n_litr, npft, n_litr),
            fr_f     = fill(1.0/n_litr, npft, n_litr),
            smpso    = fill(-66000.0, npft), smpsc = fill(-275000.0, npft),
        )
        cnfire_const = CLM.CNFireConstData()

        patch = CLM.PatchData()
        patch.itype  = [2, 10, 5, 15]
        patch.column = [1, 1, 2, 2]
        patch.wtcol  = [0.5, 0.5, 0.6, 0.4]

        col = CLM.ColumnData(); col.gridcell = [1, 1]
        grc = CLM.GridcellData(); grc.latdeg = [45.0]; grc.lat = [45.0*pi/180.0]
        dgvs = CLM.DgvsFireData(nind_patch = fill(100.0, np))

        cnveg_state = CLM.CNVegStateData()
        cnveg_state.cropf_col        = [0.0, 0.4]
        cnveg_state.farea_burned_col = [1.0e-4, 1.0e-4]
        cnveg_state.fbac1_col        = zeros(nc)
        cnveg_state.fbac_col         = zeros(nc)
        cnveg_state.baf_crop_col     = [0.0, 2.0e-5]
        cnveg_state.baf_peatf_col    = zeros(nc)
        cnveg_state.trotr1_col       = zeros(nc)
        cnveg_state.trotr2_col       = zeros(nc)
        cnveg_state.dtrotr_col       = zeros(nc)
        cnveg_state.lfc_col          = zeros(nc)
        cnveg_state.lfc2_col         = zeros(nc)

        cnveg_cs = CLM.CNVegCarbonStateData()
        cnveg_cs.leafcmax_patch           = fill(0.0, np)
        cnveg_cs.leafc_patch              = [10.0, 5.0, 8.0, 3.0]
        cnveg_cs.leafc_storage_patch      = [1.0, 0.5, 0.8, 0.3]
        cnveg_cs.leafc_xfer_patch         = [0.5, 0.25, 0.4, 0.15]
        cnveg_cs.livestemc_patch          = [20.0, 0.0, 15.0, 0.0]
        cnveg_cs.livestemc_storage_patch  = [2.0, 0.0, 1.5, 0.0]
        cnveg_cs.livestemc_xfer_patch     = [1.0, 0.0, 0.8, 0.0]
        cnveg_cs.deadstemc_patch          = [50.0, 0.0, 40.0, 0.0]
        cnveg_cs.deadstemc_storage_patch  = [3.0, 0.0, 2.5, 0.0]
        cnveg_cs.deadstemc_xfer_patch     = [1.5, 0.0, 1.2, 0.0]
        cnveg_cs.frootc_patch             = [4.0, 2.0, 3.0, 1.0]
        cnveg_cs.frootc_storage_patch     = [0.5, 0.2, 0.3, 0.1]
        cnveg_cs.frootc_xfer_patch        = [0.2, 0.1, 0.15, 0.05]
        cnveg_cs.livecrootc_patch         = [10.0, 0.0, 8.0, 0.0]
        cnveg_cs.livecrootc_storage_patch = [1.0, 0.0, 0.8, 0.0]
        cnveg_cs.livecrootc_xfer_patch    = [0.5, 0.0, 0.4, 0.0]
        cnveg_cs.deadcrootc_patch         = [25.0, 0.0, 20.0, 0.0]
        cnveg_cs.deadcrootc_storage_patch = [1.5, 0.0, 1.2, 0.0]
        cnveg_cs.deadcrootc_xfer_patch    = [0.8, 0.0, 0.6, 0.0]
        cnveg_cs.gresp_storage_patch      = [0.2, 0.1, 0.15, 0.05]
        cnveg_cs.gresp_xfer_patch         = [0.1, 0.05, 0.08, 0.03]

        cnveg_cf = CLM.CNVegCarbonFluxData()
        for f in (:m_leafc_to_fire_patch, :m_leafc_storage_to_fire_patch,
                  :m_leafc_xfer_to_fire_patch, :m_livestemc_to_fire_patch,
                  :m_livestemc_storage_to_fire_patch, :m_livestemc_xfer_to_fire_patch,
                  :m_deadstemc_to_fire_patch, :m_deadstemc_storage_to_fire_patch,
                  :m_deadstemc_xfer_to_fire_patch, :m_frootc_to_fire_patch,
                  :m_frootc_storage_to_fire_patch, :m_frootc_xfer_to_fire_patch,
                  :m_livecrootc_to_fire_patch, :m_livecrootc_storage_to_fire_patch,
                  :m_livecrootc_xfer_to_fire_patch, :m_deadcrootc_to_fire_patch,
                  :m_deadcrootc_storage_to_fire_patch, :m_deadcrootc_xfer_to_fire_patch,
                  :m_gresp_storage_to_fire_patch, :m_gresp_xfer_to_fire_patch,
                  :m_leafc_to_litter_fire_patch, :m_leafc_storage_to_litter_fire_patch,
                  :m_leafc_xfer_to_litter_fire_patch, :m_livestemc_to_litter_fire_patch,
                  :m_livestemc_storage_to_litter_fire_patch, :m_livestemc_xfer_to_litter_fire_patch,
                  :m_livestemc_to_deadstemc_fire_patch, :m_deadstemc_to_litter_fire_patch,
                  :m_deadstemc_storage_to_litter_fire_patch, :m_deadstemc_xfer_to_litter_fire_patch,
                  :m_frootc_to_litter_fire_patch, :m_frootc_storage_to_litter_fire_patch,
                  :m_frootc_xfer_to_litter_fire_patch, :m_livecrootc_to_litter_fire_patch,
                  :m_livecrootc_storage_to_litter_fire_patch, :m_livecrootc_xfer_to_litter_fire_patch,
                  :m_livecrootc_to_deadcrootc_fire_patch, :m_deadcrootc_to_litter_fire_patch,
                  :m_deadcrootc_storage_to_litter_fire_patch, :m_deadcrootc_xfer_to_litter_fire_patch,
                  :m_gresp_storage_to_litter_fire_patch, :m_gresp_xfer_to_litter_fire_patch)
            setfield!(cnveg_cf, f, zeros(np))
        end
        cnveg_cf.fire_mortality_c_to_cwdc_col   = zeros(nc, nlevdecomp)
        cnveg_cf.m_decomp_cpools_to_fire_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        cnveg_cf.m_c_to_litr_fire_col           = zeros(nc, nlevdecomp, ndecomp_pools)

        cnveg_ns = CLM.CNVegNitrogenStateData()
        cnveg_ns.leafn_patch              = [0.4, 0.2, 0.32, 0.12]
        cnveg_ns.leafn_storage_patch      = [0.04, 0.02, 0.032, 0.012]
        cnveg_ns.leafn_xfer_patch         = [0.02, 0.01, 0.016, 0.006]
        cnveg_ns.livestemn_patch          = [0.4, 0.0, 0.3, 0.0]
        cnveg_ns.livestemn_storage_patch  = [0.04, 0.0, 0.03, 0.0]
        cnveg_ns.livestemn_xfer_patch     = [0.02, 0.0, 0.016, 0.0]
        cnveg_ns.deadstemn_patch          = [1.0, 0.0, 0.8, 0.0]
        cnveg_ns.deadstemn_storage_patch  = [0.06, 0.0, 0.05, 0.0]
        cnveg_ns.deadstemn_xfer_patch     = [0.03, 0.0, 0.024, 0.0]
        cnveg_ns.frootn_patch             = [0.16, 0.08, 0.12, 0.04]
        cnveg_ns.frootn_storage_patch     = [0.02, 0.008, 0.012, 0.004]
        cnveg_ns.frootn_xfer_patch        = [0.008, 0.004, 0.006, 0.002]
        cnveg_ns.livecrootn_patch         = [0.2, 0.0, 0.16, 0.0]
        cnveg_ns.livecrootn_storage_patch = [0.02, 0.0, 0.016, 0.0]
        cnveg_ns.livecrootn_xfer_patch    = [0.01, 0.0, 0.008, 0.0]
        cnveg_ns.deadcrootn_patch         = [0.5, 0.0, 0.4, 0.0]
        cnveg_ns.deadcrootn_storage_patch = [0.03, 0.0, 0.024, 0.0]
        cnveg_ns.deadcrootn_xfer_patch    = [0.016, 0.0, 0.012, 0.0]
        cnveg_ns.retransn_patch           = [0.05, 0.02, 0.04, 0.01]

        cnveg_nf = CLM.CNVegNitrogenFluxData()
        for f in (:m_leafn_to_fire_patch, :m_leafn_storage_to_fire_patch,
                  :m_leafn_xfer_to_fire_patch, :m_livestemn_to_fire_patch,
                  :m_livestemn_storage_to_fire_patch, :m_livestemn_xfer_to_fire_patch,
                  :m_deadstemn_to_fire_patch, :m_deadstemn_storage_to_fire_patch,
                  :m_deadstemn_xfer_to_fire_patch, :m_frootn_to_fire_patch,
                  :m_frootn_storage_to_fire_patch, :m_frootn_xfer_to_fire_patch,
                  :m_livecrootn_to_fire_patch, :m_livecrootn_storage_to_fire_patch,
                  :m_livecrootn_xfer_to_fire_patch, :m_deadcrootn_to_fire_patch,
                  :m_deadcrootn_storage_to_fire_patch, :m_deadcrootn_xfer_to_fire_patch,
                  :m_retransn_to_fire_patch, :m_leafn_to_litter_fire_patch,
                  :m_leafn_storage_to_litter_fire_patch, :m_leafn_xfer_to_litter_fire_patch,
                  :m_livestemn_to_litter_fire_patch, :m_livestemn_storage_to_litter_fire_patch,
                  :m_livestemn_xfer_to_litter_fire_patch, :m_livestemn_to_deadstemn_fire_patch,
                  :m_deadstemn_to_litter_fire_patch, :m_deadstemn_storage_to_litter_fire_patch,
                  :m_deadstemn_xfer_to_litter_fire_patch, :m_frootn_to_litter_fire_patch,
                  :m_frootn_storage_to_litter_fire_patch, :m_frootn_xfer_to_litter_fire_patch,
                  :m_livecrootn_to_litter_fire_patch, :m_livecrootn_storage_to_litter_fire_patch,
                  :m_livecrootn_xfer_to_litter_fire_patch, :m_livecrootn_to_deadcrootn_fire_patch,
                  :m_deadcrootn_to_litter_fire_patch, :m_deadcrootn_storage_to_litter_fire_patch,
                  :m_deadcrootn_xfer_to_litter_fire_patch, :m_retransn_to_litter_fire_patch)
            setfield!(cnveg_nf, f, zeros(np))
        end
        cnveg_nf.fire_mortality_n_to_cwdn_col   = zeros(nc, nlevdecomp)
        cnveg_nf.m_decomp_npools_to_fire_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        cnveg_nf.m_n_to_litr_fire_col           = zeros(nc, nlevdecomp, ndecomp_pools)

        soilbgc_cf = CLM.SoilBiogeochemCarbonFluxData()
        soilbgc_cf.somc_fire_col = zeros(nc)

        decomp_cascade_con = CLM.DecompCascadeConData()
        decomp_cascade_con.is_litter = BitVector([true, true, true, false])
        decomp_cascade_con.is_cwd    = BitVector([false, false, false, true])

        leaf_prof  = fill(1.0, np, nlevdecomp)
        froot_prof = fill(1.0, np, nlevdecomp)
        croot_prof = fill(1.0, np, nlevdecomp)
        stem_prof  = fill(1.0, np, nlevdecomp)
        totsomc    = fill(5000.0, nc)
        decomp_cpools_vr = fill(100.0, nc, nlevdecomp, ndecomp_pools)
        decomp_npools_vr = fill(5.0, nc, nlevdecomp, ndecomp_pools)
        somc_fire  = zeros(nc)

        mask_soilc = trues(nc); mask_soilp = trues(np)

        return (; pftcon, cnfire_const, patch, col, grc, dgvs, cnveg_state,
                  cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf, soilbgc_cf,
                  decomp_cascade_con, leaf_prof, froot_prof, croot_prof,
                  stem_prof, totsomc, decomp_cpools_vr, decomp_npools_vr,
                  somc_fire, mask_soilc, mask_soilp,
                  np, nc, nlevdecomp, ndecomp_pools, n_litr)
    end

    function run_flux(method::Symbol)
        d = make_fire_methods_flux_data()
        out = CLM.cnfire_fluxes_dispatch!(
            method,
            d.mask_soilc, d.mask_soilp, 1:d.nc, 1:d.np,
            d.cnfire_const, d.pftcon, d.patch, d.col, d.grc,
            d.dgvs, d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf,
            d.soilbgc_cf, d.decomp_cascade_con,
            d.leaf_prof, d.froot_prof, d.croot_prof, d.stem_prof,
            d.totsomc, d.decomp_cpools_vr, d.decomp_npools_vr, d.somc_fire;
            dt=1800.0, dayspyr=365.0, nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools, i_met_lit=1, i_litr_max=d.n_litr)
        return (d=d, out=out)
    end

    # ================================================================
    # Fluxes: finite, non-negative; mass-conserving decomp-pool destination
    # ================================================================
    @testset "fluxes finite + conserving" begin
        for m in (:li2014, :li2016, :li2021, :li2024, :nofire)
            r = run_flux(m)
            d = r.d
            for p in 1:d.np
                @test isfinite(d.cnveg_cf.m_leafc_to_fire_patch[p])
                @test d.cnveg_cf.m_leafc_to_fire_patch[p] >= 0.0
                @test d.cnveg_cf.m_gresp_storage_to_fire_patch[p] >= 0.0
            end
            # Decomp-pool combustion: burned C = decomp * f * cmb_factor.
            # Litter pools use 0.5; the CWD pool uses 0.25 (Li-family value).
            f = d.cnveg_state.farea_burned_col[1]
            baf_crop = d.cnveg_state.baf_crop_col[1]
            for l in 1:3
                @test d.cnveg_cf.m_decomp_cpools_to_fire_vr_col[1,1,l] ≈ 100.0 * f * 0.5 atol=1e-15
            end
            @test d.cnveg_cf.m_decomp_cpools_to_fire_vr_col[1,1,4] ≈ 100.0 * (f - baf_crop) * 0.25 atol=1e-15
        end
    end

    # ================================================================
    # NoFire fluxes: zero burned area => zero fire C/N fluxes
    # ================================================================
    @testset "nofire fluxes are zero" begin
        d = make_fire_methods_flux_data()
        # zero the burned area so the inherited base flux routine produces nothing
        d.cnveg_state.farea_burned_col .= 0.0
        d.cnveg_state.baf_crop_col     .= 0.0
        d.cnveg_state.baf_peatf_col    .= 0.0
        CLM.cnfire_fluxes_dispatch!(
            :nofire,
            d.mask_soilc, d.mask_soilp, 1:d.nc, 1:d.np,
            d.cnfire_const, d.pftcon, d.patch, d.col, d.grc,
            d.dgvs, d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf,
            d.soilbgc_cf, d.decomp_cascade_con,
            d.leaf_prof, d.froot_prof, d.croot_prof, d.stem_prof,
            d.totsomc, d.decomp_cpools_vr, d.decomp_npools_vr, d.somc_fire;
            dt=1800.0, dayspyr=365.0, nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools, i_met_lit=1, i_litr_max=d.n_litr)

        for p in 1:d.np
            @test d.cnveg_cf.m_leafc_to_fire_patch[p] == 0.0
            @test d.cnveg_nf.m_leafn_to_fire_patch[p] == 0.0
        end
        for l in 1:d.ndecomp_pools
            @test d.cnveg_cf.m_decomp_cpools_to_fire_vr_col[1,1,l] == 0.0
        end
        @test all(d.somc_fire .== 0.0)
    end

end
