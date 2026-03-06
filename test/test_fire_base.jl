@testset "Fire Base" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures
    # ----------------------------------------------------------------
    function make_fire_base_data(; np=3, nc=2, nlevgrnd=3, nlevdecomp=1,
                                   ndecomp_pools=4, n_litr=3)
        # --- PFT constants ---
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
            lf_f     = fill(1.0 / n_litr, npft, n_litr),
            fr_f     = fill(1.0 / n_litr, npft, n_litr),
            smpso    = fill(-66000.0, npft),
            smpsc    = fill(-275000.0, npft),
        )

        # --- Fire constants ---
        cnfire_const = CLM.CNFireConstData()

        # --- Fire data ---
        fire_data = CLM.CNFireBaseData(
            btran2_patch = zeros(np),
        )

        # --- Patch ---
        patch = CLM.PatchData()
        patch.itype  = [1, 9, 4]        # PFT types (0-based Fortran indices)
        patch.column = [1, 1, 2]          # column assignments
        patch.wtcol  = [0.5, 0.3, 0.2]   # weight relative to column

        # --- Column ---
        col = CLM.ColumnData()
        col.gridcell = [1, 1]

        # --- Gridcell ---
        grc = CLM.GridcellData()
        grc.latdeg = [45.0]   # boreal (> borealat=40)

        # --- Soil state ---
        soilstate = CLM.SoilStateData()
        soilstate.watsat_col  = fill(0.45, nc, nlevgrnd)
        soilstate.rootfr_patch = fill(1.0/nlevgrnd, np, nlevgrnd)
        soilstate.sucsat_col  = fill(200.0, nc, nlevgrnd)
        soilstate.bsw_col     = fill(5.0, nc, nlevgrnd)

        # --- Water state ---
        h2osoi_vol_col = fill(0.3, nc, nlevgrnd)

        # --- DGVS ---
        dgvs = CLM.DgvsFireData(
            nind_patch = fill(100.0, np),
        )

        # --- CN Veg State ---
        cnveg_state = CLM.CNVegStateData()
        cnveg_state.cropf_col         = [0.2, 0.0]
        cnveg_state.farea_burned_col  = [1.0e-5, 2.0e-5]
        cnveg_state.fbac1_col         = [0.0, 0.0]
        cnveg_state.fbac_col          = [1.0e-5, 2.0e-5]
        cnveg_state.baf_crop_col      = [2.0e-6, 0.0]
        cnveg_state.baf_peatf_col     = [1.0e-6, 5.0e-7]
        cnveg_state.trotr1_col        = [0.0, 0.0]
        cnveg_state.trotr2_col        = [0.0, 0.0]
        cnveg_state.dtrotr_col        = [0.0, 0.0]
        cnveg_state.lfc_col           = [0.0, 0.0]
        cnveg_state.lfc2_col          = [0.0, 0.0]

        # --- Carbon state ---
        cnveg_cs = CLM.CNVegCarbonStateData()
        cnveg_cs.leafcmax_patch               = fill(0.0, np)
        cnveg_cs.leafc_patch                  = [10.0, 5.0, 8.0]
        cnveg_cs.leafc_storage_patch          = [1.0, 0.5, 0.8]
        cnveg_cs.leafc_xfer_patch             = [0.5, 0.25, 0.4]
        cnveg_cs.livestemc_patch              = [20.0, 0.0, 15.0]
        cnveg_cs.livestemc_storage_patch      = [2.0, 0.0, 1.5]
        cnveg_cs.livestemc_xfer_patch         = [1.0, 0.0, 0.8]
        cnveg_cs.deadstemc_patch              = [50.0, 0.0, 40.0]
        cnveg_cs.deadstemc_storage_patch      = [3.0, 0.0, 2.5]
        cnveg_cs.deadstemc_xfer_patch         = [1.5, 0.0, 1.2]
        cnveg_cs.frootc_patch                 = [4.0, 2.0, 3.0]
        cnveg_cs.frootc_storage_patch         = [0.5, 0.2, 0.3]
        cnveg_cs.frootc_xfer_patch            = [0.2, 0.1, 0.15]
        cnveg_cs.livecrootc_patch             = [10.0, 0.0, 8.0]
        cnveg_cs.livecrootc_storage_patch     = [1.0, 0.0, 0.8]
        cnveg_cs.livecrootc_xfer_patch        = [0.5, 0.0, 0.4]
        cnveg_cs.deadcrootc_patch             = [25.0, 0.0, 20.0]
        cnveg_cs.deadcrootc_storage_patch     = [1.5, 0.0, 1.2]
        cnveg_cs.deadcrootc_xfer_patch        = [0.8, 0.0, 0.6]
        cnveg_cs.gresp_storage_patch          = [0.2, 0.1, 0.15]
        cnveg_cs.gresp_xfer_patch             = [0.1, 0.05, 0.08]

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
        cnveg_ns.leafn_patch                  = [0.4, 0.2, 0.32]
        cnveg_ns.leafn_storage_patch          = [0.04, 0.02, 0.032]
        cnveg_ns.leafn_xfer_patch             = [0.02, 0.01, 0.016]
        cnveg_ns.livestemn_patch              = [0.4, 0.0, 0.3]
        cnveg_ns.livestemn_storage_patch      = [0.04, 0.0, 0.03]
        cnveg_ns.livestemn_xfer_patch         = [0.02, 0.0, 0.016]
        cnveg_ns.deadstemn_patch              = [1.0, 0.0, 0.8]
        cnveg_ns.deadstemn_storage_patch      = [0.06, 0.0, 0.05]
        cnveg_ns.deadstemn_xfer_patch         = [0.03, 0.0, 0.024]
        cnveg_ns.frootn_patch                 = [0.16, 0.08, 0.12]
        cnveg_ns.frootn_storage_patch         = [0.02, 0.008, 0.012]
        cnveg_ns.frootn_xfer_patch            = [0.008, 0.004, 0.006]
        cnveg_ns.livecrootn_patch             = [0.2, 0.0, 0.16]
        cnveg_ns.livecrootn_storage_patch     = [0.02, 0.0, 0.016]
        cnveg_ns.livecrootn_xfer_patch        = [0.01, 0.0, 0.008]
        cnveg_ns.deadcrootn_patch             = [0.5, 0.0, 0.4]
        cnveg_ns.deadcrootn_storage_patch     = [0.03, 0.0, 0.024]
        cnveg_ns.deadcrootn_xfer_patch        = [0.016, 0.0, 0.012]
        cnveg_ns.retransn_patch               = [0.05, 0.02, 0.04]

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
        mask_exposedveg = trues(np)
        mask_noexposedveg = falses(np)

        return (pftcon=pftcon, cnfire_const=cnfire_const, fire_data=fire_data,
                patch=patch, col=col, grc=grc, soilstate=soilstate,
                h2osoi_vol_col=h2osoi_vol_col, dgvs=dgvs,
                cnveg_state=cnveg_state, cnveg_cs=cnveg_cs, cnveg_cf=cnveg_cf,
                cnveg_ns=cnveg_ns, cnveg_nf=cnveg_nf,
                soilbgc_cf=soilbgc_cf, decomp_cascade_con=decomp_cascade_con,
                leaf_prof=leaf_prof, froot_prof=froot_prof,
                croot_prof=croot_prof, stem_prof=stem_prof,
                totsomc=totsomc, decomp_cpools_vr=decomp_cpools_vr,
                decomp_npools_vr=decomp_npools_vr, somc_fire=somc_fire,
                mask_soilc=mask_soilc, mask_soilp=mask_soilp,
                mask_exposedveg=mask_exposedveg, mask_noexposedveg=mask_noexposedveg,
                np=np, nc=nc, nlevgrnd=nlevgrnd, nlevdecomp=nlevdecomp,
                ndecomp_pools=ndecomp_pools, n_litr=n_litr)
    end

    # ================================================================
    # Test: CNFireConstData construction
    # ================================================================
    @testset "CNFireConstData defaults" begin
        c = CLM.CNFireConstData()
        @test c.borealat == 40.0
        @test c.lfuel == 75.0
        @test c.ufuel == 650.0
        @test c.cmb_cmplt_fact_litter == 0.5
        @test c.cmb_cmplt_fact_cwd == 0.25
    end

    # ================================================================
    # Test: CNFireParams construction
    # ================================================================
    @testset "CNFireParams defaults" begin
        p = CLM.CNFireParams()
        @test p.prh30 == 0.0
        @test p.ignition_efficiency == 0.0
    end

    # ================================================================
    # Test: Li2014 root wetness
    # ================================================================
    @testset "Root wetness Li2014" begin
        d = make_fire_base_data()

        CLM.cnfire_calc_fire_root_wetness_li2014!(
            d.fire_data,
            d.mask_exposedveg,
            d.mask_noexposedveg,
            1:d.np,
            d.pftcon,
            d.patch,
            d.soilstate,
            d.h2osoi_vol_col,
            d.nlevgrnd
        )

        btran2 = d.fire_data.btran2_patch
        # All patches should have btran2 in [0, 1]
        for p in 1:d.np
            @test 0.0 <= btran2[p] <= 1.0
        end
        # With h2osoi_vol = 0.3, watsat = 0.45 → s_node = 0.667
        # Should give non-trivial wetness
        @test btran2[1] > 0.0
    end

    # ================================================================
    # Test: Li2021 root wetness
    # ================================================================
    @testset "Root wetness Li2021" begin
        d = make_fire_base_data()

        CLM.cnfire_calc_fire_root_wetness_li2021!(
            d.fire_data,
            d.mask_exposedveg,
            d.mask_noexposedveg,
            1:d.np,
            d.patch,
            d.soilstate,
            d.h2osoi_vol_col,
            d.nlevgrnd
        )

        btran2 = d.fire_data.btran2_patch
        for p in 1:d.np
            @test 0.0 <= btran2[p] <= 1.0
        end
        # With s_node = 0.667 and rootfr=1/3 per layer, sum = 0.667
        @test btran2[1] ≈ 0.3 / 0.45 atol=1e-10
    end

    # ================================================================
    # Test: Li2021 non-exposed patches get zero btran2
    # ================================================================
    @testset "Root wetness Li2021 non-exposed" begin
        d = make_fire_base_data()
        # Make all patches non-exposed
        mask_noexp = trues(d.np)
        mask_exp   = falses(d.np)
        d.fire_data.btran2_patch .= 999.0

        CLM.cnfire_calc_fire_root_wetness_li2021!(
            d.fire_data,
            mask_exp,
            mask_noexp,
            1:d.np,
            d.patch,
            d.soilstate,
            d.h2osoi_vol_col,
            d.nlevgrnd
        )

        for p in 1:d.np
            @test d.fire_data.btran2_patch[p] == 0.0
        end
    end

    # ================================================================
    # Test: CNFireFluxes runs without error
    # ================================================================
    @testset "CNFireFluxes basic" begin
        d = make_fire_base_data()

        mask_actfirec, mask_actfirep = CLM.cnfire_fluxes!(
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
    # Test: Fire fluxes are proportional to burned area
    # ================================================================
    @testset "CNFireFluxes proportionality" begin
        d = make_fire_base_data()

        # Run with a known farea_burned
        d.cnveg_state.farea_burned_col[1] = 1.0e-4
        d.cnveg_state.farea_burned_col[2] = 1.0e-4
        d.cnveg_state.baf_crop_col .= 0.0
        d.cnveg_state.cropf_col .= 0.0

        CLM.cnfire_fluxes!(
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

        f = 1.0e-4
        cc_leaf = 0.4  # from pftcon

        # m_leafc_to_fire = leafc * f * cc_leaf
        expected = 10.0 * f * cc_leaf  # patch 1: leafc=10, ivt=2
        @test d.cnveg_cf.m_leafc_to_fire_patch[1] ≈ expected atol=1e-15

        # Fine root fire emission should be zero (cc = 0.0 for roots)
        @test d.cnveg_cf.m_frootc_to_fire_patch[1] == 0.0
        @test d.cnveg_cf.m_livecrootc_to_fire_patch[1] == 0.0
        @test d.cnveg_cf.m_deadcrootc_to_fire_patch[1] == 0.0
    end

    # ================================================================
    # Test: Zero fire area produces zero fluxes
    # ================================================================
    @testset "CNFireFluxes zero fire area" begin
        d = make_fire_base_data()

        # Set all fire areas to zero
        d.cnveg_state.farea_burned_col .= 0.0
        d.cnveg_state.baf_crop_col .= 0.0
        d.cnveg_state.fbac_col .= 0.0
        d.cnveg_state.cropf_col .= 0.0

        CLM.cnfire_fluxes!(
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

        # All fire emission fluxes should be zero
        for p in 1:d.np
            @test d.cnveg_cf.m_leafc_to_fire_patch[p] == 0.0
            @test d.cnveg_cf.m_livestemc_to_fire_patch[p] == 0.0
            @test d.cnveg_nf.m_leafn_to_fire_patch[p] == 0.0
        end
    end

    # ================================================================
    # Test: Peat fire (somc_fire) boreal vs non-boreal
    # ================================================================
    @testset "Peat fire somc_fire" begin
        d = make_fire_base_data()

        # Set gridcell latitude to tropical (below borealat=40)
        d.grc.latdeg[1] = 10.0
        d.cnveg_state.baf_peatf_col[1] = 1.0e-5
        d.cnveg_state.baf_peatf_col[2] = 1.0e-5

        CLM.cnfire_fluxes!(
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

        # Non-boreal: somc_fire = totsomc * baf_peatf * 6.0/33.9
        expected_nonboreal = 5000.0 * 1.0e-5 * 6.0 / 33.9
        @test d.somc_fire[1] ≈ expected_nonboreal atol=1e-10

        # Boreal (col 2 → gridcell 1, but latdeg=10 so also non-boreal)
        @test d.somc_fire[2] ≈ expected_nonboreal atol=1e-10

        # Now test boreal
        d2 = make_fire_base_data()
        d2.grc.latdeg[1] = 55.0  # boreal
        d2.cnveg_state.baf_peatf_col .= 1.0e-5

        CLM.cnfire_fluxes!(
            d2.mask_soilc, d2.mask_soilp,
            1:d2.nc, 1:d2.np,
            d2.cnfire_const, d2.pftcon, d2.patch, d2.col, d2.grc,
            d2.dgvs, d2.cnveg_state, d2.cnveg_cs, d2.cnveg_cf,
            d2.cnveg_ns, d2.cnveg_nf,
            d2.soilbgc_cf, d2.decomp_cascade_con,
            d2.leaf_prof, d2.froot_prof, d2.croot_prof, d2.stem_prof,
            d2.totsomc, d2.decomp_cpools_vr, d2.decomp_npools_vr, d2.somc_fire;
            dt=1800.0, dayspyr=365.0,
            nlevdecomp=d2.nlevdecomp, ndecomp_pools=d2.ndecomp_pools,
            i_met_lit=1, i_litr_max=d2.n_litr
        )

        # Boreal: somc_fire = baf_peatf * 2200
        expected_boreal = 1.0e-5 * 2.2e3
        @test d2.somc_fire[1] ≈ expected_boreal atol=1e-10
    end

    # ================================================================
    # Test: Decomp fire fluxes for litter and CWD pools
    # ================================================================
    @testset "Decomp pool fire fluxes" begin
        d = make_fire_base_data()

        f_burned = 5.0e-5
        d.cnveg_state.farea_burned_col .= f_burned
        d.cnveg_state.baf_crop_col .= 0.0

        CLM.cnfire_fluxes!(
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

        # Litter pool (pools 1-3): decomp_cpools_vr * f * 0.5
        for l in 1:3
            expected = 100.0 * f_burned * 0.5
            @test d.cnveg_cf.m_decomp_cpools_to_fire_vr_col[1, 1, l] ≈ expected atol=1e-12
        end

        # CWD pool (pool 4): decomp_cpools_vr * (f - baf_crop) * 0.25
        expected_cwd = 100.0 * f_burned * 0.25
        @test d.cnveg_cf.m_decomp_cpools_to_fire_vr_col[1, 1, 4] ≈ expected_cwd atol=1e-12

        # Same for N
        for l in 1:3
            expected_n = 5.0 * f_burned * 0.5
            @test d.cnveg_nf.m_decomp_npools_to_fire_vr_col[1, 1, l] ≈ expected_n atol=1e-12
        end
    end

end
