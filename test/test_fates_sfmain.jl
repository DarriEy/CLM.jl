# test_fates_sfmain.jl
# Tests for FATES Tier-F Batch-11 SPITFIRE main fire driver (SFMainMod).
#
# Builds a synthetic site with one vegetated patch carrying a fuel object + a
# carbon litter pool + one woody tree cohort + a Nesterov fire-weather state,
# then runs the fire-model sub-steps and asserts sane bounds (area-burnt fraction
# in [0,1], non-negative intensity, mortality fraction in [0,1]) plus a couple of
# formula spot-checks (rate-of-spread responds to wind & to moisture).

@testset "FATES SFMainMod (SPITFIRE main driver)" begin
    fc = CLM.fuel_classes
    nfc = CLM.num_fuel_classes
    el = CLM.carbon12_element

    # ----------------------------------------------------------------------
    # Save + restore all the module globals we touch.
    # ----------------------------------------------------------------------
    old_numpft   = CLM.numpft[]
    old_ne       = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_epos     = copy(CLM.element_pos)
    old_woody    = copy(CLM.prt_params.woody)
    old_agbfrac  = copy(CLM.prt_params.allom_agb_frac)
    old_h2cd1    = copy(CLM.prt_params.allom_h2cd1)
    old_h2cd2    = copy(CLM.prt_params.allom_h2cd2)
    old_dmode    = copy(CLM.prt_params.allom_dmode)
    old_sf       = CLM.SFParams[]
    old_ed       = CLM.EDParams[]
    old_edpft    = CLM.EDPftvarcon_inst[]
    old_global   = CLM.prt_global[]
    old_sfmode   = CLM.hlm_spitfire_mode[]
    old_nofire   = CLM.hlm_sf_nofire_def[]
    old_scalar   = CLM.hlm_sf_scalar_lightning_def[]
    old_success  = CLM.hlm_sf_successful_ignitions_def[]
    old_anthro   = CLM.hlm_sf_anthro_ignitions_def[]
    old_master   = CLM.hlm_masterproc[]

    try
        npft = 1
        CLM.numpft[] = npft
        CLM.num_elements[] = 1
        empty!(CLM.element_list); push!(CLM.element_list, el)
        CLM.element_pos[el] = 1   # carbon12 lives at litter slot 1

        # --- PRT params: a single woody PFT, linear crown depth -------------
        CLM.prt_params.woody          = [CLM.itrue]
        CLM.prt_params.allom_agb_frac = [0.6]
        CLM.prt_params.allom_h2cd1    = [0.5]   # crown depth = 0.5 * height
        CLM.prt_params.allom_h2cd2    = [1.0]
        CLM.prt_params.allom_dmode    = [1]

        # --- SPITFIRE parameters (representative defaults) ------------------
        sfp = CLM.sf_params_type()
        sfp.SF_val_fdi_alpha     = 0.000337
        sfp.SF_val_miner_total   = 0.055
        sfp.SF_val_fuel_energy   = 18000.0
        sfp.SF_val_part_dens     = 513.0
        sfp.SF_val_miner_damp    = 0.41739
        sfp.SF_val_max_durat     = 240.0
        sfp.SF_val_durat_slope   = -11.06
        sfp.SF_val_drying_ratio  = 66000.0
        sfp.SF_val_fire_threshold = 50.0
        sfp.SF_val_SAV                = [13.0, 3.58, 0.98, 0.2, 66.0, 66.0]
        sfp.SF_val_FBD                = [15.4, 16.8, 19.6, 999.0, 4.0, 4.0]
        sfp.SF_val_min_moisture       = fill(0.0, nfc)
        sfp.SF_val_mid_moisture       = fill(1.0, nfc)
        sfp.SF_val_low_moisture_Coeff = fill(1.12, nfc)
        sfp.SF_val_low_moisture_Slope = fill(1.37, nfc)
        sfp.SF_val_mid_moisture_Coeff = fill(0.62, nfc)
        sfp.SF_val_mid_moisture_Slope = fill(0.13, nfc)
        CLM.SFParams[] = sfp

        # --- EDParams (ignitions / strikes) --------------------------------
        edp = CLM.ed_params_type()
        edp.ED_val_nignitions = 15.0
        edp.cg_strikes        = 0.2
        CLM.EDParams[] = edp

        # --- EDPftvarcon (fire PFT params) ---------------------------------
        ev = CLM.EDPftvarcon_type()
        ev.fire_alpha_SH = [0.1487]
        ev.bark_scaler   = [0.03]
        ev.crown_kill    = [0.775]
        CLM.EDPftvarcon_inst[] = ev

        # --- hlm flags: scalar-lightning SPITFIRE mode ---------------------
        CLM.hlm_masterproc[]                  = CLM.ifalse
        CLM.hlm_sf_nofire_def[]               = 0
        CLM.hlm_sf_scalar_lightning_def[]     = 1
        CLM.hlm_sf_successful_ignitions_def[] = 2
        CLM.hlm_sf_anthro_ignitions_def[]     = 3
        CLM.hlm_spitfire_mode[]               = 1   # scalar lightning

        # --- prt_global: carbon-only with leaf/sapw/struct organs ----------
        g = CLM.prt_global_type()
        CLM.ZeroGlobal!(g)
        g.hyp_name = "test_sfmain_carbon"
        g.num_vars     = 3
        g.num_bc_in    = 0
        g.num_bc_out   = 0
        g.num_bc_inout = 0
        g.state_descriptor = [CLM.state_descriptor_type() for _ in 1:g.num_vars]
        CLM.RegisterVarInGlobal!(g, 1, "leaf carbon",      "leaf_c", CLM.leaf_organ,   el, 1)
        CLM.RegisterVarInGlobal!(g, 2, "sapwood carbon",   "sapw_c", CLM.sapw_organ,   el, 1)
        CLM.RegisterVarInGlobal!(g, 3, "structure carbon", "str_c",  CLM.struct_organ, el, 1)
        CLM.prt_global[] = g

        # ------------------------------------------------------------------
        # Helper: build a vegetated patch carrying fuel + a carbon litter pool.
        # ------------------------------------------------------------------
        function make_patch(; area = 1000.0, tree_area = 800.0,
                              leaf_litter = 0.2, twig = 0.3, smb = 0.2, lgb = 0.1,
                              trunk = 0.5, grass = 0.1)
            p = CLM.fates_patch_type()
            p.area              = area
            p.total_tree_area   = tree_area
            p.total_grass_area  = 0.0
            p.nocomp_pft_label  = 1
            p.patchno           = 1
            p.livegrass         = grass
            p.fire              = 0
            p.frac_burnt        = 0.0
            p.tau_l             = 0.0
            p.tfc_ros           = 0.0
            p.ros_front         = 0.0
            p.ros_back          = 0.0
            p.fi                = 0.0
            p.fd                = 0.0
            p.scorch_ht         = zeros(Float64, CLM.maxpft)
            p.fuel              = CLM.fuel_type()

            # one carbon litter pool
            litt = CLM.litter_type()
            litt.ag_cwd     = [twig, smb, lgb, trunk]
            litt.leaf_fines = [leaf_litter]
            p.litter = [litt]

            # populate the fuel loadings + fractions from the litter
            CLM.update_loading!(p.fuel, leaf_litter, twig, smb, lgb, trunk, grass)
            CLM.sum_loading!(p.fuel)
            CLM.calculate_fractional_loading!(p.fuel)
            return p
        end

        # ------------------------------------------------------------------
        # Build a site with a single vegetated patch + a tree cohort.
        # ------------------------------------------------------------------
        site = CLM.ed_site_type()
        site.wind          = 120.0   # m/min
        site.NF            = 0.0
        site.NF_successful = 0.0
        site.fdi           = 0.0
        nw = CLM.nesterov_index()
        nw.fire_weather_index  = 60.0
        nw.effective_windspeed = 100.0  # m/min
        site.fireWeather = nw

        patch = make_patch()
        site.oldest_patch   = patch
        site.youngest_patch = patch

        # a woody tree cohort with set carbon state
        coh = CLM.fates_cohort_type(height = 20.0, pft = 1)
        coh.n   = 0.1
        coh.dbh = 30.0
        prt = CLM.prt_vartypes()
        CLM.InitPRTVartype!(prt)
        CLM.SetState!(prt, CLM.leaf_organ,   el, 2.0)
        CLM.SetState!(prt, CLM.sapw_organ,   el, 10.0)
        CLM.SetState!(prt, CLM.struct_organ, el, 20.0)
        coh.prt = prt
        patch.tallest = coh; patch.shortest = coh

        # ==================================================================
        # Run the sub-steps directly (skip update_fire_weather!, which needs
        # bc_in running-mean temperature; fireWeather state set above).
        # ==================================================================
        CLM.update_fuel_characteristics!(site)
        @test patch.fuel.non_trunk_loading > 0.0
        @test patch.fuel.bulk_density_notrunks > 0.0
        @test patch.fuel.SAV_notrunks > 0.0

        CLM.rate_of_spread!(site)
        @test isfinite(patch.ros_front)
        @test patch.ros_front >= 0.0
        @test patch.ros_back >= 0.0
        # backward ROS <= forward ROS (wind-driven attenuation, site.wind > 0)
        @test patch.ros_back <= patch.ros_front + 1e-12

        CLM.ground_fuel_consumption!(site)
        @test all(0.0 .<= patch.fuel.frac_burnt .<= 1.0)
        @test patch.fuel.frac_burnt[CLM.live_grass(fc)] <= 0.8 + 1e-12  # grass capped
        @test patch.tau_l >= 0.0
        @test patch.tau_l <= 8.0 + 1e-12   # residence time capped at 8 min
        @test patch.tfc_ros >= 0.0

        bc = CLM.bc_in_type()
        bc.lightning24 = [5.0]
        bc.pop_density = [10.0]
        bc.precip24_pa = [0.0]
        bc.relhumid24_pa = [40.0]
        bc.wind24_pa = [2.0]
        CLM.area_burnt_intensity!(site, bc)
        @test site.fdi >= 0.0 && site.fdi <= 1.0
        @test site.NF >= 0.0
        @test site.NF_successful >= 0.0
        @test isfinite(patch.fi) && patch.fi >= 0.0      # non-negative intensity
        @test patch.frac_burnt >= 0.0 && patch.frac_burnt <= 1.0  # area-burnt frac in [0,1]
        @test patch.fire == 0 || patch.fire == 1

        # Force a fire so the effects sub-steps have something to act on.
        patch.fire = 1
        patch.fi   = 5000.0      # strong intensity -> tall scorch height
        patch.tau_l = 5.0        # lethal-heating duration

        CLM.crown_scorching!(site)
        @test patch.scorch_ht[1] >= 0.0
        @test isfinite(patch.scorch_ht[1])

        CLM.crown_damage!(site)
        @test coh.fraction_crown_burned >= 0.0 && coh.fraction_crown_burned <= 1.0

        CLM.cambial_damage_kill!(site)
        @test coh.cambial_mort >= 0.0 && coh.cambial_mort <= 1.0

        CLM.post_fire_mortality!(site)
        @test coh.crownfire_mort >= 0.0 && coh.crownfire_mort <= 1.0
        @test coh.fire_mort >= 0.0 && coh.fire_mort <= 1.0  # mortality frac in [0,1]

        # ==================================================================
        # Formula spot-check #1: ROS increases with effective wind speed.
        # ==================================================================
        # Dry fuel (high Nesterov index) so the reaction intensity is nonzero.
        site_lo = CLM.ed_site_type(); site_lo.wind = 120.0
        nw_lo = CLM.nesterov_index(); nw_lo.fire_weather_index = 12000.0
        nw_lo.effective_windspeed = 10.0
        site_lo.fireWeather = nw_lo
        p_lo = make_patch(); site_lo.oldest_patch = p_lo; site_lo.youngest_patch = p_lo
        CLM.update_fuel_characteristics!(site_lo); CLM.rate_of_spread!(site_lo)

        site_hi = CLM.ed_site_type(); site_hi.wind = 120.0
        nw_hi = CLM.nesterov_index(); nw_hi.fire_weather_index = 12000.0
        nw_hi.effective_windspeed = 300.0
        site_hi.fireWeather = nw_hi
        p_hi = make_patch(); site_hi.oldest_patch = p_hi; site_hi.youngest_patch = p_hi
        CLM.update_fuel_characteristics!(site_hi); CLM.rate_of_spread!(site_hi)

        @test p_hi.ros_front > p_lo.ros_front   # more wind -> faster spread

        # ==================================================================
        # Formula spot-check #2: ROS decreases with higher fuel moisture
        # (drier fire weather index -> drier fuel -> faster spread).
        # ==================================================================
        site_dry = CLM.ed_site_type(); site_dry.wind = 120.0
        nw_dry = CLM.nesterov_index(); nw_dry.fire_weather_index = 15000.0  # very dry
        nw_dry.effective_windspeed = 100.0
        site_dry.fireWeather = nw_dry
        p_dry = make_patch(); site_dry.oldest_patch = p_dry; site_dry.youngest_patch = p_dry
        CLM.update_fuel_characteristics!(site_dry); CLM.rate_of_spread!(site_dry)

        site_wet = CLM.ed_site_type(); site_wet.wind = 120.0
        nw_wet = CLM.nesterov_index(); nw_wet.fire_weather_index = 8000.0  # damper
        nw_wet.effective_windspeed = 100.0
        site_wet.fireWeather = nw_wet
        p_wet = make_patch(); site_wet.oldest_patch = p_wet; site_wet.youngest_patch = p_wet
        CLM.update_fuel_characteristics!(site_wet); CLM.rate_of_spread!(site_wet)

        @test p_dry.fuel.average_moisture_notrunks < p_wet.fuel.average_moisture_notrunks
        @test p_dry.ros_front >= p_wet.ros_front   # drier fuel -> faster (or equal) spread

        # ==================================================================
        # No-fire mode: fire_model! leaves frac_burnt zeroed, runs no physics.
        # ==================================================================
        CLM.hlm_spitfire_mode[] = 0   # == nofire_def -> not > nofire
        site_nf = CLM.ed_site_type()
        p_nf = make_patch(); p_nf.frac_burnt = 0.7; p_nf.fire = 1
        site_nf.oldest_patch = p_nf; site_nf.youngest_patch = p_nf
        CLM.fire_model!(site_nf, bc)
        @test p_nf.frac_burnt == 0.0   # zeroed by the driver preamble
        @test p_nf.fire == 0

    finally
        CLM.numpft[] = old_numpft
        CLM.num_elements[] = old_ne
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
        copyto!(CLM.element_pos, old_epos)
        CLM.prt_params.woody          = old_woody
        CLM.prt_params.allom_agb_frac = old_agbfrac
        CLM.prt_params.allom_h2cd1    = old_h2cd1
        CLM.prt_params.allom_h2cd2    = old_h2cd2
        CLM.prt_params.allom_dmode    = old_dmode
        CLM.SFParams[]            = old_sf
        CLM.EDParams[]           = old_ed
        CLM.EDPftvarcon_inst[]   = old_edpft
        CLM.prt_global[]         = old_global
        CLM.hlm_spitfire_mode[]               = old_sfmode
        CLM.hlm_sf_nofire_def[]               = old_nofire
        CLM.hlm_sf_scalar_lightning_def[]     = old_scalar
        CLM.hlm_sf_successful_ignitions_def[] = old_success
        CLM.hlm_sf_anthro_ignitions_def[]     = old_anthro
        CLM.hlm_masterproc[]                  = old_master
    end
end
