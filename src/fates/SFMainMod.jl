# SFMainMod.jl
# Julia port of FATES src/fates/fire/SFMainMod.F90
#
# The SPITFIRE main fire routine. Code originally developed by Allan Spessa &
# Rosie Fisher as part of the NERC-QUEST project. Runs the daily fire model on a
# FATES site: per-patch fuel characterization -> fire weather / danger index ->
# rate of spread (Rothermel 1972 / Thonicke et al. 2010) -> ground fuel
# consumption -> fire intensity & area burned -> fire effects (crown scorch,
# crown damage, cambial damage / cambial kill, post-fire tree mortality).
#
# Julia port notes:
#   * Fortran subroutines -> Julia bang functions with the SAME names
#     (fire_model!, update_fire_weather!, update_fuel_characteristics!,
#     rate_of_spread!, ground_fuel_consumption!, area_burnt_intensity!,
#     crown_scorching!, crown_damage!, cambial_damage_kill!,
#     post_fire_mortality!).
#   * Operates on the merged ed_site_type / fates_patch_type / fates_cohort_type
#     plus the `fuel_type` (FatesFuelMod) and the `fire_weather` /
#     `nesterov_index` objects (SFFireWeatherMod / SFNesterovMod).
#   * Patch linked list walked via .older / .younger; cohort list via
#     .tallest / .shorter, exactly as the Fortran pointer walks.
#   * Parameters are read from the live singletons: sf_params() (SFParamsMod),
#     ed_params() (EDParamsMod), edpftvarcon_inst() (EDPftvarcon), prt_params
#     (PRTParametersMod). `fates_r8` -> Float64.
#   * The SPITFIRE spread / intensity / mortality formulas are preserved verbatim
#     (Thonicke et al. 2010 equation numbering kept in comments).
#   * Standalone — NOT added to CLMInstances or any dual-copied struct.
#
# Deps: EDTypesMod (ed_site_type, area, CalculateTreeGrassAreaSite!),
# FatesPatchMod (fates_patch_type, UpdateLiveGrass!), FatesCohortMod
# (fates_cohort_type), FatesFuelMod (fuel_type + fuel bang-functions),
# SFFireWeatherMod / SFNesterovMod (fire_weather, nesterov_index, update_index!,
# update_effective_windspeed!), FatesFuelClassesMod (fuel_classes + accessors,
# num_fuel_classes), SFParamsMod (sf_params), EDParamsMod (ed_params, maxpft),
# EDPftvarcon (edpftvarcon_inst), PRTParametersMod (prt_params), PRTGenericMod
# (element_pos, carbon12_element, leaf_organ, sapw_organ, struct_organ, GetState),
# FatesAllometryMod (CrownDepth), FatesConstantsMod (constants),
# FatesInterfaceTypesMod (hlm_* flags, bc_in_type, numpft).

# ---------------------------------------------------------------------------
# fire_model!  — top-level daily fire driver
# ---------------------------------------------------------------------------
"""
    fire_model!(currentSite::ed_site_type, bc_in::bc_in_type)

Run the daily SPITFIRE fire model on a site (Fortran `fire_model`). Zeroes the
per-patch fire diagnostics, then — if SPITFIRE is enabled
(`hlm_spitfire_mode[] > hlm_sf_nofire_def[]`) — runs the full chain: fire
weather, fuel characteristics, rate of spread, ground fuel consumption, area
burnt + intensity, crown scorching, crown damage, cambial damage / kill, and
post-fire mortality.
"""
function fire_model!(currentSite::ed_site_type, bc_in::bc_in_type)
    # zero fire things
    currentPatch = currentSite.youngest_patch
    while currentPatch !== nothing
        currentPatch.frac_burnt = 0.0
        currentPatch.fire = 0
        currentPatch = currentPatch.older
    end

    if hlm_spitfire_mode[] > hlm_sf_nofire_def[]
        update_fire_weather!(currentSite, bc_in)
        update_fuel_characteristics!(currentSite)
        rate_of_spread!(currentSite)
        ground_fuel_consumption!(currentSite)
        area_burnt_intensity!(currentSite, bc_in)
        crown_scorching!(currentSite)
        crown_damage!(currentSite)
        cambial_damage_kill!(currentSite)
        post_fire_mortality!(currentSite)
    end

    return currentSite
end

# ---------------------------------------------------------------------------
# update_fire_weather!  (Fortran UpdateFireWeather)
# ---------------------------------------------------------------------------
"""
    update_fire_weather!(currentSite::ed_site_type, bc_in::bc_in_type)

Update the site's fire-weather index and calculate the effective wind speed based
on vegetation characteristics (Fortran `UpdateFireWeather`). Uses tree and grass
fractions averaged over the whole site to prevent extreme divergence. Boundary
conditions (temperature, precipitation, relative humidity, wind) are taken from
the oldest vegetated patch.
"""
function update_fire_weather!(currentSite::ed_site_type, bc_in::bc_in_type)
    tfrz = t_water_freeze_k_1atm

    currentPatch = currentSite.oldest_patch

    # If the oldest patch is a bareground patch (i.e. nocomp mode is on) use the
    # first vegetated patch for the iofp index (i.e. the next younger patch)
    if currentPatch.nocomp_pft_label == nocomp_bareground
        currentPatch = currentPatch.younger
    end

    iofp = currentPatch.patchno
    temp_C = GetMean(currentPatch.tveg24) - tfrz
    precip = bc_in.precip24_pa[iofp] * sec_per_day
    rh = bc_in.relhumid24_pa[iofp]
    wind = bc_in.wind24_pa[iofp]

    # convert to m/min
    currentSite.wind = wind * sec_per_min

    # update fire weather index
    update_index!(currentSite.fireWeather, temp_C, precip, rh, wind)

    # calculate site-level tree, grass, and bare fraction
    tree_fraction, grass_fraction, bare_fraction =
        CalculateTreeGrassAreaSite!(currentSite)

    # update effective wind speed
    update_effective_windspeed!(currentSite.fireWeather, wind * sec_per_min,
                                tree_fraction, grass_fraction, bare_fraction)

    return currentSite
end

# ---------------------------------------------------------------------------
# update_fuel_characteristics!  (Fortran UpdateFuelCharacteristics)
# ---------------------------------------------------------------------------
"""
    update_fuel_characteristics!(currentSite::ed_site_type)

Update fuel characteristics on each (vegetated) patch of the site (Fortran
`UpdateFuelCharacteristics`): live grass, fuel loading from litter, fractional
loading, fuel moisture, and the non-trunk geometric properties (bulk density and
surface-area-to-volume ratio).
"""
function update_fuel_characteristics!(currentSite::ed_site_type)
    p = sf_params()

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        if currentPatch.nocomp_pft_label != nocomp_bareground
            # calculate live grass [kgC/m2]
            UpdateLiveGrass!(currentPatch)

            # update fuel loading [kgC/m2]
            litter = currentPatch.litter[element_pos[carbon12_element]]
            update_loading!(currentPatch.fuel, sum(litter.leaf_fines),
                litter.ag_cwd[1], litter.ag_cwd[2], litter.ag_cwd[3],
                litter.ag_cwd[4], currentPatch.livegrass)

            # sum up fuel classes and calculate fractional loading for each
            sum_loading!(currentPatch.fuel)
            calculate_fractional_loading!(currentPatch.fuel)

            # calculate fuel moisture [m3/m3]
            update_fuel_moisture!(currentPatch.fuel, p.SF_val_SAV,
                p.SF_val_drying_ratio, currentSite.fireWeather)

            # calculate geometric properties
            average_bulk_density_notrunks!(currentPatch.fuel, p.SF_val_FBD)
            average_sav_notrunks!(currentPatch.fuel, p.SF_val_SAV)
        end
        currentPatch = currentPatch.younger
    end

    return currentSite
end

# ---------------------------------------------------------------------------
# rate_of_spread!  (Fortran rate_of_spread)
# ---------------------------------------------------------------------------
"""
    rate_of_spread!(currentSite::ed_site_type)

Compute the forward and backward rate of spread (Rothermel 1972 / Thonicke et al.
2010) for each patch and store on `patch.ros_front` / `patch.ros_back` [m/min]
(Fortran `rate_of_spread`).
"""
function rate_of_spread!(currentSite::ed_site_type)
    p = sf_params()
    q_dry = 581.0  # heat of pre-ignition of dry fuels (kJ/kg)

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        if currentPatch.nocomp_pft_label != nocomp_bareground &&
           currentPatch.fuel.non_trunk_loading > nearzero

            fuel = currentPatch.fuel

            # remove mineral content from net fuel load per Thonicke 2010 for ir
            fuel.non_trunk_loading = fuel.non_trunk_loading * (1.0 - p.SF_val_miner_total)

            # ----start spreading---

            # beta = packing ratio (unitless): fraction of fuel array volume
            # occupied by fuel, or compactness of the fuel bed
            beta = fuel.bulk_density_notrunks / p.SF_val_part_dens

            # Equation A6 in Thonicke et al. 2010 - packing ratio (unitless)
            if fuel.SAV_notrunks < nearzero
                beta_op = 0.0
            else
                beta_op = 0.200395 * (fuel.SAV_notrunks^(-0.8189))
            end

            if beta_op < nearzero
                beta_ratio = 0.0
            else
                beta_ratio = beta / beta_op  # unitless
            end

            # ---heat of pre-ignition--- Equation A4 in Thonicke et al. 2010
            # q_ig in kJ/kg
            q_ig = q_dry + 2594.0 * fuel.average_moisture_notrunks

            # ---effective heating number--- Equation A3 in Thonicke et al. 2010
            eps = exp(-4.528 / fuel.SAV_notrunks)
            # Equation A7 (Rothermel 1972 eqn 49)
            b = 0.15988 * (fuel.SAV_notrunks^0.54)
            # Equation A8 (Rothermel 1972 eqn 48)
            c = 7.47 * (exp(-0.8711 * (fuel.SAV_notrunks^0.55)))
            # Equation A9 (Rothermel 1972 eqn 50)
            e = 0.715 * (exp(-0.01094 * fuel.SAV_notrunks))

            # Equation A5 in Thonicke et al. 2010 - phi_wind (unitless)
            # convert effective windspeed from m/min to ft/min for Rothermel ROS
            phi_wind = c * ((3.281 * currentSite.fireWeather.effective_windspeed)^b) *
                       (beta_ratio^(-e))

            # ---propagating flux--- Equation A2 Thonicke / Eq.42 Rothermel 1972
            xi = (exp((0.792 + 3.7597 * (fuel.SAV_notrunks^0.5)) * (beta + 0.1))) /
                 (192.0 + 7.9095 * fuel.SAV_notrunks)

            # ---reaction intensity--- Equation in table A1 Thonicke et al. 2010
            a = 8.9033 * (fuel.SAV_notrunks^(-0.7913))
            a_beta = exp(a * (1.0 - beta_ratio))  # dummy for reaction_v_opt

            # reaction velocity (per min): max = Rothermel 1972 eqn 36, opt = eqn 38
            reaction_v_max = 1.0 / (0.0591 + 2.926 * (fuel.SAV_notrunks^(-1.5)))
            reaction_v_opt = reaction_v_max * (beta_ratio^a) * a_beta

            # mw_weight = relative fuel moisture / fuel moisture of extinction
            mw_weight = fuel.average_moisture_notrunks / fuel.MEF_notrunks

            # moist_damp (unitless) - Equation in table A1 Thonicke et al. 2010
            moist_damp = max(0.0, (1.0 - (2.59 * mw_weight) + (5.11 * (mw_weight^2.0)) -
                                   (3.52 * (mw_weight^3.0))))

            # ir = reaction intensity in kJ/m2/min. non_trunk_loading converted
            # from kgC/m2 to kgBiomass/m2 for ir calculation
            ir = reaction_v_opt * (fuel.non_trunk_loading / 0.45) *
                 p.SF_val_fuel_energy * moist_damp * p.SF_val_miner_damp

            if (fuel.bulk_density_notrunks <= 0.0) || (eps <= 0.0) || (q_ig <= 0.0)
                currentPatch.ros_front = 0.0
            else  # Equation 9 Thonicke et al. 2010 - forward ROS in m/min
                currentPatch.ros_front = (ir * xi * (1.0 + phi_wind)) /
                    (fuel.bulk_density_notrunks * eps * q_ig)
            end

            # Equation 10 in Thonicke et al. 2010 - backward ROS (Can FBP 1992),
            # m/min; backward ROS wind not changed by vegetation
            currentPatch.ros_back = currentPatch.ros_front *
                exp(-0.012 * currentSite.wind)
        end
        currentPatch = currentPatch.younger
    end

    return currentSite
end

# ---------------------------------------------------------------------------
# ground_fuel_consumption!  (Fortran ground_fuel_consumption)
# ---------------------------------------------------------------------------
"""
    ground_fuel_consumption!(currentSite::ed_site_type)

Calculate the hypothetic fuel consumed by the fire (Fortran
`ground_fuel_consumption`): per-fuel-class burnt fraction (Eq. B1 of Thonicke et
al. 2010), the lethal-heating residence time `patch.tau_l` [min] (Peterson & Ryan
1986), and the total fuel consumed by the spreading flaming front `patch.TFC_ROS`
(stored on `patch.tfc_ros`) [kgC/m2 of burned area].
"""
function ground_fuel_consumption!(currentSite::ed_site_type)
    p = sf_params()

    tr_sf = trunks(fuel_classes)
    tw_sf = twigs(fuel_classes)
    dl_sf = dead_leaves(fuel_classes)
    lg_sf = live_grass(fuel_classes)

    tau_b = zeros(Float64, num_fuel_classes)      # lethal heating rates per class [min]
    fc_ground = zeros(Float64, num_fuel_classes)  # fuel consumed per area burned [kgC/m2]

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        if currentPatch.nocomp_pft_label != nocomp_bareground
            fuel = currentPatch.fuel

            fuel.frac_burnt[1:num_fuel_classes] .= 1.0

            # Calculate fraction of litter burnt for all classes.
            # Equation B1 in Thonicke et al. 2010
            for cc in 1:num_fuel_classes  # all pools, even if those pools don't exist
                moist = fuel.effective_moisture[cc]
                # 1. Very dry litter
                if moist <= p.SF_val_min_moisture[cc]
                    fuel.frac_burnt[cc] = 1.0
                end
                # 2. Low to medium moistures
                if moist > p.SF_val_min_moisture[cc] && moist <= p.SF_val_mid_moisture[cc]
                    fuel.frac_burnt[cc] = max(0.0, min(1.0,
                        p.SF_val_low_moisture_Coeff[cc] -
                        p.SF_val_low_moisture_Slope[cc] * moist))
                else
                    # For medium to high moistures.
                    if moist > p.SF_val_mid_moisture[cc] && moist <= 1.0
                        fuel.frac_burnt[cc] = max(0.0, min(1.0,
                            p.SF_val_mid_moisture_Coeff[cc] -
                            p.SF_val_mid_moisture_Slope[cc] * moist))
                    end
                end
                # Very wet litter
                if moist >= 1.0  # this shouldn't happen?
                    fuel.frac_burnt[cc] = 0.0
                end
            end

            # we can't ever kill -all- of the grass.
            fuel.frac_burnt[lg_sf] = min(0.8, fuel.frac_burnt[lg_sf])

            # reduce burnt amount for mineral content.
            fuel.frac_burnt[1:num_fuel_classes] .=
                fuel.frac_burnt[1:num_fuel_classes] .* (1.0 - p.SF_val_miner_total)

            # ---Calculate amount of fuel burnt.---
            litt_c = currentPatch.litter[element_pos[carbon12_element]]
            for cc in tw_sf:tr_sf
                fc_ground[cc] = fuel.frac_burnt[cc] * litt_c.ag_cwd[cc]
            end
            fc_ground[dl_sf] = fuel.frac_burnt[dl_sf] * sum(litt_c.leaf_fines)
            fc_ground[lg_sf] = fuel.frac_burnt[lg_sf] * currentPatch.livegrass

            # Cambial-kill residence time follows Peterson & Ryan (1986). tau_b is
            # the duration of the lethal heating. The /10 converts kgC/m2 to
            # gC/cm2 as in the P&R paper.
            for cc in 1:num_fuel_classes
                tau_b[cc] = 39.4 *
                    (fuel.frac_loading[cc] * fuel.non_trunk_loading / 0.45 / 10.0) *
                    (1.0 - ((1.0 - fuel.frac_burnt[cc])^0.5))
            end
            tau_b[tr_sf] = 0.0
            # Cap the residence time to 8 mins (P&R 1986).
            currentPatch.tau_l = min(8.0, sum(tau_b))

            # ---overall fuel consumed by spreading fire--- ignore 1000hr (trunk)
            # fuels; just interested in fuels affecting ROS
            currentPatch.tfc_ros = sum(fc_ground) - fc_ground[tr_sf]
        end
        currentPatch = currentPatch.younger
    end

    return currentSite
end

# ---------------------------------------------------------------------------
# area_burnt_intensity!  (Fortran area_burnt_intensity)
# ---------------------------------------------------------------------------
"""
    area_burnt_intensity!(currentSite::ed_site_type, bc_in::bc_in_type)

Calculate the fire danger index, number of ignitions, fire duration, area burnt
and fire intensity for each patch (Fortran `area_burnt_intensity`). Sets
`site.fdi`, `site.NF`, `site.NF_successful`, and per patch `patch.fd` [min],
`patch.frac_burnt` [0-1], `patch.fi` [kW/m] and `patch.fire` [0/1].
"""
function area_burnt_intensity!(currentSite::ed_site_type, bc_in::bc_in_type)
    p = sf_params()
    edp = ed_params()

    pot_hmn_ign_counts_alpha = 0.0035  # Li et al. 2012 (#/person/month)
    km2_to_m2 = 1000000.0              # km2 -> m2
    m_per_min__to__km_per_hour = 0.06  # m/min -> km/hr
    forest_grassland_lengthtobreadth_threshold = 0.55  # tree-cover threshold

    # ---initialize site parameters to zero---
    currentSite.NF_successful = 0.0

    # Equation 7 Venevsky et al. 2002 / Eq.8 Thonicke et al. 2010
    if hlm_spitfire_mode[] == hlm_sf_successful_ignitions_def[]
        currentSite.fdi = 1.0  # "SUCCESSFUL IGNITION" data: force extreme potential
        cloud_to_ground_strikes = 1.0  # use 100% incoming observed ignitions
    else  # USING LIGHTNING DATA
        currentSite.fdi = 1.0 - exp(-p.SF_val_fdi_alpha *
                                    currentSite.fireWeather.fire_weather_index)
        cloud_to_ground_strikes = edp.cg_strikes
    end

    currentPatch = currentSite.oldest_patch

    # If the oldest patch is a bareground patch use the first vegetated patch
    if currentPatch.nocomp_pft_label == nocomp_bareground
        currentPatch = currentPatch.younger
    end

    # NF = number of lightning strikes per day per km2 scaled by cg strikes
    iofp = currentPatch.patchno
    if hlm_spitfire_mode[] == hlm_sf_scalar_lightning_def[]
        currentSite.NF = edp.ED_val_nignitions * years_per_day * cloud_to_ground_strikes
    else  # use external daily lightning ignition data
        currentSite.NF = bc_in.lightning24[iofp] * cloud_to_ground_strikes
    end

    # Calculate anthropogenic ignitions according to Li et al. (2012)
    if hlm_spitfire_mode[] == hlm_sf_anthro_ignitions_def[]
        anthro_ign_count = pot_hmn_ign_counts_alpha * 6.8 *
            bc_in.pop_density[iofp]^0.43 / 30.0
        currentSite.NF = currentSite.NF + anthro_ign_count
    end

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        if currentPatch.nocomp_pft_label != nocomp_bareground

            # ---initialize patch parameters to zero---
            currentPatch.fi = 0.0
            currentPatch.fire = 0
            currentPatch.fd = 0.0
            currentPatch.frac_burnt = 0.0

            if currentSite.NF > 0.0
                # Equation 14 Thonicke et al. 2010 - fire duration in minutes
                currentPatch.fd = (p.SF_val_max_durat + 1.0) /
                    (1.0 + p.SF_val_max_durat * exp(p.SF_val_durat_slope * currentSite.fdi))

                tree_fraction_patch = currentPatch.total_tree_area / currentPatch.area

                ews = currentSite.fireWeather.effective_windspeed
                if (ews * m_per_min__to__km_per_hour) < 1.0  # 16.67m/min = 1km/hr
                    lb = 1.0
                else
                    if tree_fraction_patch > forest_grassland_lengthtobreadth_threshold
                        # EQ 79 forest fuels (CFFBPS Ont.Inf.Rep. ST-X-3, 1992)
                        lb = (1.0 + (8.729 *
                            ((1.0 - (exp(-0.03 * m_per_min__to__km_per_hour * ews)))^2.155)))
                    else  # EQ 80 grass fuels (CFFBPS, w/ Wotton et al. 2009 errata)
                        lb = (1.1 * ((m_per_min__to__km_per_hour * ews)^0.464))
                    end
                end

                # ---- length of major axis----
                db = currentPatch.ros_back * currentPatch.fd  # m
                df = currentPatch.ros_front * currentPatch.fd  # m

                # --- area burnt---
                if lb > 0.0
                    # size of fire = Eq.14 Arora & Boer JGR 2005 (area of ellipse)
                    size_of_fire = ((pi_const / (4.0 * lb)) * ((df + db)^2.0))

                    # AB = daily area burnt (m2 per km2 per day)
                    AB = size_of_fire * currentSite.NF * currentSite.fdi

                    # frac_burnt: area burned per area patch per day
                    currentPatch.frac_burnt = (min(0.99, AB / km2_to_m2))
                else
                    currentPatch.frac_burnt = 0.0
                end

                ROS = currentPatch.ros_front / 60.0   # m/min to m/sec
                W = currentPatch.tfc_ros / 0.45       # kgC/m2 -> kgbiomass/m2

                # EQ 15 Thonicke et al. 2010 - fire intensity kW/m
                currentPatch.fi = p.SF_val_fuel_energy * W * ROS

                # 'decide_fire': track fires above the kW/m energy threshold
                if currentPatch.fi > p.SF_val_fire_threshold
                    currentPatch.fire = 1  # Fire...    :D
                    currentSite.NF_successful = currentSite.NF_successful +
                        currentSite.NF * currentSite.fdi * currentPatch.area / area
                else
                    currentPatch.fire = 0  # No fire... :-/
                    currentPatch.fd = 0.0
                    currentPatch.frac_burnt = 0.0
                end
            end  # NF ignitions check
        end  # nocomp_pft_label check
        currentPatch = currentPatch.younger
    end

    return currentSite
end

# ---------------------------------------------------------------------------
# crown_scorching!  (Fortran crown_scorching)
# ---------------------------------------------------------------------------
"""
    crown_scorching!(currentSite::ed_site_type)

Calculate the scorch height `patch.scorch_ht(pft)` [m] for all PFTs on each
burning patch (Fortran `crown_scorching`). Equation 16 in Thonicke et al. 2010
(Van Wagner 1973 EQ8 / 2/3 Byram 1959).
"""
function crown_scorching!(currentSite::ed_site_type)
    edpft = edpftvarcon_inst()

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        if currentPatch.nocomp_pft_label != nocomp_bareground
            tree_ag_biomass = 0.0
            if currentPatch.fire == 1
                currentCohort = currentPatch.tallest
                while currentCohort !== nothing
                    if prt_params.woody[currentCohort.pft] == itrue  # trees only
                        leaf_c = GetState(currentCohort.prt, leaf_organ, carbon12_element)
                        sapw_c = GetState(currentCohort.prt, sapw_organ, carbon12_element)
                        struct_c = GetState(currentCohort.prt, struct_organ, carbon12_element)

                        tree_ag_biomass = tree_ag_biomass +
                            currentCohort.n * (leaf_c +
                            prt_params.allom_agb_frac[currentCohort.pft] * (sapw_c + struct_c))
                    end
                    currentCohort = currentCohort.shorter
                end

                for i_pft in 1:numpft[]
                    if tree_ag_biomass > 0.0 && prt_params.woody[i_pft] == itrue
                        # Equation 16 in Thonicke et al. 2010
                        currentPatch.scorch_ht[i_pft] =
                            edpft.fire_alpha_SH[i_pft] * (currentPatch.fi^0.667)
                    else
                        currentPatch.scorch_ht[i_pft] = 0.0
                    end
                end
            end
        end
        currentPatch = currentPatch.younger
    end

    return currentSite
end

# ---------------------------------------------------------------------------
# crown_damage!  (Fortran crown_damage)
# ---------------------------------------------------------------------------
"""
    crown_damage!(currentSite::ed_site_type)

Update each tree cohort's `cohort.fraction_crown_burned` [0-1], the proportion of
crown affected by fire (Fortran `crown_damage`). Equation 17 in Thonicke et al.
2010.
"""
function crown_damage!(currentSite::ed_site_type)
    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        if currentPatch.nocomp_pft_label != nocomp_bareground
            if currentPatch.fire == 1
                currentCohort = currentPatch.tallest
                while currentCohort !== nothing
                    currentCohort.fraction_crown_burned = 0.0
                    if prt_params.woody[currentCohort.pft] == itrue  # trees only
                        crown_depth = CrownDepth(currentCohort.height, currentCohort.pft)

                        if currentPatch.scorch_ht[currentCohort.pft] <
                           (currentCohort.height - crown_depth)
                            # Flames lower than bottom of canopy.
                            currentCohort.fraction_crown_burned = 0.0
                        else
                            # Flames part of way up canopy. Eq.17 Thonicke 2010.
                            if (currentCohort.height > 0.0) &&
                               (currentPatch.scorch_ht[currentCohort.pft] >=
                                (currentCohort.height - crown_depth))
                                currentCohort.fraction_crown_burned =
                                    (currentPatch.scorch_ht[currentCohort.pft] -
                                     (currentCohort.height - crown_depth)) / crown_depth
                            else
                                # Flames over top of canopy.
                                currentCohort.fraction_crown_burned = 1.0
                            end
                        end
                        # Check for strange values.
                        currentCohort.fraction_crown_burned =
                            min(1.0, max(0.0, currentCohort.fraction_crown_burned))
                    end
                    currentCohort = currentCohort.shorter
                end
            end
        end
        currentPatch = currentPatch.younger
    end

    return currentSite
end

# ---------------------------------------------------------------------------
# cambial_damage_kill!  (Fortran cambial_damage_kill)
# ---------------------------------------------------------------------------
"""
    cambial_damage_kill!(currentSite::ed_site_type)

Calculate the probability that each tree dies due to cambial char
(`cohort.cambial_mort` [0-1]) using the lethal-stem-heating duration
`patch.tau_l` [min] (Fortran `cambial_damage_kill`). Equations 19-21 in Thonicke
et al. 2010.
"""
function cambial_damage_kill!(currentSite::ed_site_type)
    edpft = edpftvarcon_inst()

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        if currentPatch.nocomp_pft_label != nocomp_bareground
            if currentPatch.fire == 1
                currentCohort = currentPatch.tallest
                while currentCohort !== nothing
                    if prt_params.woody[currentCohort.pft] == itrue  # trees only
                        # Equation 21 in Thonicke et al. 2010 - bark thickness [cm]
                        bt = edpft.bark_scaler[currentCohort.pft] * currentCohort.dbh
                        # Equation 20 - time to kill cambium (min)
                        tau_c = 2.9 * bt^2.0
                        # Equation 19 in Thonicke et al. 2010
                        if (currentPatch.tau_l / tau_c) >= 2.0
                            currentCohort.cambial_mort = 1.0
                        else
                            if (currentPatch.tau_l / tau_c) > 0.22
                                currentCohort.cambial_mort =
                                    (0.563 * (currentPatch.tau_l / tau_c)) - 0.125
                            else
                                currentCohort.cambial_mort = 0.0
                            end
                        end
                    end
                    currentCohort = currentCohort.shorter
                end
            end
        end
        currentPatch = currentPatch.younger
    end

    return currentSite
end

# ---------------------------------------------------------------------------
# post_fire_mortality!  (Fortran post_fire_mortality)
# ---------------------------------------------------------------------------
"""
    post_fire_mortality!(currentSite::ed_site_type)

Calculate each tree cohort's post-fire mortality `cohort.fire_mort` [0-1] as the
joint probability of crown-scorch (`cohort.crownfire_mort`) and cambial-char
(`cohort.cambial_mort`) death, assuming the two are independent (Fortran
`post_fire_mortality`). Equations 18 & 22 in Thonicke et al. 2010.
"""
function post_fire_mortality!(currentSite::ed_site_type)
    edpft = edpftvarcon_inst()

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        if currentPatch.nocomp_pft_label != nocomp_bareground
            if currentPatch.fire == 1
                currentCohort = currentPatch.tallest
                while currentCohort !== nothing
                    currentCohort.fire_mort = 0.0
                    currentCohort.crownfire_mort = 0.0
                    if prt_params.woody[currentCohort.pft] == itrue
                        # Equation 22 in Thonicke et al. 2010
                        currentCohort.crownfire_mort =
                            edpft.crown_kill[currentCohort.pft] *
                            currentCohort.fraction_crown_burned^3.0
                        # Equation 18 in Thonicke et al. 2010 - joint probability
                        currentCohort.fire_mort = max(0.0, min(1.0,
                            currentCohort.crownfire_mort + currentCohort.cambial_mort -
                            (currentCohort.crownfire_mort * currentCohort.cambial_mort)))
                    else
                        # Grass mode of death is removal of leaves.
                        currentCohort.fire_mort = 0.0
                    end
                    currentCohort = currentCohort.shorter
                end
            end
        end
        currentPatch = currentPatch.younger
    end

    return currentSite
end
