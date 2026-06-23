# FatesSoilBGCFluxMod.jl
# Julia port of FATES src/fates/biogeochem/FatesSoilBGCFluxMod.F90
#
# The FATES <-> host soil-biogeochem flux coupling. This module:
#   * UnPackNutrientAquisitionBCs : receives the per-competitor N (NH4/NO3) and P
#     uptake-flux boundary conditions returned by the host soil-BGC model and
#     parses them out to each cohort's daily uptake/demand (or applies a
#     prescribed-uptake fraction of demand when not coupled). Zeroes the input BCs
#     afterward so they re-integrate.
#   * PrepCH4BCs : prepares the output boundary conditions for the host's methane
#     model (above/below-ground NPP, fine-root C, root profile, root respiration,
#     woody-fraction-for-aerenchyma).
#   * PrepNutrientAquisitionBCs : builds the per-competitor fine-root carbon
#     (veg_rootc) and decomposer-microbial-biomass (decompmicc) boundary
#     conditions the host BGC needs, depending on the competition scheme
#     (RD/ECA) and competitor scaling (trivial/coupled).
#   * EffluxIntoLitterPools : transfers per-cohort root exudation/efflux into the
#     patch's root-fines fragmentation (labile) pool, distributed by root depth.
#   * FluxIntoLitterPools : ports the fragmenting-litter fluxes (CWD, leaf/root
#     fines, decaying seeds) into the host's decomposing litter pools, partitioned
#     into cellulose/lignin/labile chemical fractions and distributed over the
#     decomposition layers (surface additions via an exponential depth profile).
#     Also computes the MIMICS ligninC/totalN litter-quality boundary condition.
#
# Translation notes:
#   * fates_r8 -> Float64.
#   * Fortran `bc_in(s)%...`/`bc_out%...` -> the merged bc_in_type / bc_out_type.
#   * Fortran pointer linked lists (oldest_patch/younger, tallest/shorter) ->
#     `while p !== nothing` walks over the merged patch/cohort types.
#   * `ccohort%prt%GetState(...)` -> the free function `GetState(ccohort.prt, ...)`
#     from PRTGenericMod.
#   * Module-global `area`/`area_inv` -> EDTypesMod's `area` (= 10000.0) and
#     `area_inv`. `AREA`/`AREA_INV` in the Fortran are aliases of these.
#   * `n_uptake_mode`/`p_uptake_mode`, `ED_val_cwd_fcel`/`ED_val_cwd_flig` are
#     fields of the EDParams singleton -> accessed via `ed_params()`.
#   * `fates_np_comp_scaling` is a Ref{Int}.
#   * `hlm_nu_com`/`hlm_decomp` are Ref{String}; `hlm_parteh_mode`/`hlm_use_ch4`
#     are Ref{Int}.
#   * `EDPftvarcon_inst`/`prt_params` accessed via `EDPftvarcon_inst[]` / the
#     `prt_params` singleton.
#
# The flux-partitioning + nutrient-demand math is preserved verbatim and is
# mass-conserving.

# =====================================================================================

"""
    UnPackNutrientAquisitionBCs(sites, bc_in)

Receive the nutrient uptake-flux boundary conditions and parse them out to the
cohorts. Called before FATES dynamics (before PARTEH). The uptake fluxes are
assumed to be an integrated daily quantity, accumulated each short-BGC step over
the day; the input BCs are zeroed at the end so they can re-integrate.

`sites` is a vector of `ed_site_type`; `bc_in` is a matching vector of
`bc_in_type` (one per site).
"""
function UnPackNutrientAquisitionBCs(sites::AbstractVector{ed_site_type},
                                     bc_in::AbstractVector{bc_in_type})

    nsites = length(sites)

    # We can exit if this is a c-only simulation
    if hlm_parteh_mode[] == prt_carbon_allom_hyp
        # These can now be zero'd
        for s in 1:nsites
            bc_in[s].plant_nh4_uptake_flux .= 0.0
            bc_in[s].plant_no3_uptake_flux .= 0.0
            bc_in[s].plant_p_uptake_flux   .= 0.0
        end
        return nothing
    end

    edpft = EDPftvarcon_inst[]
    npmode = ed_params().n_uptake_mode
    ppmode = ed_params().p_uptake_mode

    for s in 1:nsites

        # If the plant is in "prescribed uptake mode" then we are not coupling
        # with the soil bgc model. In this case, the bc_in structure is
        # meaningless. Instead, we give the plants a parameterized fraction of
        # their demand.

        if npmode == prescribed_n_uptake

            cpatch = sites[s].oldest_patch
            while cpatch !== nothing
                ccohort = cpatch.tallest
                while ccohort !== nothing
                    pft = ccohort.pft
                    fnrt_c = GetState(ccohort.prt, fnrt_organ, carbon12_element)
                    ccohort.daily_n_demand = fnrt_c *
                        (edpft.vmax_nh4[pft] + edpft.vmax_no3[pft]) * sec_per_day
                    ccohort.daily_nh4_uptake = fnrt_c * edpft.vmax_nh4[pft] *
                        edpft.prescribed_nuptake[pft] * sec_per_day
                    ccohort.daily_no3_uptake = fnrt_c * edpft.vmax_no3[pft] *
                        edpft.prescribed_nuptake[pft] * sec_per_day
                    ccohort = ccohort.shorter
                end
                cpatch = cpatch.younger
            end

        elseif npmode == coupled_n_uptake

            icomp = 0
            cpatch = sites[s].oldest_patch
            while cpatch !== nothing
                ccohort = cpatch.tallest
                while ccohort !== nothing
                    icomp += 1
                    pft = ccohort.pft
                    fnrt_c = GetState(ccohort.prt, fnrt_organ, carbon12_element)
                    ccohort.daily_n_demand = fnrt_c *
                        (edpft.vmax_nh4[pft] + edpft.vmax_no3[pft]) * sec_per_day
                    # N Uptake:  Convert g/m2/day -> kg/plant/day
                    ccohort.daily_nh4_uptake =
                        bc_in[s].plant_nh4_uptake_flux[icomp, 1] * kg_per_g * area / ccohort.n
                    ccohort.daily_no3_uptake =
                        bc_in[s].plant_no3_uptake_flux[icomp, 1] * kg_per_g * area / ccohort.n
                    ccohort = ccohort.shorter
                end
                cpatch = cpatch.younger
            end

        end

        if ppmode == prescribed_p_uptake
            cpatch = sites[s].oldest_patch
            while cpatch !== nothing
                ccohort = cpatch.tallest
                while ccohort !== nothing
                    pft = ccohort.pft
                    fnrt_c = GetState(ccohort.prt, fnrt_organ, carbon12_element)
                    ccohort.daily_p_demand = fnrt_c * edpft.vmax_p[pft] * sec_per_day
                    ccohort.daily_p_gain = fnrt_c * edpft.vmax_p[pft] * sec_per_day *
                        edpft.prescribed_nuptake[pft]
                    ccohort = ccohort.shorter
                end
                cpatch = cpatch.younger
            end

        elseif ppmode == coupled_p_uptake

            icomp = 0
            cpatch = sites[s].oldest_patch
            while cpatch !== nothing
                ccohort = cpatch.tallest
                while ccohort !== nothing
                    icomp += 1
                    pft = ccohort.pft
                    fnrt_c = GetState(ccohort.prt, fnrt_organ, carbon12_element)
                    ccohort.daily_p_demand = fnrt_c * edpft.vmax_p[pft] * sec_per_day
                    # P Uptake:  Convert g/m2/day -> kg/plant/day
                    ccohort.daily_p_gain =
                        bc_in[s].plant_p_uptake_flux[icomp, 1] * kg_per_g * area / ccohort.n
                    ccohort = ccohort.shorter
                end
                cpatch = cpatch.younger
            end

        end

        # These can now be zero'd
        bc_in[s].plant_nh4_uptake_flux .= 0.0
        bc_in[s].plant_no3_uptake_flux .= 0.0
        bc_in[s].plant_p_uptake_flux   .= 0.0

    end

    return nothing
end

# =====================================================================================

"""
    PrepCH4BCs(csite, bc_in, bc_out)

Prepare the output boundary conditions for the host's methane (CH4) calculations:
per-patch above/below-ground NPP, fine-root C, the depth-normalized root profile,
root respiration over depth, and the woody crown-area fraction (for aerenchyma).
"""
function PrepCH4BCs(csite::ed_site_type, bc_in::bc_in_type, bc_out::bc_out_type)

    ema_npp_tscale = 10.0  # 10 day (kept for parity; unused in current form)

    # Initialize to zero
    bc_out.annavg_agnpp_pa    .= 0.0
    bc_out.annavg_bgnpp_pa    .= 0.0
    bc_out.annsum_npp_pa      .= 0.0
    bc_out.rootfr_pa          .= 0.0
    bc_out.frootc_pa          .= 0.0
    bc_out.root_resp          .= 0.0
    bc_out.woody_frac_aere_pa .= 0.0
    bc_out.ema_npp = 0.0

    p = prt_params

    fp = 0
    cpatch = csite.oldest_patch
    while cpatch !== nothing
        if cpatch.nocomp_pft_label != nocomp_bareground
            # Patch ordering when passing boundary conditions always goes from
            # oldest to youngest, following EDPatchDynamics::set_patchno().
            fp += 1

            agnpp = 0.0
            bgnpp = 0.0
            woody_area = 0.0
            plant_area = 0.0

            ccohort = cpatch.tallest
            while ccohort !== nothing

                # For consistency, only apply calculations to non-new cohorts.
                # New cohorts will not have respiration rates at this point.
                if !ccohort.isnew

                    pft = ccohort.pft

                    set_root_fraction(csite.rootfrac_scr, pft, csite.zi_soil;
                                      max_nlevroot=bc_in.max_rooting_depth_index_col)

                    fnrt_c = GetState(ccohort.prt, fnrt_organ, carbon12_element)

                    # [kgC/day]
                    sapw_net_alloc   = GetNetAlloc(ccohort.prt, sapw_organ, carbon12_element) * days_per_sec
                    store_net_alloc  = GetNetAlloc(ccohort.prt, store_organ, carbon12_element) * days_per_sec
                    leaf_net_alloc   = GetNetAlloc(ccohort.prt, leaf_organ, carbon12_element) * days_per_sec
                    fnrt_net_alloc   = GetNetAlloc(ccohort.prt, fnrt_organ, carbon12_element) * days_per_sec
                    struct_net_alloc = GetNetAlloc(ccohort.prt, struct_organ, carbon12_element) * days_per_sec
                    repro_net_alloc  = GetNetAlloc(ccohort.prt, repro_organ, carbon12_element) * days_per_sec

                    # [kgC/plant/day] -> [gC/m2/s]
                    agnpp = agnpp + ccohort.n / cpatch.area * (leaf_net_alloc + repro_net_alloc +
                        p.allom_agb_frac[pft] * (sapw_net_alloc + store_net_alloc + struct_net_alloc)) * g_per_kg

                    # [kgC/plant/day] -> [gC/m2/s]
                    bgnpp = bgnpp + ccohort.n / cpatch.area * (fnrt_net_alloc +
                        (1.0 - p.allom_agb_frac[pft]) * (sapw_net_alloc + store_net_alloc + struct_net_alloc)) * g_per_kg

                    if hlm_use_ch4[] == itrue

                        nlev = bc_in.nlevsoil

                        # Fine root fraction over depth
                        for j in 1:nlev
                            bc_out.rootfr_pa[fp, j] += csite.rootfrac_scr[j]
                        end

                        # Fine root carbon, convert [kg/plant] -> [g/m2]
                        bc_out.frootc_pa[fp] += fnrt_c * ccohort.n / cpatch.area * g_per_kg

                        # (gC/m2/s) root respiration (fine root MR + total root GR)
                        # (kgC/indiv/yr) -> gC/m2/s
                        for j in 1:nlev
                            bc_out.root_resp[j] += ccohort.resp_acc_hold * years_per_day * g_per_kg *
                                days_per_sec * ccohort.n * area_inv *
                                (1.0 - p.allom_agb_frac[pft]) * csite.rootfrac_scr[j]
                        end

                    end

                    if p.woody[pft] == itrue
                        woody_area = woody_area + ccohort.c_area
                    end
                    plant_area = plant_area + ccohort.c_area

                end

                ccohort = ccohort.shorter
            end

            if hlm_use_ch4[] == itrue
                nlev = bc_in.nlevsoil
                rootfr_sum = 0.0
                for j in 1:nlev
                    rootfr_sum += bc_out.rootfr_pa[fp, j]
                end
                if rootfr_sum > nearzero
                    for j in 1:nlev
                        bc_out.rootfr_pa[fp, j] = bc_out.rootfr_pa[fp, j] / rootfr_sum
                    end
                end

                # gC/m2/s
                bc_out.annavg_agnpp_pa[fp] = agnpp
                bc_out.annavg_bgnpp_pa[fp] = bgnpp
                # gc/m2/yr
                bc_out.annsum_npp_pa[fp] = (bgnpp + agnpp) * days_per_year * sec_per_day

                if plant_area > nearzero
                    bc_out.woody_frac_aere_pa[fp] = woody_area / plant_area
                end
            end
        end
        cpatch = cpatch.younger
    end

    bc_out.ema_npp = csite.ema_npp

    return nothing
end

# =====================================================================================

"""
    PrepNutrientAquisitionBCs(csite, bc_in, bc_out)

Generate the boundary-condition output structures for the host soil-BGC nutrient
competition, depending on (1) the active competition method (`hlm_nu_com`,
RD/ECA) and (2) the competitor scaling type (`fates_np_comp_scaling`,
trivial/coupled). Fills per-competitor fine-root carbon `veg_rootc`, the
competitor PFT index map `ft_index`, the number of plant competitors
`num_plant_comps`, and (ECA only) the decomposer microbial biomass `decompmicc`.
"""
function PrepNutrientAquisitionBCs(csite::ed_site_type, bc_in::bc_in_type,
                                   bc_out::bc_out_type)

    decompmicc_lambda = 2.5     # Depth attenuation exponent for decomposer biomass
    decompmicc_zmax   = 7.0e-2  # Depth of maximum decomposer biomass

    edpft = EDPftvarcon_inst[]

    # Whether this is a trivial or coupled run, the following variables get
    # initialized in the same way
    bc_out.veg_rootc .= 0.0
    bc_out.ft_index  .= -1
    if hlm_nu_com[] == "ECA"
        bc_out.decompmicc .= 0.0
        bc_out.cn_scalar  .= 1.0
        bc_out.cp_scalar  .= 1.0
    end

    if fates_np_comp_scaling[] == trivial_np_comp_scaling
        if hlm_nu_com[] == "RD"
            bc_out.num_plant_comps = 1
            bc_out.ft_index[1] = 1
            return nothing
        end
    end

    # For both the trivial case with ECA, and the coupled case we still need to
    # calculate the root biomass and decompmicc arrays (the former for the latter
    # when trivial). So we don't differentiate.

    icomp = 0
    cpatch = csite.oldest_patch
    while cpatch !== nothing
        ccohort = cpatch.tallest
        while ccohort !== nothing

            if fates_np_comp_scaling[] == coupled_np_comp_scaling
                icomp += 1
            else
                icomp = 1
            end

            pft = ccohort.pft
            bc_out.ft_index[icomp] = pft

            set_root_fraction(csite.rootfrac_scr, pft, csite.zi_soil;
                              max_nlevroot=bc_in.max_rooting_depth_index_col)

            fnrt_c = GetState(ccohort.prt, fnrt_organ, carbon12_element)

            # Map the soil layers to the decomposition layers (which may be synonomous)
            # veg_rootc in units: [gC/m3] = [kgC/plant] * [plant/ha] * [ha/10k m2] * [1000 g/kg] * [1/m]
            for j in 1:bc_in.nlevdecomp
                id = bc_in.decomp_id[j]  # Map from soil layer to decomp layer
                veg_rootc = fnrt_c * ccohort.n * csite.rootfrac_scr[j] * area_inv *
                    g_per_kg / csite.dz_soil[j]

                bc_out.veg_rootc[icomp, id] += veg_rootc

                if hlm_nu_com[] == "ECA"
                    # 2-parameter exponential attenuation function for decomposer biomass
                    decompmicc_layer = edpft.decompmicc[pft] *
                        exp(-decompmicc_lambda * abs(csite.z_soil[j] - decompmicc_zmax))

                    bc_out.decompmicc[id] += decompmicc_layer * veg_rootc
                end
            end
            ccohort = ccohort.shorter
        end
        cpatch = cpatch.younger
    end

    # Calculate decomposer microbial biomass by weighting with root biomass.
    # This is just the normalization step.
    if hlm_nu_com[] == "ECA"
        for id in 1:bc_in.nlevdecomp
            colsum = 0.0
            for ic in 1:size(bc_out.veg_rootc, 1)
                colsum += bc_out.veg_rootc[ic, id]
            end
            bc_out.decompmicc[id] = bc_out.decompmicc[id] / max(nearzero, colsum)
        end
    end

    if fates_np_comp_scaling[] == coupled_np_comp_scaling
        bc_out.num_plant_comps = icomp
    else
        bc_out.num_plant_comps = 1
    end

    return nothing
end

# =====================================================================================

"""
    EffluxIntoLitterPools(csite, cpatch, ccohort, bc_in)

Transfer one cohort's root exudation/efflux (C/N/P) into the patch's root-fines
fragmentation (labile) pool, distributed over depth by the cohort's root profile.
Reuses `root_fines_frag` to save memory and reuse mass-balance/disturbance
machinery. Does NOT increment the site-level `frag_out` mass-flux check (that is
incremented later in the call sequence).
"""
function EffluxIntoLitterPools(csite::ed_site_type, cpatch::fates_patch_type,
                               ccohort::fates_cohort_type, bc_in::bc_in_type)

    set_root_fraction(csite.rootfrac_scr, ccohort.pft, csite.zi_soil;
                      max_nlevroot=bc_in.max_rooting_depth_index_col)

    # Loop over the different elements.
    for el in 1:num_elements[]

        elid = element_list[el]
        efflux = if elid == carbon12_element
            ccohort.daily_c_efflux
        elseif elid == nitrogen_element
            ccohort.daily_n_efflux
        elseif elid == phosphorus_element
            ccohort.daily_p_efflux
        else
            0.0
        end

        litt = cpatch.litter[el]

        for j in 1:csite.nlevsoil
            # kg/m2/day
            litt.root_fines_frag[ilabile, j] += efflux * ccohort.n * area_inv *
                csite.rootfrac_scr[j]
        end

    end

    return nothing
end

# =====================================================================================

"""
    FluxIntoLitterPools(csite, bc_in, bc_out)

Take the flux out of the fragmenting litter pools and port it into the
decomposing litter pools, for each PARTEH element. Fragmenting pools physically
fragment but do not respire; decomposing (host BGC) pools are humus/non-flammable.

The fragmenting CWD, leaf/root fines, and decaying seeds are partitioned into
cellulose/lignin/labile chemical fractions and distributed over decomposition
layers (above-ground/surface additions via an exponential depth profile,
below-ground via the soil-to-decomp layer map), then normalized over layer depth
and converted kg/m2/day -> g/m3/s. Also incorporates the
SoilBiogeochemVerticalProfile functionality.

When coupled with MIMICS, computes the litter-quality boundary condition
`litt_flux_ligc_per_n` (ligninC per total N): from the N flux directly if N is
tracked, else estimated from per-organ stoichiometry as a proxy.

Created by Charlie Koven and Rosie Fisher, 2014-2015.
"""
function FluxIntoLitterPools(csite::ed_site_type, bc_in::bc_in_type,
                             bc_out::bc_out_type)

    surfprof_exp = 10.0  # 1/e-folding depth (1/m) for surface components

    edpft = EDPftvarcon_inst[]
    p = prt_params
    ED_val_cwd_fcel = ed_params().ED_val_cwd_fcel
    ED_val_cwd_flig = ed_params().ED_val_cwd_flig

    # Number of effective soil layers to transfer from
    nlev_eff_soil = max(bc_in.max_rooting_depth_index_col, 1)

    # The decomposition layers are most likely the exact same layers as the soil
    # layers (same depths also), unless single-layer (nlevdecomp = 1).
    nlev_eff_decomp = min(bc_in.nlevdecomp, nlev_eff_soil)

    # Define a single shallow surface profile for surface additions (leaves,
    # stems, N deposition). Sends above-ground mass into the soil pools using an
    # exponential depth decay function. Since it sends an absolute mass [kg] into
    # variable layer widths, multiply the profile by the layer width so wider
    # layers get proportionally more. After the masses are sent, each layer will
    # normalize by depth.
    surface_prof = zeros(Float64, bc_in.nlevsoil)
    z_decomp = 0.0
    for id in 1:nlev_eff_decomp
        z_decomp = z_decomp + 0.5 * bc_in.dz_decomp_sisl[id]
        surface_prof[id] = exp(-surfprof_exp * z_decomp) * bc_in.dz_decomp_sisl[id]
        z_decomp = z_decomp + 0.5 * bc_in.dz_decomp_sisl[id]
    end
    surface_prof_tot = sum(@view surface_prof[1:nlev_eff_decomp])
    for id in 1:nlev_eff_decomp
        surface_prof[id] = surface_prof[id] / surface_prof_tot
    end

    # Loop over the different elements.
    for el in 1:num_elements[]

        elid = element_list[el]

        # Make handles to the cellulose/labile/lignin flux partitions, and zero
        # the boundary flux arrays.
        local flux_cel_si::Vector{Float64}
        local flux_lab_si::Vector{Float64}
        local flux_lig_si::Vector{Float64}
        if elid == carbon12_element
            bc_out.litt_flux_cel_c_si .= 0.0
            bc_out.litt_flux_lig_c_si .= 0.0
            bc_out.litt_flux_lab_c_si .= 0.0
            flux_cel_si = bc_out.litt_flux_cel_c_si
            flux_lab_si = bc_out.litt_flux_lab_c_si
            flux_lig_si = bc_out.litt_flux_lig_c_si
        elseif elid == nitrogen_element
            bc_out.litt_flux_cel_n_si .= 0.0
            bc_out.litt_flux_lig_n_si .= 0.0
            bc_out.litt_flux_lab_n_si .= 0.0
            flux_cel_si = bc_out.litt_flux_cel_n_si
            flux_lab_si = bc_out.litt_flux_lab_n_si
            flux_lig_si = bc_out.litt_flux_lig_n_si
        elseif elid == phosphorus_element
            bc_out.litt_flux_cel_p_si .= 0.0
            bc_out.litt_flux_lig_p_si .= 0.0
            bc_out.litt_flux_lab_p_si .= 0.0
            flux_cel_si = bc_out.litt_flux_cel_p_si
            flux_lab_si = bc_out.litt_flux_lab_p_si
            flux_lig_si = bc_out.litt_flux_lig_p_si
        else
            continue
        end

        currentPatch = csite.oldest_patch
        while currentPatch !== nothing

            # Pointer to the litter object for the current element on the patch
            litt = currentPatch.litter[el]
            area_frac = currentPatch.area / area

            for ic in 1:ncwd

                for id in 1:nlev_eff_decomp
                    flux_cel_si[id] += litt.ag_cwd_frag[ic] * ED_val_cwd_fcel * area_frac * surface_prof[id]
                    flux_lig_si[id] += litt.ag_cwd_frag[ic] * ED_val_cwd_flig * area_frac * surface_prof[id]
                end

                for j in 1:nlev_eff_soil
                    id = bc_in.decomp_id[j]  # Map from soil layer to decomp layer
                    flux_cel_si[id] += litt.bg_cwd_frag[ic, j] * ED_val_cwd_fcel * area_frac
                    flux_lig_si[id] += litt.bg_cwd_frag[ic, j] * ED_val_cwd_flig * area_frac
                end
            end

            # leaf fragmentation fluxes
            for id in 1:nlev_eff_decomp
                flux_lab_si[id] += litt.leaf_fines_frag[ilabile]    * area_frac * surface_prof[id]
                flux_cel_si[id] += litt.leaf_fines_frag[icellulose] * area_frac * surface_prof[id]
                flux_lig_si[id] += litt.leaf_fines_frag[ilignin]    * area_frac * surface_prof[id]
            end

            # decaying seeds from the litter pool
            for ipft in 1:numpft[]
                for id in 1:nlev_eff_decomp
                    seedflux = litt.seed_decay[ipft] + litt.seed_germ_decay[ipft]
                    flux_lab_si[id] += seedflux * edpft.lf_flab[ipft] * area_frac * surface_prof[id]
                    flux_cel_si[id] += seedflux * edpft.lf_fcel[ipft] * area_frac * surface_prof[id]
                    flux_lig_si[id] += seedflux * edpft.lf_flig[ipft] * area_frac * surface_prof[id]
                end
            end

            # root fragmentation fluxes
            for j in 1:nlev_eff_soil
                id = bc_in.decomp_id[j]
                flux_lab_si[id] += litt.root_fines_frag[ilabile, j]    * area_frac
                flux_cel_si[id] += litt.root_fines_frag[icellulose, j] * area_frac
                flux_lig_si[id] += litt.root_fines_frag[ilignin, j]    * area_frac
            end

            currentPatch = currentPatch.younger
        end

        # Normalize all masses over the decomposition layer's depth.
        # Convert from kg/m2/day -> g/m3/s.
        for id in 1:nlev_eff_decomp
            flux_cel_si[id] = days_per_sec * g_per_kg * flux_cel_si[id] / bc_in.dz_decomp_sisl[id]
            flux_lig_si[id] = days_per_sec * g_per_kg * flux_lig_si[id] / bc_in.dz_decomp_sisl[id]
            flux_lab_si[id] = days_per_sec * g_per_kg * flux_lab_si[id] / bc_in.dz_decomp_sisl[id]
        end

    end  # do elements

    # If coupled with MIMICS, we need an assessment of litter quality, ie
    # ligC/totalN. If we are not tracking N in the litter flux (ie C-only model)
    # we approximate by estimating mean C:N ratios of each plant organ and
    # multiplying by the C fluxes to get an approximate N flux. Note: in C-only,
    # we will not capture any re-absorption.
    if hlm_decomp[] == "MIMICS"

        local sum_N::Float64

        # If we track nitrogen (cnp or other), diagnose the c-lig/n ratio
        # directly from the pools.
        if element_pos[nitrogen_element] > 0

            # Sum totalN fluxes over depth [g/m2]
            sum_N = 0.0
            for j in 1:nlev_eff_soil
                sum_N += (bc_out.litt_flux_cel_n_si[j] +
                          bc_out.litt_flux_lig_n_si[j] +
                          bc_out.litt_flux_lab_n_si[j]) * bc_in.dz_sisl[j]
            end

        else

            # Carbon-Only: use the stoichiometry parameters to estimate the C:N
            # of live vegetation and the seedbank, and use that as a proxy for
            # the C:N of the litter flux.
            sum_N = 0.0

            currentPatch = csite.oldest_patch
            while currentPatch !== nothing

                litt = currentPatch.litter[element_pos[carbon12_element]]
                area_frac = currentPatch.area * area_inv

                tot_leaf_c = 0.0
                tot_leaf_n = 0.0
                tot_fnrt_c = 0.0
                tot_fnrt_n = 0.0
                tot_wood_c = 0.0
                tot_wood_n = 0.0

                ccohort = currentPatch.tallest
                while ccohort !== nothing
                    ipft = ccohort.pft
                    leaf_c   = ccohort.n * area_inv * GetState(ccohort.prt, leaf_organ, carbon12_element)
                    sapw_c   = ccohort.n * area_inv * GetState(ccohort.prt, sapw_organ, carbon12_element)
                    fnrt_c   = ccohort.n * area_inv * GetState(ccohort.prt, fnrt_organ, carbon12_element)
                    struct_c = ccohort.n * area_inv * GetState(ccohort.prt, struct_organ, carbon12_element)
                    leaf_n   = leaf_c   * p.nitr_stoich_p1[ipft, p.organ_param_id[leaf_organ]]
                    sapw_n   = sapw_c   * p.nitr_stoich_p1[ipft, p.organ_param_id[sapw_organ]]
                    fnrt_n   = fnrt_c   * p.nitr_stoich_p1[ipft, p.organ_param_id[fnrt_organ]]
                    struct_n = struct_c * p.nitr_stoich_p1[ipft, p.organ_param_id[struct_organ]]
                    tot_leaf_c += leaf_c
                    tot_leaf_n += leaf_n
                    tot_fnrt_c += fnrt_c
                    tot_fnrt_n += fnrt_n
                    tot_wood_c += sapw_c + struct_c
                    tot_wood_n += sapw_n + struct_n
                    ccohort = ccohort.shorter
                end

                if tot_wood_c > nearzero
                    sum_N += area_frac * sum(litt.ag_cwd_frag) * (tot_wood_n / tot_wood_c)
                    sum_N += area_frac * sum(litt.bg_cwd_frag) * (tot_wood_n / tot_wood_c)
                end
                if tot_leaf_c > nearzero
                    sum_N += area_frac * sum(litt.leaf_fines_frag) * (tot_leaf_n / tot_leaf_c)
                end
                if tot_fnrt_c > nearzero
                    sum_N += area_frac * sum(litt.root_fines_frag) * (tot_fnrt_n / tot_fnrt_c)
                end
                for ipft in 1:numpft[]
                    sum_N += area_frac * currentPatch.nitr_repro_stoich[ipft] *
                        (litt.seed_decay[ipft] + litt.seed_germ_decay[ipft])
                end

                currentPatch = currentPatch.younger
            end

            # Convert from kg/m2/day -> g/m2/s
            sum_N = sum_N * days_per_sec * g_per_kg

        end

        # Sum over layers and multiply by depth g/m3/s * m -> g/m2/s
        sum_ligC = 0.0
        for j in 1:nlev_eff_soil
            sum_ligC += bc_out.litt_flux_lig_c_si[j] * bc_in.dz_sisl[j]
        end

        if sum_N > nearzero
            bc_out.litt_flux_ligc_per_n = sum_ligC / sum_N
        else
            bc_out.litt_flux_ligc_per_n = 0.0
        end

    end

    return nothing
end
