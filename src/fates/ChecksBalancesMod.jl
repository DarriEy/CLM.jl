# ChecksBalancesMod.jl
# Julia port of FATES src/fates/main/ChecksBalancesMod.F90
#
# Site-level carbon / mass-balance checking + diagnostics. Sums the site carbon
# (or N / P) stocks and fluxes and verifies conservation:
#   * SiteMassStock           — sum biomass + seed + litter stock over all
#     patches of a site for a given element.
#   * PatchMassStock          — sum the live-plant, seed, and litter stocks for a
#     single patch and element.
#   * CheckLitterPools        — sanity check that no litter pool went negative
#     (a debugging aid, not a carbon-balance check).
#   * CheckIntegratedMassPools — accumulate the time-integrated net fluxes in /
#     out of vegetation and litter and compare against the instantaneous state
#     (writes the err_liveveg / err_litter diagnostics).
#
# Translation notes:
#   * fates_r8 -> Float64. The Fortran `pointer` patch/cohort linked lists become
#     the Union{...,Nothing} walks already used by the merged FATES site/patch/
#     cohort types (oldest_patch -> younger ; tallest -> shorter).
#   * Operates on the merged ed_site_type / fates_patch_type / fates_cohort_type
#     + the EDTypes helper types (site_ifluxbal_type / site_fluxdiags_type /
#     elem_diag_type / site_massbal_type) + the PRTGenericMod element/organ
#     consts (carbon12/nitrogen/phosphorus_element, leaf/fnrt/sapw/store/repro/
#     struct_organ, num_elements, element_list) + FatesLitterMod (litter_type,
#     ncwd, ndcmpy).
#   * `currentPatch%litter(el)` -> `currentPatch.litter[el]`; `area`/`area_inv`
#     are the EDTypesMod module consts.
#   * The two endrun balance-comparison blocks are commented out in the Fortran
#     (the tolerance check is disabled there) — preserved as comments here.
#   * The stock-summation, flux-accumulation, and balance math are byte-for-byte
#     faithful to the Fortran.

"""
    SiteMassStock(currentSite, el) -> (total_stock, biomass_stock, litter_stock, seed_stock)

Sum the total carbon (or N / P) stock on a site for element index `el` (the FATES
element position, not the global PARTEH id), by walking the age-ordered patch
linked list (oldest_patch -> younger) and accumulating each patch's
[`PatchMassStock`](@ref). All stocks are in [kg].

Returns the tuple `(total_stock, biomass_stock, litter_stock, seed_stock)` where
`total_stock = biomass_stock + seed_stock + litter_stock`.
"""
function SiteMassStock(currentSite::ed_site_type, el::Integer)

    litter_stock  = 0.0
    biomass_stock = 0.0
    seed_stock    = 0.0

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing

        patch_biomass, patch_seed, patch_litter = PatchMassStock(currentPatch, el)
        litter_stock  = litter_stock  + patch_litter
        biomass_stock = biomass_stock + patch_biomass
        seed_stock    = seed_stock    + patch_seed

        currentPatch = currentPatch.younger
    end  # end patch loop

    total_stock = biomass_stock + seed_stock + litter_stock

    return total_stock, biomass_stock, litter_stock, seed_stock
end

# =====================================================================================

"""
    PatchMassStock(currentPatch, el) -> (live_stock, seed_stock, litter_stock)

Sum the mass of the different stocks on a single patch for element index `el`:
non-seed litter (CWD + fines), viable + germinated seed, and the live-plant
mass (the six organs summed over all cohorts, weighted by cohort number density
`n`). All stocks in [kg].
"""
function PatchMassStock(currentPatch::fates_patch_type, el::Integer)

    litt       = currentPatch.litter[el]
    element_id = element_list[el]

    # Total non-seed litter in [kg]
    litter_stock = currentPatch.area *
        (sum(litt.ag_cwd)     +
         sum(litt.bg_cwd)     +
         sum(litt.leaf_fines) +
         sum(litt.root_fines))

    # Total mass of viable seeds in [kg]
    seed_stock = currentPatch.area *
        (sum(litt.seed) + sum(litt.seed_germ))

    # Total mass on living plants
    live_stock = 0.0
    currentCohort = currentPatch.tallest
    while currentCohort !== nothing
        live_stock = live_stock +
            (GetState(currentCohort.prt, struct_organ, element_id) +
             GetState(currentCohort.prt, sapw_organ,   element_id) +
             GetState(currentCohort.prt, leaf_organ,   element_id) +
             GetState(currentCohort.prt, fnrt_organ,   element_id) +
             GetState(currentCohort.prt, store_organ,  element_id) +
             GetState(currentCohort.prt, repro_organ,  element_id)) *
            currentCohort.n
        currentCohort = currentCohort.shorter
    end  # end cohort loop

    return live_stock, seed_stock, litter_stock
end

# =====================================================================================

"""
    CheckLitterPools(currentSite, bc_in)

Check that the litter pools do not have weird (negative) values. This is a
debugging aid, NOT a carbon-balance check; it errors (via `error`, mirroring the
Fortran `endrun`) on the first negative pool it finds, reporting the pool index,
soil layer / PFT / decomposability index, element id, value, and lat/lon.

Mirrors the Fortran `do_litter_checks` flag (default `true`); set it to `false`
to skip the scan entirely.
"""
function CheckLitterPools(currentSite::ed_site_type, bc_in::bc_in_type)

    # We only really run this scheme if we think things are really broken.
    # The balance checks should be our first line of defense that are always on.
    do_litter_checks = true

    do_litter_checks || return nothing

    numlevsoil = bc_in.nlevsoil

    for el in 1:num_elements[]

        currentPatch = currentSite.oldest_patch
        while currentPatch !== nothing

            litt       = currentPatch.litter[el]
            element_id = litt.element_id

            for c in 1:ncwd
                if litt.ag_cwd[c] < 0.0
                    @info "In pool: $c"
                    @info "Element id: $element_id"
                    @info "Negative AG CWD: $(litt.ag_cwd[c])"
                    @info "lat/lon: $(currentSite.lat) $(currentSite.lon)"
                    error("ChecksBalancesMod: negative AG CWD pool")
                end

                for ilyr in 1:numlevsoil
                    if litt.bg_cwd[c, ilyr] < 0.0
                        @info "In pool: $c"
                        @info "Soil layer: $ilyr"
                        @info "Element id: $element_id"
                        @info "Negative BG CWD: $(litt.bg_cwd[c, ilyr])"
                        @info "lat/lon: $(currentSite.lat) $(currentSite.lon)"
                        error("ChecksBalancesMod: negative BG CWD pool")
                    end
                end
            end

            for pft in 1:numpft[]
                if litt.seed[pft] < 0.0
                    @info "For PFT: $pft"
                    @info "Element id: $element_id"
                    @info "Negative seed pool: $(litt.seed[pft])"
                    @info "lat/lon: $(currentSite.lat) $(currentSite.lon)"
                    error("ChecksBalancesMod: negative seed pool")
                end
            end

            for dcmpy in 1:ndcmpy
                if litt.leaf_fines[dcmpy] < 0.0
                    @info "Element id: $element_id"
                    @info "Negative leaf fine litter: $(litt.leaf_fines[dcmpy])"
                    @info "lat/lon: $(currentSite.lat) $(currentSite.lon)"
                    error("ChecksBalancesMod: negative leaf fine litter")
                end

                for ilyr in 1:numlevsoil
                    if litt.root_fines[dcmpy, ilyr] < 0.0
                        @info "For PFT: $dcmpy"
                        @info "Soil layer: $ilyr"
                        @info "Element id: $element_id"
                        @info "Negative root fine litter: $(litt.root_fines[dcmpy, ilyr])"
                        @info "lat/lon: $(currentSite.lat) $(currentSite.lon)"
                        error("ChecksBalancesMod: negative root fine litter")
                    end
                end
            end

            currentPatch = currentPatch.older
        end
    end

    return nothing
end

# ==================================================================

"""
    CheckIntegratedMassPools(site)

Check that the time-integrated net fluxes in / out of the vegetation and the
litter match the quantity in those pools. For each element, the instantaneous
stock is assessed via [`SiteMassStock`](@ref), the integrated daily fluxes are
accumulated into `site.iflux_balance[el]`, and the discrepancy is written to the
flux diagnostics (`err_liveveg` / `err_litter`).

Fluxes are in [kg/m2/day] applied at /day frequency, so they integrate to
[kg/m2] without a unit conversion. `area_inv` converts the per-site [kg/site]
stocks to per-area [kg/m2]. The Fortran's tolerance-comparison endrun blocks are
disabled (commented out) there and remain so here.

Mirrors the Fortran `check_iflux_bal` flag (default `true`).
"""
function CheckIntegratedMassPools(site::ed_site_type)

    check_iflux_bal = true

    if check_iflux_bal

        # For carbon balance checks, we need to initialize the total carbon stock
        for el in 1:num_elements[]
            total_stock, biomass_stock, litter_stock, seed_stock =
                SiteMassStock(site, el)

            ibal      = site.iflux_balance[el]
            ediag     = site.flux_diags.elem[el]
            diag      = site.flux_diags
            site_mass = site.mass_balance[el]

            # Initialize the integrated flux balance diagnostics
            # No need to initialize the instantaneous states, those are re-calculated
            ibal.state_liveveg = (biomass_stock + seed_stock) * area_inv
            ibal.state_litter  = litter_stock * area_inv

            # Flux for live veg: net uptake (either NPP or net nutrient uptake) +
            #                    net spatial seed flux -
            #                    veg turnover to litter -
            #                    veg loss to fire (SEEDS DONT BURN) -
            #                    veg loss from exported harvest -
            #                    seed turnover
            tot_litter_input = (sum(ediag.surf_fine_litter_input) +
                                sum(ediag.root_litter_input) +
                                sum(ediag.cwd_ag_input) +
                                sum(ediag.cwd_bg_input)) * area_inv

            element = element_list[el]
            if element == carbon12_element
                net_uptake = diag.npp + site_mass.net_root_uptake * area_inv
            elseif element == nitrogen_element
                net_uptake = site_mass.net_root_uptake * area_inv
            elseif element == phosphorus_element
                net_uptake = site_mass.net_root_uptake * area_inv
            else
                @info "FATES: an invalid chemical species was detected"
                error("ChecksBalancesMod: invalid chemical species")
            end

            # Fluxes are in [kg/m2/day] so they can be added, the frequency is
            # /day so no conversion necessary to integrate to [kg/m2]
            ibal.iflux_liveveg = ibal.iflux_liveveg +
                (net_uptake
                 - tot_litter_input
                 - ediag.burned_liveveg
                 - ediag.exported_harvest
                 + site_mass.seed_in  * area_inv
                 - site_mass.seed_out * area_inv
                 - ediag.tot_seed_turnover)

            # Flux for litter: veg turnover + seed turnover - "these are tot_litter_input"
            #                  fragmentation -
            #                  burned litter
            ibal.iflux_litter = ibal.iflux_litter +
                tot_litter_input -
                (site_mass.frag_out * area_inv - ediag.tot_seed_turnover) -
                (site_mass.burn_flux_to_atm * area_inv - ediag.burned_liveveg)

            ediag.err_liveveg = ibal.iflux_liveveg - ibal.state_liveveg
            ediag.err_litter  = ibal.iflux_litter - ibal.state_litter

            # Perform the comparison between integrated flux and state
            # (DISABLED in the Fortran — tolerance endrun blocks commented out)
            #
            # if abs(ediag.err_liveveg) > iflux_tol[el]
            #     ... endrun ...
            # end
            # if abs(ibal.state_litter - ibal.iflux_litter) > iflux_tol[el]
            #     ... endrun ...
            # end
        end
    end

    return nothing
end
