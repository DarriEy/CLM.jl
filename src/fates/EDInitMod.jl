# EDInitMod.jl
# FATES (Tier F) Batch 17 — FATES cold-start state initialization.
#
# Faithful port of CTSM FATES `main/EDInitMod.F90` (1341 lines).  This is the
# DEFAULT cold-start path that builds the FATES demographic state at startup
# (the counterpart to FatesInventoryInitMod, Batch 16, which builds it from
# external forest-inventory files):
#
#   * `init_site_vars!`     — allocate the site-level (size x pft / damage / soil /
#                             flux-diagnostic / seed) arrays from the boundary
#                             conditions and attach the fire-weather object.
#   * `zero_site!`          — zero/unset every site-level accumulator and state.
#   * `set_site_properties!`— set the initial site properties on a cold start
#                             (phenology dates/status, GDD, water/temp memory, fire
#                             ignitions, the fixed-biogeog area_pft mapping +
#                             nocomp PFT trimming + renormalization, and the
#                             min_allowed_landuse_fraction).
#   * `init_patches!`       — create the starting patch(es): a bareground patch
#                             (nocomp+fixed-biogeog) plus one or more vegetated
#                             patches per land-use state / nocomp PFT, area-account
#                             them to the site `area`, seed each with cohorts,
#                             initialize the mass-balance stock, and renumber.
#   * `init_cohorts!`       — seed a patch with the initial cohorts (one small
#                             cohort of each "used" PFT) at the prescribed initial
#                             density (or dbh, in nocomp), sized by allometry +
#                             PARTEH, then fuse/sort.
#
# Reuses the already-ported constructors (NOT reimplemented):
#   * patch `Create!` (FatesPatchMod), litter `InitConditions!`,
#   * `create_cohort` / `InitPRTObject!` / `sort_cohorts` / `fuse_cohorts`
#     (EDCohortDynamicsMod, Batch 12),
#   * `set_patchno` (EDPatchDynamicsMod, Batch 14),
#   * `SiteMassStock` (ChecksBalancesMod), `calculate_SP_properties` (EDPhysiologyMod),
#   * allometry h2d/h/bagw/bbgw/bleaf/bfineroot/bsap/bdead/bstore/carea
#     (FatesAllometryMod), `SetState!` / `CheckInitialConditions` (PRTGenericMod),
#   * `StorageNutrientTarget`, `get_age_class_index`, the site zero-helpers
#     `ZeroMassBalState!` / `ZeroMassBalFlux!` / `ZeroFluxDiags!` (EDTypesMod),
#   * `nesterov_index` / `init_fire_weather!` (SFNesterovMod).
#
# Upstream-Fortran quirks preserved (see inline comments):
#   * The Julia `calculate_SP_properties` RETURNS (leaf_c, dbh, cohort_n, c_area)
#     rather than using out-args; the Fortran's hardcoded htop=0.5, tlai=0.2,
#     tsai=0.1, canopy_layer=1 SP-init values are kept.
#   * `cohort_n = initd*area`, then in nocomp multiplied by `sum(use_this_pft)`
#     so that nocomp runs (one PFT per patch) are comparable to competition runs.
#   * The bareground patch is created only in nocomp + fixed-biogeog mode, given
#     non-trivial bareground area.
#   * `area_error` (skipped tiny patches) is folded back into the surviving patch
#     areas via `area = area * AREA/(AREA - area_error)`; then a precision-level
#     mismatch with `AREA` is renormalized, a larger one is a fatal error.
#   * The SP-mode crown-area correction after patch resizing scales cohort `n`.
#   * Plant-hydraulics rhizosphere init (`updateSizeDepRhizHydProps`) is gated
#     behind `hlm_use_planthydro` (off by default) — stubbed as a TODO.
#   * The fixed-biogeog area_pft mapping reads `EDPftvarcon_inst[].hlm_pft_map`
#     and `bc_in.pft_areafrac(_lu)`; the NaN-guard makes everything bareground.
#   * `use_this_pft = itrue` for ALL pfts unless fixed-biogeog filters by area.

# ============================================================================
# init_site_vars!
# ============================================================================

"""
    init_site_vars!(site_in::ed_site_type, bc_in, bc_out)

Allocate the site-level arrays from the boundary conditions and attach the
fire-weather object.  Mirrors the Fortran `init_site_vars`.  Array dimensions
follow the Fortran (term/fmort/imort size-x-pft, damage classes gated on
`hlm_use_tree_damage`, soil layers from `bc_in.nlevsoil`, seed dispersal, SP
targets, flux diagnostics, area_pft / land-use / age arrays).
"""
function init_site_vars!(site_in::ed_site_type, bc_in, bc_out=nothing)
    npft = numpft[]
    nsc  = nlevsclass[]
    nmt  = n_term_mort_types

    site_in.term_nindivs_canopy = Array{Float64,3}(undef, nmt, nsc, npft)
    site_in.term_nindivs_ustory = Array{Float64,3}(undef, nmt, nsc, npft)
    site_in.demotion_rate       = Vector{Float64}(undef, nsc)
    site_in.promotion_rate      = Vector{Float64}(undef, nsc)
    site_in.imort_rate          = Matrix{Float64}(undef, nsc, npft)
    site_in.fmort_rate_canopy   = Matrix{Float64}(undef, nsc, npft)
    site_in.fmort_rate_ustory   = Matrix{Float64}(undef, nsc, npft)
    site_in.fmort_rate_cambial  = Matrix{Float64}(undef, nsc, npft)
    site_in.fmort_rate_crown    = Matrix{Float64}(undef, nsc, npft)
    site_in.growthflux_fusion   = Matrix{Float64}(undef, nsc, npft)
    site_in.mass_balance        = [site_massbal_type() for _ in 1:num_elements[]]
    site_in.iflux_balance       = [site_ifluxbal_type() for _ in 1:num_elements[]]

    if hlm_use_tree_damage[] == itrue
        nld = nlevdamage[]
        site_in.term_nindivs_canopy_damage = Array{Float64,3}(undef, nld, nsc, npft)
        site_in.term_nindivs_ustory_damage = Array{Float64,3}(undef, nld, nsc, npft)
        site_in.imort_rate_damage           = Array{Float64,3}(undef, nld, nsc, npft)
        site_in.imort_cflux_damage          = Matrix{Float64}(undef, nld, nsc)
        site_in.term_cflux_canopy_damage    = Matrix{Float64}(undef, nld, nsc)
        site_in.term_cflux_ustory_damage    = Matrix{Float64}(undef, nld, nsc)
        site_in.fmort_rate_canopy_damage    = Array{Float64,3}(undef, nld, nsc, npft)
        site_in.fmort_rate_ustory_damage    = Array{Float64,3}(undef, nld, nsc, npft)
        site_in.fmort_cflux_canopy_damage   = Matrix{Float64}(undef, nld, nsc)
        site_in.fmort_cflux_ustory_damage   = Matrix{Float64}(undef, nld, nsc)
    else
        site_in.term_nindivs_canopy_damage = Array{Float64,3}(undef, 1, 1, 1)
        site_in.term_nindivs_ustory_damage = Array{Float64,3}(undef, 1, 1, 1)
        site_in.imort_rate_damage           = Array{Float64,3}(undef, 1, 1, 1)
        site_in.imort_cflux_damage          = Matrix{Float64}(undef, 1, 1)
        site_in.term_cflux_canopy_damage    = Matrix{Float64}(undef, 1, 1)
        site_in.term_cflux_ustory_damage    = Matrix{Float64}(undef, 1, 1)
        site_in.fmort_rate_canopy_damage    = Array{Float64,3}(undef, 1, 1, 1)
        site_in.fmort_rate_ustory_damage    = Array{Float64,3}(undef, 1, 1, 1)
        site_in.fmort_cflux_canopy_damage   = Matrix{Float64}(undef, 1, 1)
        site_in.fmort_cflux_ustory_damage   = Matrix{Float64}(undef, 1, 1)
    end

    site_in.term_carbonflux_canopy = Matrix{Float64}(undef, nmt, npft)
    site_in.term_carbonflux_ustory = Matrix{Float64}(undef, nmt, npft)
    site_in.imort_carbonflux        = Vector{Float64}(undef, npft)
    site_in.fmort_carbonflux_canopy = Vector{Float64}(undef, npft)
    site_in.fmort_carbonflux_ustory = Vector{Float64}(undef, npft)

    site_in.term_abg_flux  = Matrix{Float64}(undef, nsc, npft)
    site_in.imort_abg_flux = Matrix{Float64}(undef, nsc, npft)
    site_in.fmort_abg_flux = Matrix{Float64}(undef, nsc, npft)

    site_in.nlevsoil    = bc_in.nlevsoil
    site_in.rootfrac_scr = Vector{Float64}(undef, site_in.nlevsoil)
    site_in.zi_soil      = Vector{Float64}(undef, site_in.nlevsoil + 1)  # 0:nlevsoil
    site_in.dz_soil      = Vector{Float64}(undef, site_in.nlevsoil)
    site_in.z_soil       = Vector{Float64}(undef, site_in.nlevsoil)

    site_in.area_PFT             = Matrix{Float64}(undef, npft, n_landuse_cats)
    site_in.landuse_vector_gt_min = Vector{Bool}(undef, n_landuse_cats)

    site_in.use_this_pft = Vector{Int}(undef, npft)
    site_in.area_by_age  = Vector{Float64}(undef, nlevage[])

    # For CNP dynamics, track the mean l2fr of recruits for different pfts/canopy
    # positions.
    site_in.rec_l2fr = Matrix{Float64}(undef, npft, nclmax)

    # SP mode
    site_in.sp_tlai = Vector{Float64}(undef, npft)
    site_in.sp_tsai = Vector{Float64}(undef, npft)
    site_in.sp_htop = Vector{Float64}(undef, npft)

    # Site-level flux diagnostics
    site_in.flux_diags = site_fluxdiags_type()
    site_in.flux_diags.elem = [elem_diag_type() for _ in 1:num_elements[]]
    for el in 1:num_elements[]
        site_in.flux_diags.elem[el].surf_fine_litter_input = Vector{Float64}(undef, npft)
        site_in.flux_diags.elem[el].root_litter_input      = Vector{Float64}(undef, npft)
    end
    site_in.flux_diags.nh4_uptake_scpf = Vector{Float64}(undef, npft * nsc)
    site_in.flux_diags.no3_uptake_scpf = Vector{Float64}(undef, npft * nsc)
    site_in.flux_diags.sym_nfix_scpf   = Vector{Float64}(undef, npft * nsc)
    site_in.flux_diags.n_efflux_scpf   = Vector{Float64}(undef, npft * nsc)
    site_in.flux_diags.p_uptake_scpf   = Vector{Float64}(undef, npft * nsc)
    site_in.flux_diags.p_efflux_scpf   = Vector{Float64}(undef, npft * nsc)

    # Initialize the static soil arrays from the boundary (initial) condition.
    copyto!(site_in.zi_soil, bc_in.zi_sisl)
    copyto!(site_in.dz_soil, bc_in.dz_sisl)
    copyto!(site_in.z_soil,  bc_in.z_sisl)

    # Seed dispersal
    site_in.seed_in  = Vector{Float64}(undef, npft)
    site_in.seed_out = Vector{Float64}(undef, npft)

    # Fire weather object (Nesterov index is the one ported option).
    site_in.fireWeather = nesterov_index()
    init_fire_weather!(site_in.fireWeather)

    return nothing
end

# ============================================================================
# zero_site!
# ============================================================================

"""
    zero_site!(site_in::ed_site_type)

Zero / unset every site-level accumulator and phenology / fire / disturbance /
mass-balance / flux-diagnostic state.  Mirrors the Fortran `zero_site`.
`fates_unset_int`/`NaN` are used exactly where the Fortran uses `fates_unset_int`
/ `nan`.
"""
function zero_site!(site_in::ed_site_type)
    site_in.oldest_patch   = nothing
    site_in.youngest_patch = nothing

    # PHENOLOGY
    site_in.cstatus       = fates_unset_int
    fill!(site_in.dstatus, fates_unset_int)
    site_in.grow_deg_days = NaN
    site_in.snow_depth    = NaN
    site_in.nchilldays    = fates_unset_int
    site_in.ncolddays     = fates_unset_int
    site_in.cleafondate   = fates_unset_int
    site_in.cleafoffdate  = fates_unset_int
    fill!(site_in.dleafondate,  fates_unset_int)
    fill!(site_in.dleafoffdate, fates_unset_int)
    site_in.cndaysleafon  = fates_unset_int
    site_in.cndaysleafoff = fates_unset_int
    fill!(site_in.dndaysleafon,  fates_unset_int)
    fill!(site_in.dndaysleafoff, fates_unset_int)
    fill!(site_in.elong_factor, NaN)

    fill!(site_in.liqvol_memory, NaN)
    fill!(site_in.smp_memory, NaN)
    fill!(site_in.vegtemp_memory, NaN)

    site_in.phen_model_date = fates_unset_int

    # Disturbance rates tracking
    site_in.primary_land_patchfusion_error = 0.0
    fill!(site_in.disturbance_rates, 0.0)
    fill!(site_in.landuse_transition_matrix, 0.0)

    # FIRE
    site_in.fdi           = 0.0   # daily fire danger index (0-1)
    site_in.NF            = 0.0   # daily lightning strikes per km2
    site_in.NF_successful = 0.0   # daily successful ignitions per km2

    for el in 1:num_elements[]
        ZeroMassBalState!(site_in.mass_balance[el])
        ZeroMassBalFlux!(site_in.mass_balance[el])
    end

    ZeroFluxDiags!(site_in.flux_diags)

    # This will be initialized in FatesSoilBGCFluxMod:PrepCH4BCs(); below -9000
    # triggers the smoother's first-value path.
    site_in.ema_npp = -9999.9

    # termination and recruitment info
    fill!(site_in.term_nindivs_canopy, 0.0)
    fill!(site_in.term_nindivs_ustory, 0.0)
    site_in.term_crownarea_canopy = 0.0
    site_in.term_crownarea_ustory = 0.0
    site_in.imort_crownarea       = 0.0
    site_in.fmort_crownarea_canopy = 0.0
    site_in.fmort_crownarea_ustory = 0.0
    fill!(site_in.term_carbonflux_canopy, 0.0)
    fill!(site_in.term_carbonflux_ustory, 0.0)
    fill!(site_in.recruitment_rate, 0.0)
    fill!(site_in.imort_rate, 0.0)
    fill!(site_in.imort_carbonflux, 0.0)
    fill!(site_in.fmort_rate_canopy, 0.0)
    fill!(site_in.fmort_rate_ustory, 0.0)
    fill!(site_in.fmort_carbonflux_canopy, 0.0)
    fill!(site_in.fmort_carbonflux_ustory, 0.0)
    fill!(site_in.fmort_rate_cambial, 0.0)
    fill!(site_in.fmort_rate_crown, 0.0)
    fill!(site_in.term_abg_flux, 0.0)
    fill!(site_in.imort_abg_flux, 0.0)
    fill!(site_in.fmort_abg_flux, 0.0)

    # fusion-induced growth flux of individuals
    fill!(site_in.growthflux_fusion, 0.0)

    # demotion/promotion info
    fill!(site_in.demotion_rate, 0.0)
    site_in.demotion_carbonflux = 0.0
    fill!(site_in.promotion_rate, 0.0)
    site_in.promotion_carbonflux = 0.0

    # damage transition info
    fill!(site_in.imort_rate_damage, 0.0)
    fill!(site_in.term_nindivs_canopy_damage, 0.0)
    fill!(site_in.term_nindivs_ustory_damage, 0.0)
    fill!(site_in.imort_cflux_damage, 0.0)
    fill!(site_in.term_cflux_canopy_damage, 0.0)
    fill!(site_in.term_cflux_ustory_damage, 0.0)
    site_in.crownarea_canopy_damage = 0.0
    site_in.crownarea_ustory_damage = 0.0
    fill!(site_in.fmort_rate_canopy_damage, 0.0)
    fill!(site_in.fmort_rate_ustory_damage, 0.0)
    fill!(site_in.fmort_cflux_canopy_damage, 0.0)
    fill!(site_in.fmort_cflux_ustory_damage, 0.0)

    # Resources management (logging/harvesting, etc)
    site_in.resources_management.harvest_debt     = 0.0
    site_in.resources_management.harvest_debt_sec = 0.0
    site_in.resources_management.trunk_product_site = 0.0

    # canopy spread
    site_in.spread = 0.0

    fill!(site_in.area_PFT, 0.0)
    site_in.area_bareground = 0.0

    # Seed dispersal
    fill!(site_in.seed_in, 0.0)
    fill!(site_in.seed_out, 0.0)

    fill!(site_in.use_this_pft, fates_unset_int)
    fill!(site_in.area_by_age, 0.0)

    site_in.transition_landuse_from_off_to_on = false

    return nothing
end

# ============================================================================
# set_site_properties!
# ============================================================================

"""
    set_site_properties!(nsites, sites, bc_in)

Set the initial site properties on a cold start (`hlm_is_restart == ifalse`):
phenology dates/status/GDD, water + veg-temperature memory, fire ignition rates,
the recruit l2fr, and — in fixed-biogeography mode — the per-PFT-per-land-use
area mapping (with the nocomp PFT-count trimming + per-land-use renormalization)
and `min_allowed_landuse_fraction`.  On a restart this routine is a no-op (the
values are read in afterward).  Mirrors the Fortran `set_site_properties`.
"""
function set_site_properties!(nsites::Integer, sites::AbstractVector, bc_in::AbstractVector)
    # If this is a restart, leave everything unset and let the restart read fill
    # the values after this routine.
    if hlm_is_restart[] != ifalse
        return nothing
    end

    GDD        = 30.0
    cleafon    = 100
    cleafoff   = 300
    cndleafon  = 0
    cndleafoff = 0
    cstat      = phen_cstat_notcold      # Leaves are on
    dstat      = phen_dstat_moiston      # Leaves are on
    dleafoff   = 300
    dleafon    = 100
    dndleafon  = 0
    dndleafoff = 0
    liqvolmem  = 0.5
    smpmem     = 0.0
    elong_fac  = 1.0

    npft = numpft[]

    for s in 1:nsites
        site = sites[s]
        site.nchilldays = 0
        site.ncolddays  = 0   # recalculated in phenology immediately, so memory-less,
                              # but needed for first history-file value
        site.phen_model_date = 0
        site.cleafondate     = cleafon  - hlm_day_of_year[]
        site.cleafoffdate    = cleafoff - hlm_day_of_year[]
        site.cndaysleafon    = cndleafon
        site.cndaysleafoff   = cndleafoff
        for ft in 1:npft
            site.dleafoffdate[ft]  = dleafoff - hlm_day_of_year[]
            site.dleafondate[ft]   = dleafon  - hlm_day_of_year[]
            site.dndaysleafon[ft]  = dndleafon
            site.dndaysleafoff[ft] = dndleafoff
        end
        site.grow_deg_days = GDD

        for ft in 1:npft
            for w in 1:numWaterMem
                site.liqvol_memory[w, ft] = liqvolmem
                site.smp_memory[w, ft]    = smpmem
            end
        end
        for t in 1:num_vegtemp_mem
            site.vegtemp_memory[t] = 0.0
        end

        site.cstatus = cstat
        for ft in 1:npft
            site.dstatus[ft]      = dstat
            site.elong_factor[ft] = elong_fac
        end

        site.NF            = 0.0
        site.NF_successful = 0.0
        fill!(site.area_PFT, 0.0)

        for ft in 1:npft
            site.rec_l2fr[ft, :] .= prt_params.allom_l2fr[ft]
        end

        # Difficult to choose a reasonable starting smoothing value; cold-start to -1.
        site.ema_npp = -9999.0

        if hlm_use_fixed_biogeog[] == itrue
            _set_site_area_pft!(site, bc_in[s])
        end

        for ft in 1:npft
            # Setting to true ensures all pfts are used for nocomp w/ no biogeog.
            site.use_this_pft[ft] = itrue
            if hlm_use_fixed_biogeog[] == itrue
                if any(site.area_PFT[ft, :] .> 0.0)
                    site.use_this_pft[ft] = itrue
                else
                    site.use_this_pft[ft] = ifalse
                end
            end
        end

        # Minimum allowable land-use fraction: a function of the minimum allowable
        # patch size, and (in nocomp) the bare-ground + min pft fractions.
        if hlm_use_nocomp[] == itrue
            if (1.0 - site.area_bareground) > nearzero
                site.min_allowed_landuse_fraction =
                    min_patch_area_forced /
                    (area * min_nocomp_pftfrac_perlanduse * (1.0 - site.area_bareground))
            else
                # All bare ground; doesn't matter, set to one to ignore land use.
                site.min_allowed_landuse_fraction = 1.0
            end
        else
            site.min_allowed_landuse_fraction = min_patch_area_forced / area
        end
    end

    return nothing
end

# Helper: the fixed-biogeog area_pft mapping block of set_site_properties.
# Maps HLM PFTs onto FATES PFTs (per land-use type, with/without LUH2), guards
# NaNs to bareground, removes tiny patches, trims to max nocomp pfts, and
# renormalizes per land-use type.
function _set_site_area_pft!(site::ed_site_type, bc_in)
    npft  = numpft[]
    edp   = ed_params()
    evcon = EDPftvarcon_inst[]

    if hlm_use_luh[] == itrue
        # MAPPING of FATES PFTs onto HLM PFTs WITH land use.  pft_areafrac_lu is the
        # area in each HLM PFT and land-use type (surface dataset); hlm_pft_map maps
        # that into each FATES PFT (param file).
        if !(any(isnan, bc_in.pft_areafrac_lu) || isnan(bc_in.baregroundfrac))
            for i_lu in 1:n_landuse_cats
                if !is_crop[i_lu]
                    for hlm_pft in 1:size(evcon.hlm_pft_map, 2)
                        for fates_pft in 1:npft
                            site.area_PFT[fates_pft, i_lu] +=
                                evcon.hlm_pft_map[fates_pft, hlm_pft] *
                                bc_in.pft_areafrac_lu[hlm_pft, i_lu]
                        end
                    end
                else
                    # Crops: pft_areafrac_lu only exists for natural PFTs, so use the
                    # crop_lu_pft_vector mapping and set that PFT's area to 1.
                    site.area_PFT[edp.crop_lu_pft_vector[i_lu], i_lu] = 1.0
                end
            end
            site.area_bareground = bc_in.baregroundfrac
        else
            # NaN guard: make everything bare ground.
            site.area_bareground = 1.0
            fill!(site.area_PFT, 0.0)
            fates_log_println("Nan values for pftareafrac. dumping site info.")
            dump_site(site)
        end
    else
        # MAPPING of FATES PFTs onto HLM PFTs (no land use).  area is mapped into the
        # primaryland land-use index.
        for hlm_pft in 1:size(evcon.hlm_pft_map, 2)
            for fates_pft in 1:npft
                site.area_PFT[fates_pft, primaryland] +=
                    evcon.hlm_pft_map[fates_pft, hlm_pft] *
                    bc_in.pft_areafrac[hlm_pft]
            end
            # NOTE (Fortran quirk): area_bareground is assigned inside the hlm_pft loop
            # from pft_areafrac(0), so it takes the loop's last value.  bc_in.pft_areafrac
            # is 1-based here; the FATES "PFT 0" bareground fraction is carried in
            # bc_in.baregroundfrac, which we read instead.
            site.area_bareground = bc_in.baregroundfrac
        end
    end

    # Handle edge cases.
    for i_lu in 1:n_landuse_cats
        for ft in 1:npft
            # Remove tiny patches to prevent numerical errors in terminate patches.
            if site.area_PFT[ft, i_lu] < min_nocomp_pftfrac_perlanduse &&
               site.area_PFT[ft, i_lu] > nearzero
                site.area_PFT[ft, i_lu] = 0.0
            end
            # Negative area -> end run.
            if site.area_PFT[ft, i_lu] < 0.0
                fates_endrun("negative area ft=$(ft) i_lu=$(i_lu): $(site.area_PFT[ft, i_lu])")
            end
        end
    end

    # In nocomp mode, if the number of nocomp PFTs of a land-use type exceeds the
    # maximum number of patches allowed, keep only that many PFTs, starting with the
    # largest-area PFTs and working down.
    if hlm_use_nocomp[] == itrue
        for i_lu in 1:n_landuse_cats
            maxnp = edp.max_nocomp_pfts_by_landuse[i_lu]
            if count(>(0.0), @view site.area_PFT[:, i_lu]) > maxnp
                temp_vec = zeros(npft)
                for _ in 1:maxnp
                    # MAXLOC: index of the current largest-area PFT.
                    imax = argmax(@view site.area_PFT[:, i_lu])
                    temp_vec[imax] = site.area_PFT[imax, i_lu]
                    site.area_PFT[imax, i_lu] = 0.0
                end
                site.area_PFT[:, i_lu] .= temp_vec
            end
        end
    end

    # Re-normalize PFT area so it sums to one for each (active) land-use type; track
    # bare ground separately in nocomp.
    for i_lu in 1:n_landuse_cats
        sumarea = sum(@view site.area_PFT[:, i_lu])
        if sumarea > nearzero
            site.area_PFT[:, i_lu] ./= sumarea
        else
            # No PFT area in primary lands -> set bare-ground fraction to one.
            if i_lu == primaryland
                site.area_bareground = 1.0
                site.area_PFT[:, i_lu] .= 0.0
            end
        end
    end

    return nothing
end

# ============================================================================
# init_patches!
# ============================================================================

"""
    init_patches!(nsites, sites, bc_in)

Create the starting patch(es) for each site (the near-bare-ground cold-start),
seed each vegetated patch with cohorts, area-account the patches to the site
`area`, initialize the mass-balance stock, and renumber.  This is the default
(non-inventory) path; the inventory path delegates to
[`initialize_sites_by_inventory!`](@ref) (Batch 16).  Mirrors the Fortran
`init_patches`.

Inventory init is signalled by `hlm_use_inventory_init == itrue`, in which case
the caller must supply the inventory `sitelist_lines` keyword (+ pss/css
providers); the default cold-start needs none.
"""
function init_patches!(nsites::Integer, sites::AbstractVector, bc_in::AbstractVector;
                       sitelist_lines=nothing,
                       pss_provider=nothing, css_provider=nothing)

    age = 0.0   # notional age of a cold-start patch

    if hlm_use_inventory_init[] == itrue
        # ----- Inventory based cold-start -----
        for s in 1:nsites
            # Closed-canopy forest inventories likely have smaller spread factors
            # than bare ground.
            sites[s].spread = init_spread_inventory
        end

        initialize_sites_by_inventory!(nsites, sites, bc_in;
                                       sitelist_lines = sitelist_lines,
                                       pss_provider = pss_provider,
                                       css_provider = css_provider)

        # Initialize the total carbon stock for mass-balance checks.
        for s in 1:nsites
            for el in 1:num_elements[]
                _, biomass_stock, litter_stock, seed_stock = SiteMassStock(sites[s], el)
                sites[s].mass_balance[el].old_stock =
                    biomass_stock + litter_stock + seed_stock
                sites[s].iflux_balance[el].iflux_liveveg = (biomass_stock + seed_stock) * area_inv
                sites[s].iflux_balance[el].iflux_litter  = litter_stock * area_inv
            end
            set_patchno(sites[s])
        end
    else
        # ----- Near-Bare-Ground (NBG) cold-start -----
        num_nocomp_pfts = (hlm_use_nocomp[] == itrue) ? numpft[] : 1

        for s in 1:nsites
            site = sites[s]
            fill!(site.sp_tlai, 0.0)
            fill!(site.sp_tsai, 0.0)
            fill!(site.sp_htop, 0.0)

            site.spread = init_spread_near_bare_ground

            # Read LUH state data to determine the initial land-use types.
            state_vector = zeros(n_landuse_cats)
            if hlm_use_luh[] == itrue
                # TODO Batch NN: GetLUHStatedata (FatesLandUseChangeMod) is not yet
                # ported; LUH2 init is off by default (hlm_use_luh defaults to unset/
                # ifalse).  When ported, fill state_vector + landuse_vector_gt_min here.
                state_vector .= GetLUHStatedata(bc_in[s])
                for i_lu in 1:n_landuse_cats
                    site.landuse_vector_gt_min[i_lu] =
                        state_vector[i_lu] > site.min_allowed_landuse_fraction
                end
            else
                # No LUH2: initialize with primarylands (index 1) only.
                state_vector[primaryland] = 1.0
            end

            # Confirm state vector sums to 1.
            if abs(sum(state_vector) - 1.0) > rsnbl_math_prec
                fates_endrun("state vector must sum to 1: $(sum(state_vector))")
            end

            is_first_patch = true
            area_error = 0.0

            # First, a bare-ground patch (nocomp + fixed-biogeog only).
            if hlm_use_nocomp[] == itrue && hlm_use_fixed_biogeog[] == itrue
                newparea = area * site.area_bareground
                if newparea > min_patch_area_forced
                    newp = fates_patch_type()
                    Create!(newp, age, newparea, nocomp_bareground_land, nocomp_bareground,
                            num_swb, numpft[], site.nlevsoil, hlm_current_tod[],
                            ed_params().regeneration_model)
                    newp.patchno = 1
                    newp.younger = nothing
                    newp.older   = nothing
                    site.youngest_patch = newp
                    site.oldest_patch   = newp
                    is_first_patch = false

                    litt_init = (hlm_use_sp[] == itrue) ? fates_unset_r8 : 0.0
                    for el in 1:num_elements[]
                        InitConditions!(newp.litter[el], litt_init, litt_init,
                                        litt_init, litt_init, litt_init, litt_init)
                    end
                else
                    area_error += newparea
                end
            end

            end_landuse_idx = (hlm_use_luh[] == itrue) ? n_landuse_cats : 1

            # Create the non-bareground (vegetated) patches when either: we are not
            # doing both nocomp & fixed-biogeog, OR there is non-zero bareground area.
            if (1.0 - site.area_bareground) > nearzero ||
               !(hlm_use_nocomp[] == itrue && hlm_use_fixed_biogeog[] == itrue)

                for i_lu_state in 1:end_landuse_idx
                    if state_vector[i_lu_state] > nearzero
                        for n in 1:num_nocomp_pfts
                            nocomp_pft = (hlm_use_nocomp[] == itrue) ? n : fates_unset_int

                            if hlm_use_nocomp[] == itrue
                                if hlm_use_fixed_biogeog[] == itrue
                                    newparea = site.area_PFT[nocomp_pft, i_lu_state] *
                                               area * state_vector[i_lu_state] *
                                               (1.0 - site.area_bareground)
                                else
                                    newparea = area * state_vector[i_lu_state] / numpft[]
                                end
                            else
                                # Default: one patch with the area of the whole site.
                                newparea = area * state_vector[i_lu_state]
                            end

                            # Stop tiny patches being initialized (PFT not present in nocomp).
                            if newparea > min_patch_area_forced
                                newp = fates_patch_type()
                                Create!(newp, age, newparea, i_lu_state, nocomp_pft,
                                        num_swb, numpft[], site.nlevsoil, hlm_current_tod[],
                                        ed_params().regeneration_model)

                                if is_first_patch
                                    newp.patchno = 1
                                    newp.younger = nothing
                                    newp.older   = nothing
                                    site.youngest_patch = newp
                                    site.oldest_patch   = newp
                                    is_first_patch = false
                                else
                                    # N>1 patches (only when nocomp or land use is on).
                                    # The new patch is the 'youngest', arbitrarily.
                                    newp.patchno = nocomp_pft + (i_lu_state - 1) * numpft[]
                                    newp.older   = site.youngest_patch
                                    newp.younger = nothing
                                    site.youngest_patch.younger = newp
                                    site.youngest_patch = newp
                                end

                                litt_init = (hlm_use_sp[] == itrue) ? fates_unset_r8 : 0.0
                                for el in 1:num_elements[]
                                    InitConditions!(newp.litter[el], litt_init, litt_init,
                                                    litt_init, litt_init, litt_init, litt_init)
                                end

                                init_cohorts!(site, newp, bc_in[s])
                            else
                                area_error += newparea
                            end
                        end
                    end
                end
            end

            # If we skipped small patches above, resize surviving patches.
            if area_error > nearzero
                newp = site.oldest_patch
                while newp !== nothing
                    newp.area = newp.area * area / (area - area_error)
                    newp = newp.younger
                end
            end

            # Check the total area adds to the site area.
            total = 0.0
            newp = site.oldest_patch
            while newp !== nothing
                total += newp.area
                newp = newp.younger
            end

            area_diff = total - area
            if abs(area_diff) > nearzero
                if abs(area_diff) < area_error_4
                    # Precision error: renormalize all patch areas to sum to AREA.
                    newp = site.oldest_patch
                    while newp !== nothing
                        newp.area = newp.area * (area / total)
                        newp = newp.younger
                    end
                else
                    # A big error, not just a precision error.
                    fates_endrun(string("issue with patch area in EDinit ",
                        area_diff, " ", total, " lat=", site.lat, " lon=", site.lon))
                end
            end

            # We might have messed up crown areas now — correct if in SP mode.
            if hlm_use_sp[] == itrue
                newp = site.oldest_patch
                while newp !== nothing
                    cohort = newp.tallest
                    while cohort !== nothing
                        if abs(cohort.c_area - newp.area) < area_error_3
                            old_carea = cohort.c_area
                            cohort.c_area = cohort.c_area - (cohort.c_area - newp.area)
                            cohort.n = cohort.n * (cohort.c_area / old_carea)
                        end
                        cohort = cohort.shorter
                    end
                    newp = newp.younger
                end
            end

            # Initialize the total carbon stock for mass-balance checks.
            for el in 1:num_elements[]
                _, biomass_stock, litter_stock, seed_stock = SiteMassStock(site, el)
                site.mass_balance[el].old_stock = biomass_stock + litter_stock + seed_stock
                site.iflux_balance[el].iflux_liveveg = (biomass_stock + seed_stock) * area_inv
                site.iflux_balance[el].iflux_litter  = litter_stock * area_inv
            end

            set_patchno(site)
        end
    end

    # Zero all the patch fire variables for the first timestep.
    for s in 1:nsites
        currentPatch = sites[s].youngest_patch
        while currentPatch !== nothing
            currentPatch.livegrass  = 0.0
            currentPatch.ros_front  = 0.0
            currentPatch.tau_l      = 0.0
            currentPatch.tfc_ros    = 0.0
            currentPatch.fi         = 0.0
            currentPatch.fire       = 0
            currentPatch.fd         = 0.0
            currentPatch.ros_back   = 0.0
            fill!(currentPatch.scorch_ht, 0.0)
            currentPatch.frac_burnt = 0.0
            currentPatch = currentPatch.older
        end
    end

    # Set the rhizosphere shells based on the plant initialization.  The hydraulics
    # state was set inside init_cohorts!->create_cohort.
    if hlm_use_planthydro[] == itrue
        for s in 1:nsites
            _init_patches_rhiz_hydro!(sites[s], bc_in[s])
        end
    end

    # Make sure there are no very tiny patches.
    for s in 1:nsites
        currentPatch = sites[s].youngest_patch
        while currentPatch !== nothing
            if currentPatch.area < min_patch_area_forced
                fates_endrun(string("edinit somehow making tiny patches ",
                    currentPatch.land_use_label, " ", currentPatch.nocomp_pft_label,
                    " ", currentPatch.area))
            end
            currentPatch = currentPatch.older
        end
    end

    return nothing
end

# Plant-hydraulics rhizosphere init (only entered when hlm_use_planthydro==itrue).
function _init_patches_rhiz_hydro!(site::ed_site_type, bc_in)
    # Set the rhizosphere shell geometry/conductances from the plant init state.
    UpdateSizeDepRhizHydProps(site, bc_in)
    return nothing
end

# ============================================================================
# init_cohorts!
# ============================================================================

"""
    init_cohorts!(site_in::ed_site_type, patch_in::fates_patch_type, bc_in)

Seed a patch with the initial cohorts on bare ground: one small cohort of each
"used" PFT (determined by the fixed-biogeog site filter x the nocomp patch
filter), sized from `EDPftvarcon_inst.initd` (positive => initial density;
negative => initial dbh, nocomp only), with phenology / elongation factors set
from the PFT deciduous type + site cold/drought status, biomass from allometry,
PARTEH per-element state, and inserted via `create_cohort`.  Then (non-SP) the
cohorts are fused and sorted.  Mirrors the Fortran `init_cohorts`.
"""
function init_cohorts!(site_in::ed_site_type, patch_in::fates_patch_type, bc_in)
    npft  = numpft[]
    evcon = EDPftvarcon_inst[]

    patch_in.tallest  = nothing
    patch_in.shortest = nothing

    # If any pft starts at a large (dbh-specified) size, the whole site needs spread 0.
    for pft in 1:npft
        if evcon.initd[pft] < 0.0
            site_in.spread = init_spread_inventory
        end
    end

    # Manage interactions of fixed-biogeog (site filter) and nocomp (patch filter):
    #   1. biogeog=false, nocomp=false: all PFTs on (DEFAULT)
    #   2. biogeog=true,  nocomp=false: site-level filter
    #   3. biogeog=false, nocomp=true : patch-level filter
    #   4. biogeog=true,  nocomp=true : patch + site filter
    use_pft_local = fill(itrue, npft)
    for pft in 1:npft
        use_pft_local[pft] = itrue   # Case 1
        if hlm_use_fixed_biogeog[] == itrue
            use_pft_local[pft] = site_in.use_this_pft[pft]   # Case 2
            if hlm_use_nocomp[] == itrue && pft != patch_in.nocomp_pft_label
                use_pft_local[pft] = ifalse   # Case 3
            end
        else
            if hlm_use_nocomp[] == itrue && pft != patch_in.nocomp_pft_label
                use_pft_local[pft] = ifalse   # Case 4
            end
        end
    end

    recruitstatus = 0    # newly created cohorts are initialized, not recruited
    zero_co_age   = 0.0  # the age of a newly recruited cohort is zero

    for pft in 1:npft
        use_pft_local[pft] == itrue || continue

        l2fr         = prt_params.allom_l2fr[pft]
        canopy_trim  = 1.0
        crown_damage = 1   # assume undamaged

        fnrt_drop_fraction = prt_params.phen_fnrt_drop_fraction[pft]
        stem_drop_fraction = prt_params.phen_stem_drop_fraction[pft]

        # Phenology variables.
        local efleaf_coh::Float64, effnrt_coh::Float64, efstem_coh::Float64, leaf_status::Int
        if hlm_use_sp[] == itrue
            # Satellite phenology: do not override SP values with built-in phenology.
            efleaf_coh = 1.0
            effnrt_coh = 1.0
            efstem_coh = 1.0
            leaf_status = leaves_on
        else
            if prt_params.season_decid[pft] == itrue &&
               (site_in.cstatus == phen_cstat_nevercold || site_in.cstatus == phen_cstat_iscold)
                # Cold deciduous, off season: complete abscission.
                efleaf_coh = 0.0
                effnrt_coh = 1.0 - fnrt_drop_fraction
                efstem_coh = 1.0 - stem_drop_fraction
                leaf_status = leaves_off
            elseif prt_params.stress_decid[pft] == ihard_stress_decid ||
                   prt_params.stress_decid[pft] == isemi_stress_decid
                # Drought deciduous: tissues other than leaves use a combination of
                # the elongation factor (e) and the drop fraction (x) so the remaining
                # biomass = e when x=1, and original when x=0.
                efleaf_coh = site_in.elong_factor[pft]
                effnrt_coh = 1.0 - (1.0 - efleaf_coh) * fnrt_drop_fraction
                efstem_coh = 1.0 - (1.0 - efleaf_coh) * stem_drop_fraction
                if efleaf_coh > 0.0
                    leaf_status = leaves_on
                else
                    leaf_status = leaves_off
                end
            else
                # Evergreens, or deciduous during growing season: fully flushed.
                efleaf_coh = 1.0
                effnrt_coh = 1.0
                efstem_coh = 1.0
                leaf_status = leaves_on
            end
        end

        local dbh::Float64, height::Float64, cohort_n::Float64, c_area::Float64, c_leaf::Float64

        # Positive initd => initial density (calc diameter). Negative => initial dbh
        # (calc density, nocomp only).
        if evcon.initd[pft] > nearzero
            cohort_n = evcon.initd[pft] * patch_in.area
            if hlm_use_nocomp[] == itrue
                # nocomp has only one PFT per patch (vs numpft); scale density up by
                # the number of used PFTs so runs are comparable to competition mode.
                cohort_n = cohort_n * sum(site_in.use_this_pft)
            end
            height = evcon.hgt_min[pft]

            if hlm_use_sp[] == itrue
                # At this point we don't know the bc_in tlai/tsai/htop, so this is an
                # arbitrary first-timestep value (Fortran htop=0.5, tlai=0.2, tsai=0.1).
                height = 0.5
                c_leaf, dbh, cohort_n, c_area =
                    calculate_SP_properties(height, 0.2, 0.1, patch_in.area, pft,
                                            crown_damage, 1, evcon.vcmax25top[pft, 1])
            else
                dbh, _ = h2d_allom(height, pft)
                c_leaf, _ = bleaf(dbh, pft, crown_damage, canopy_trim, efleaf_coh)
                c_area = 0.0   # (computed inside create_cohort/Create for non-SP)
            end
        else
            # Interpret as initial diameter and calculate density (nocomp only).
            if hlm_use_nocomp[] == itrue
                dbh = abs(evcon.initd[pft])
                # Crown area of a single plant.
                _, c_area = carea_allom(dbh, 1.0, init_spread_inventory, pft, crown_damage)
                # Initial density required to close the canopy.
                cohort_n = patch_in.area / c_area
                c_leaf, _ = bleaf(dbh, pft, crown_damage, canopy_trim, efleaf_coh)
                # Crown area of the cohort.
                _, c_area = carea_allom(dbh, cohort_n, init_spread_inventory, pft, crown_damage)
                height, _ = h_allom(dbh, pft)
            else
                fates_endrun("Negative fates_recruit_init_density can only be used in no comp mode")
                dbh = 0.0; height = 0.0; cohort_n = 0.0; c_area = 0.0; c_leaf = 0.0
            end
        end

        # Allometric biomass pools (carbon).
        c_agw, _  = bagw_allom(dbh, pft, crown_damage, efstem_coh)
        c_bgw, _  = bbgw_allom(dbh, pft, efstem_coh)
        c_fnrt, _ = bfineroot(dbh, pft, canopy_trim, l2fr, effnrt_coh)
        _, c_sapw, _ = bsap_allom(dbh, pft, crown_damage, canopy_trim, efstem_coh)
        c_struct, _  = bdead_allom(c_agw, c_bgw, c_sapw, pft)
        c_store, _   = bstore_allom(dbh, pft, crown_damage, canopy_trim)

        # Initialize the mass of every element in every organ.
        prt = InitPRTObject!()

        for el in 1:num_elements[]
            element_id = element_list[el]
            if element_id == carbon12_element
                m_struct = c_struct
                m_leaf   = c_leaf
                m_fnrt   = c_fnrt
                m_sapw   = c_sapw
                m_store  = c_store
                m_repro  = 0.0
            elseif element_id == nitrogen_element
                m_struct = c_struct * prt_params.nitr_stoich_p1[pft, prt_params.organ_param_id[struct_organ]]
                m_leaf   = c_leaf   * prt_params.nitr_stoich_p1[pft, prt_params.organ_param_id[leaf_organ]]
                m_fnrt   = c_fnrt   * prt_params.nitr_stoich_p1[pft, prt_params.organ_param_id[fnrt_organ]]
                m_sapw   = c_sapw   * prt_params.nitr_stoich_p1[pft, prt_params.organ_param_id[sapw_organ]]
                m_repro  = 0.0
                m_store  = StorageNutrientTarget(pft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
            elseif element_id == phosphorus_element
                m_struct = c_struct * prt_params.phos_stoich_p1[pft, prt_params.organ_param_id[struct_organ]]
                m_leaf   = c_leaf   * prt_params.phos_stoich_p1[pft, prt_params.organ_param_id[leaf_organ]]
                m_fnrt   = c_fnrt   * prt_params.phos_stoich_p1[pft, prt_params.organ_param_id[fnrt_organ]]
                m_sapw   = c_sapw   * prt_params.phos_stoich_p1[pft, prt_params.organ_param_id[sapw_organ]]
                m_repro  = 0.0
                m_store  = StorageNutrientTarget(pft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
            else
                m_struct = 0.0; m_leaf = 0.0; m_fnrt = 0.0
                m_sapw = 0.0; m_store = 0.0; m_repro = 0.0
            end

            if hlm_parteh_mode[] == prt_carbon_allom_hyp ||
               hlm_parteh_mode[] == prt_cnp_flex_allom_hyp
                # All leaf mass into the first bin; the rest of the bins zeroed.
                SetState!(prt, leaf_organ, element_id, m_leaf, 1)
                for iage in 2:nleafage[]
                    SetState!(prt, leaf_organ, element_id, 0.0, iage)
                end
                SetState!(prt, fnrt_organ,   element_id, m_fnrt)
                SetState!(prt, sapw_organ,   element_id, m_sapw)
                SetState!(prt, store_organ,  element_id, m_store)
                SetState!(prt, struct_organ, element_id, m_struct)
                SetState!(prt, repro_organ,  element_id, m_repro)
            else
                fates_endrun("Unspecified PARTEH module during create_cohort")
            end
        end

        CheckInitialConditions(prt)

        create_cohort(site_in, patch_in, pft, cohort_n, height, zero_co_age, dbh, prt,
                      efleaf_coh, effnrt_coh, efstem_coh, leaf_status, recruitstatus,
                      canopy_trim, c_area, 1, crown_damage, site_in.spread, bc_in)
    end

    if hlm_use_sp[] == ifalse
        fuse_cohorts(site_in, patch_in, bc_in)
        sort_cohorts(patch_in)
    end

    return nothing
end

# Small log helper (mirrors Fortran write(fates_log(),*) ...).
fates_log_println(msg) = println(fates_log(), msg)
