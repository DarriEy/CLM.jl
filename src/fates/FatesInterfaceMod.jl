# FatesInterfaceMod.jl
# Julia port of FATES src/fates/main/FatesInterfaceMod.F90  (Tier F, Batch 18).
#
# THE host<->FATES coupling module — the seam where FATES meets clm_driver / the
# host land model. It:
#   * sets the FATES control parameters from the HLM (set_fates_ctrlparms — reads
#     all the hlm_* flags carried in FatesInterfaceTypesMod),
#   * sizes the FATES global dimensions (SetFatesGlobalElements1/2),
#   * allocates / zeroes / fills the boundary-condition structures
#     (allocate_bcin/bcout/bcpconst, zero_bcs, set_bcs, set_bcpconst),
#   * sets the FATES clock (SetFatesTime),
#   * defines the running-mean time-averaging globals (InitTimeAveragingGlobals),
#   * does the per-timestep running-mean update (UpdateFatesRMeansTStep),
#   * builds the grid-cell-neighbor list for seed dispersal
#     (DetermineGridCellNeighbors).
#
# ------------------------------------------------------------------------------
# Translation notes / conventions (CLAUDE.md):
#   * fates_r8 -> Float64. Fortran subroutines -> bang functions. Fortran names
#     preserved for traceability.
#   * The Fortran module-global `save` state lives in FatesInterfaceTypesMod as
#     `const Ref{...}` (the hlm_* flags, the fates_maxElements* dims, the
#     fates_hdim_* history maps, the hlm_current_* clock, numpft/nlevsclass/...).
#     We MUTATE those Refs here exactly where the Fortran assigns the bare module
#     variable. We do NOT redefine them.
#   * The EDParamsMod-side globals (dinc_vai, dlower_vai, maxpatch_total,
#     maxpatches_by_landuse, ED_val_history_*_bin_edges, regeneration_model,
#     bgc_soil_salinity, eca_plant_escalar, n/p_uptake_mode) live on the
#     `ed_params()` parameter-object instance — accessed/mutated via `ed_params()`.
#   * EDPftvarcon_inst -> `edpftvarcon_inst()`; prt_params is a module instance.
#   * num_elements/element_list/element_pos (PRTGenericMod) are Refs/arrays.
#   * fates_np_comp_scaling (FatesConstantsMod) is a Ref.
#   * bc_in/bc_out/bc_pconst allocation is delegated to the already-ported
#     allocate_bcin!/allocate_bcout!/allocate_bcpconst! in FatesInterfaceTypesMod
#     (sized identically to the Fortran allocate routines). Here we add the
#     interface-level orchestration: zero_bcs, set_bcs, set_bcpconst, and the
#     fates_interface_type container.
#
# ------------------------------------------------------------------------------
# fates_hist:  The Fortran does `use FatesHistoryInterfaceMod, only : fates_hist`
# (line 110) but NEVER references `fates_hist` anywhere in the body — the import
# is vestigial. `fates_history_maps` populates the `fates_hdim_*` module globals
# (in FatesInterfaceTypesMod), NOT `fates_hist`. So there is NOTHING to guard:
# the symbol `fates_hist` does not appear in this Julia port at all. The parent
# integrates FatesHistoryInterfaceMod (which DOES define fates_hist) earlier in
# the include order; this module simply doesn't need it.
#
# Standalone: nothing here is added to CLMInstances or any dual-copied struct.
# Deps: FatesInterfaceTypesMod (the bc/flag/dim Refs + allocators), EDParamsMod
# (ed_params + Register/Receive + Report), EDPftvarcon (edpftvarcon_inst +
# Report/Check), PRTParametersMod (prt_params), PRTGenericMod (element_list/pos/
# num_elements + organ/element ids + parteh hyp ids), PRTParamsFATESMod
# (PRT*Params), SFParamsMod (SpitFire*Params), FatesConstantsMod, FatesGlobals,
# FatesRadiationMemMod (num_swb/ivis/inir), FatesLitterMod (ncwd/ndcmpy),
# FatesHydraulicsMemMod (nshell/nlevsoi_hyd_max), FatesFuelClassesMod
# (num_fuel_classes), FatesRunningMeanMod, FatesDispersalMod, FatesUtilsMod
# (GetNeighborDistance), FatesPlantHydraulicsMod (InitHydroGlobals!),
# FatesParametersInterface (the reader + params type).

# for debugging this module (Fortran module var `debug`)
const _fates_interface_debug = Ref{Bool}(false)

# ===========================================================================
# fates_interface_type — the root of the FATES hierarchy of instantaneous
# state held by the HLM. Carries the site list and the per-site bc_in/bc_out
# vectors plus a single bc_pconst.
# ===========================================================================

"""
    fates_interface_type

Root container the host land model allocates per thread (Fortran
`fates_interface_type`). Holds the site list, the per-site boundary-condition
in/out vectors, and the (single) parameter-constant structure.

Fields match the Fortran type:
- `nsites::Int`
- `sites::Vector{ed_site_type}`     — FATES sites (Fortran pointer array)
- `bc_in::Vector{bc_in_type}`       — host -> FATES boundaries (per site)
- `bc_out::Vector{bc_out_type}`     — FATES -> host boundaries (per site)
- `bc_pconst::bc_pconst_type`       — parameter constants (one instance)
"""
Base.@kwdef mutable struct fates_interface_type <: AbstractFatesInterface
    nsites::Int = 0
    sites::Vector{ed_site_type}   = ed_site_type[]
    bc_in::Vector{bc_in_type}     = bc_in_type[]
    bc_out::Vector{bc_out_type}   = bc_out_type[]
    bc_pconst::bc_pconst_type     = bc_pconst_type()
    # FATES history-output interface (registry + buffers). Mirrors the Fortran
    # host's `fates_hist` (a `fates_history_interface_type` owned by clm_fates).
    # Instantiated + Init'd by clm_fates_init! when history is wired; the daily
    # (update_history_dyn1!) and per-timestep (update_history_hifrq1!) fills write
    # into its buffers. Write-only diagnostics — never read back into the sim.
    hist::Union{fates_history_interface_type,Nothing} = nothing
end

# ===========================================================================
# FatesInterfaceInit — initialize the FATES globals (log unit + verbosity).
# ===========================================================================

"""
    FatesInterfaceInit(log_unit, global_verbose)

Initialize the FATES global logging state. Mirrors Fortran
`FatesInterfaceInit` (delegates to `FatesGlobalsInit`).
"""
function FatesInterfaceInit(log_unit::Integer, global_verbose::Bool)
    FatesGlobalsInit(log_unit, global_verbose)
    return nothing
end

# ===========================================================================
# allocate_bcpconst / set_bcpconst — the parameter-constant structure.
# The sizing delegates to allocate_bcpconst! (FatesInterfaceTypesMod). set
# copies the ECA nutrient-competition params off EDPftvarcon_inst + the global
# eca_plant_escalar.
# ===========================================================================

"""
    allocate_bcpconst(bc_pconst, nlevdecomp)

Allocate the parameter-constant arrays (PFT-dimensioned ECA vectors + the
decomposition->uptake `j_uptake` map). Mirrors Fortran `allocate_bcpconst`.
Delegates to `allocate_bcpconst!`, sized by `numpft[]` and `nlevdecomp`.
"""
function allocate_bcpconst(bc_pconst::bc_pconst_type, nlevdecomp::Integer)
    allocate_bcpconst!(bc_pconst; num_pft = numpft[], nlevdecomp = nlevdecomp)
    return nothing
end

"""
    set_bcpconst(bc_pconst, nlevdecomp)

Fill the parameter-constant arrays from `EDPftvarcon_inst` and the global
`eca_plant_escalar`. Mirrors Fortran `set_bcpconst`. (`nlevdecomp` retained for
signature fidelity; not used in the Fortran body either.)
"""
function set_bcpconst(bc_pconst::bc_pconst_type, nlevdecomp::Integer)
    np  = numpft[]
    pft = edpftvarcon_inst()

    bc_pconst.vmax_nh4[1:np] .= pft.vmax_nh4[1:np]
    bc_pconst.vmax_no3[1:np] .= pft.vmax_no3[1:np]
    bc_pconst.vmax_p[1:np]   .= pft.vmax_p[1:np]

    bc_pconst.eca_km_nh4[1:np]      .= pft.eca_km_nh4[1:np]
    bc_pconst.eca_km_no3[1:np]      .= pft.eca_km_no3[1:np]
    bc_pconst.eca_km_p[1:np]        .= pft.eca_km_p[1:np]
    bc_pconst.eca_km_ptase[1:np]    .= pft.eca_km_ptase[1:np]
    bc_pconst.eca_vmax_ptase[1:np]  .= pft.eca_vmax_ptase[1:np]
    bc_pconst.eca_alpha_ptase[1:np] .= pft.eca_alpha_ptase[1:np]
    bc_pconst.eca_lambda_ptase[1:np] .= pft.eca_lambda_ptase[1:np]
    bc_pconst.eca_plant_escalar       = ed_params().eca_plant_escalar

    return nothing
end

# ===========================================================================
# zero_bcs — zero a site's input + output boundary conditions.
# ===========================================================================

"""
    zero_bcs(fates, s)

Zero the boundary conditions for site `s` of `fates`. Mirrors Fortran
`zero_bcs`. Conditional blocks (salinity, plant-hydro, parteh mode, LUH) are
gated on the same flags / globals the Fortran uses.
"""
function zero_bcs(fates::fates_interface_type, s::Integer)
    bc_in  = fates.bc_in[s]
    bc_out = fates.bc_out[s]

    # ---- Input boundaries ----
    fill!(bc_in.lightning24,   0.0)
    fill!(bc_in.pop_density,   0.0)
    fill!(bc_in.precip24_pa,   0.0)
    fill!(bc_in.relhumid24_pa, 0.0)
    fill!(bc_in.wind24_pa,     0.0)

    fill!(bc_in.solad_parb,     0.0)
    fill!(bc_in.solai_parb,     0.0)
    fill!(bc_in.smp_sl,         0.0)
    fill!(bc_in.eff_porosity_sl, 0.0)
    fill!(bc_in.watsat_sl,      0.0)
    fill!(bc_in.tempk_sl,       0.0)
    fill!(bc_in.h2o_liqvol_sl,  0.0)
    fill!(bc_in.filter_vegzen_pa, false)
    fill!(bc_in.coszen_pa,      0.0)
    fill!(bc_in.fcansno_pa,     0.0)
    fill!(bc_in.albgr_dir_rb,   0.0)
    fill!(bc_in.albgr_dif_rb,   0.0)
    bc_in.max_rooting_depth_index_col = 0
    bc_in.tot_het_resp   = 0.0
    bc_in.tot_somc       = 0.0
    bc_in.tot_litc       = 0.0
    bc_in.snow_depth_si  = 0.0
    bc_in.frac_sno_eff_si = 0.0
    fill!(bc_in.w_scalar_sisl, 0.0)
    fill!(bc_in.t_scalar_sisl, 0.0)

    if do_fates_salinity
        fill!(bc_in.salinity_sl, 0.0)
    end

    if hlm_use_planthydro[] == itrue
        fill!(bc_in.qflx_transp_pa, 0.0)
        fill!(bc_in.swrad_net_pa,   0.0)
        fill!(bc_in.lwrad_net_pa,   0.0)
        fill!(bc_in.watsat_sisl,    0.0)
        fill!(bc_in.watres_sisl,    0.0)
        fill!(bc_in.sucsat_sisl,    0.0)
        fill!(bc_in.bsw_sisl,       0.0)
        fill!(bc_in.hksat_sisl,     0.0)
    end

    # ---- Output boundaries ----
    fill!(bc_out.active_suction_sl, false)
    fill!(bc_out.fsun_pa,   0.0)
    fill!(bc_out.laisun_pa, 0.0)
    fill!(bc_out.laisha_pa, 0.0)
    fill!(bc_out.rootr_pasl, 0.0)
    fill!(bc_out.btran_pa,  0.0)

    # MIMICS litter quality, always initialize to unset
    bc_out.litt_flux_ligc_per_n = fates_unset_r8

    # Fates -> BGC fragmentation mass fluxes (gated on parteh hypothesis)
    pm = hlm_parteh_mode[]
    if pm == prt_carbon_allom_hyp
        fill!(bc_out.litt_flux_cel_c_si, 0.0)
        fill!(bc_out.litt_flux_lig_c_si, 0.0)
        fill!(bc_out.litt_flux_lab_c_si, 0.0)
        # Yes, zero out N flux arrays for c-only runs (CLM keeps them on+zero).
        fill!(bc_out.litt_flux_cel_n_si, 0.0)
        fill!(bc_out.litt_flux_lig_n_si, 0.0)
        fill!(bc_out.litt_flux_lab_n_si, 0.0)
    elseif pm == prt_cnp_flex_allom_hyp
        fill!(bc_in.plant_nh4_uptake_flux, 0.0)
        fill!(bc_in.plant_no3_uptake_flux, 0.0)
        fill!(bc_in.plant_p_uptake_flux,   0.0)
        fill!(bc_out.source_p,   0.0)
        fill!(bc_out.source_nh4, 0.0)
        fill!(bc_out.litt_flux_cel_c_si, 0.0)
        fill!(bc_out.litt_flux_lig_c_si, 0.0)
        fill!(bc_out.litt_flux_lab_c_si, 0.0)
        fill!(bc_out.litt_flux_cel_n_si, 0.0)
        fill!(bc_out.litt_flux_lig_n_si, 0.0)
        fill!(bc_out.litt_flux_lab_n_si, 0.0)
        fill!(bc_out.litt_flux_cel_p_si, 0.0)
        fill!(bc_out.litt_flux_lig_p_si, 0.0)
        fill!(bc_out.litt_flux_lab_p_si, 0.0)
    else
        println(fates_log(), "An unknown parteh hypothesis was passed")
        println(fates_log(), "while zeroing output boundary conditions")
        println(fates_log(), "hlm_parteh_mode: ", pm)
        fates_endrun("zero_bcs: unknown parteh hypothesis")
    end

    fill!(bc_out.rssun_pa, 0.0)
    fill!(bc_out.rssha_pa, 0.0)

    fill!(bc_out.albd_parb, 0.0)
    fill!(bc_out.albi_parb, 0.0)
    fill!(bc_out.fabd_parb, 0.0)
    fill!(bc_out.fabi_parb, 0.0)
    fill!(bc_out.ftdd_parb, 0.0)
    fill!(bc_out.ftid_parb, 0.0)
    fill!(bc_out.ftii_parb, 0.0)

    fill!(bc_out.elai_pa, 0.0)
    fill!(bc_out.esai_pa, 0.0)
    fill!(bc_out.tlai_pa, 0.0)
    fill!(bc_out.tsai_pa, 0.0)
    fill!(bc_out.htop_pa, 0.0)
    fill!(bc_out.hbot_pa, 0.0)
    fill!(bc_out.displa_pa, 0.0)
    fill!(bc_out.z0m_pa, 0.0)
    fill!(bc_out.dleaf_pa, 0.0)
    fill!(bc_out.nocomp_pft_label_pa, 0)

    fill!(bc_out.canopy_fraction_pa, 0.0)
    fill!(bc_out.frac_veg_nosno_alb_pa, 0.0)

    if hlm_use_planthydro[] == itrue
        fill!(bc_out.qflx_soil2root_sisl, 0.0)
        fill!(bc_out.qflx_ro_sisl,        0.0)
    end
    bc_out.plant_stored_h2o_si = 0.0

    # Land Use related
    bc_out.gpp_site = 0.0
    bc_out.ar_site  = 0.0
    bc_out.hrv_deadstemc_to_prod10c  = 0.0
    bc_out.hrv_deadstemc_to_prod100c = 0.0

    if hlm_use_luh[] == itrue
        fill!(bc_in.hlm_luh_states,      0.0)
        fill!(bc_in.hlm_luh_transitions, 0.0)
    end

    return nothing
end

# ===========================================================================
# allocate_bcin / allocate_bcout — delegate to the FatesInterfaceTypesMod
# allocators, threading the flag-driven dimension choices the Fortran makes.
# ===========================================================================

"""
    allocate_bcin(bc_in, nlevsoil_in, nlevdecomp_in, num_lu_harvest_cats,
                  num_luh2_states, num_luh2_transitions, surfpft_lb, surfpft_ub)

Allocate + initialize the host->FATES boundary-condition vectors. Mirrors
Fortran `allocate_bcin` — including the layer-consistency checks (vertsoilc
=> nlevdecomp == nlevsoil; otherwise nlevdecomp == 1) and the flag-driven
sizing of the nutrient-uptake / harvest / LUH / fixed-biogeog / SP arrays.
Delegates the bulk array sizing to `allocate_bcin!`.
"""
function allocate_bcin(bc_in::bc_in_type, nlevsoil_in::Integer, nlevdecomp_in::Integer,
                       num_lu_harvest_cats::Integer, num_luh2_states::Integer,
                       num_luh2_transitions::Integer, surfpft_lb::Integer, surfpft_ub::Integer)
    bc_in.nlevsoil = nlevsoil_in

    if nlevsoil_in > numlevsoil_max
        println(fates_log(), "The number of soil layers imposed by the host model")
        println(fates_log(), "is larger than what we have allocated in our static")
        println(fates_log(), "arrays. Please increase the size of numlevsoil_max")
        fates_endrun("allocate_bcin: nlevsoil_in > numlevsoil_max")
    end

    bc_in.nlevdecomp = nlevdecomp_in

    if hlm_use_vertsoilc[] == itrue
        if bc_in.nlevdecomp != bc_in.nlevsoil
            println(fates_log(), "The host has signaled a vertically resolved")
            println(fates_log(), "soil decomposition model. Therefore, the ")
            println(fates_log(), "total number of soil layers should equal the")
            println(fates_log(), "total number of decomposition layers.")
            println(fates_log(), "nlevdecomp: ", bc_in.nlevdecomp)
            println(fates_log(), "nlevsoil: ",   bc_in.nlevsoil)
            fates_endrun("allocate_bcin: vertsoilc nlevdecomp != nlevsoil")
        end
    else
        if bc_in.nlevdecomp != 1
            println(fates_log(), "The host has signaled a non-vertically resolved")
            println(fates_log(), "soil decomposition model. Therefore, the ")
            println(fates_log(), "total number of decomposition layers should be 1.")
            println(fates_log(), "nlevdecomp: ", bc_in.nlevdecomp)
            fates_endrun("allocate_bcin: non-vertsoilc nlevdecomp != 1")
        end
    end

    # Plant nutrient aquisition: if CNP (upscaled to competitors) pass back
    # competitor x 1, otherwise size 1x1.  (See Fortran note.)
    max_comp = (hlm_parteh_mode[] == prt_cnp_flex_allom_hyp) ? max_comp_per_site[] : 1

    # Fixed-biogeog uses different dims for LUH vs not.
    num_pft_alloc = 0
    num_lu_alloc  = 0
    if hlm_use_fixed_biogeog[] == itrue
        if hlm_use_luh[] == itrue
            # pft_areafrac_lu(size(hlm_pft_map,2), n_landuse_cats - n_crop_lu_types)
            num_pft_alloc = size(edpftvarcon_inst().hlm_pft_map, 2)
            num_lu_alloc  = n_landuse_cats - n_crop_lu_types
        else
            # pft_areafrac(surfpft_lb:surfpft_ub)
            num_pft_alloc = surfpft_ub - surfpft_lb + 1
        end
    end

    # SP-mode surface arrays use the surface-pft span.
    num_sp_alloc = (hlm_use_sp[] == itrue) ? (surfpft_ub - surfpft_lb + 1) : 0

    # Land-use harvest categories (Fortran allocates 0 when harvest is off).
    num_harvest = (hlm_use_lu_harvest[] > 0) ? num_lu_harvest_cats : 0

    # LUH state/transition arrays only when use_luh.
    nluh_states      = (hlm_use_luh[] == itrue) ? num_luh2_states : 0
    nluh_transitions = (hlm_use_luh[] == itrue) ? num_luh2_transitions : 0

    allocate_bcin!(bc_in;
                   npatches = maxpatch_total(), nlevsoil = nlevsoil_in,
                   nlevdecomp = nlevdecomp_in, num_rad_bands = num_swb,
                   max_comp = max_comp,
                   num_harvest_cats = num_harvest,
                   num_luh2_states = nluh_states,
                   num_luh2_transitions = nluh_transitions,
                   num_pft = num_pft_alloc, num_lu = num_lu_alloc,
                   num_sp = num_sp_alloc)

    return nothing
end

"""
    allocate_bcout(bc_out, nlevsoil_in, nlevdecomp_in)

Allocate + initialize the FATES->host boundary-condition vectors. Mirrors
Fortran `allocate_bcout`. Delegates the bulk array sizing to `allocate_bcout!`,
with the radiation arrays sized over `fates_maxPatchesPerSite` and the rest over
`maxpatch_total` as in the Fortran. The competitor dimension uses
`max_comp_per_site`.
"""
function allocate_bcout(bc_out::bc_out_type, nlevsoil_in::Integer, nlevdecomp_in::Integer)
    allocate_bcout!(bc_out;
                    npatches = maxpatch_total(), nlevsoil = nlevsoil_in,
                    nlevdecomp = nlevdecomp_in, num_rad_bands = num_swb,
                    max_comp = max_comp_per_site[])

    # Give the bare-ground root fractions a nominal fraction of unity over depth
    # (Fortran: bc_out%rootfr_pa(0,1:nlevsoil) = 1/nlevsoil). In the 1-based
    # Julia layout the bare-ground patch maps to row 1.
    if nlevsoil_in > 0 && size(bc_out.rootfr_pa, 1) >= 1
        bc_out.rootfr_pa[1, 1:nlevsoil_in] .= 1.0 / nlevsoil_in
    end

    return nothing
end

# ===========================================================================
# set_bcs — set boundary conditions not yet functional from the HLM. Currently
# just the (constant) soil salinity. MUST be called once at init, after the
# PFT parameter file is read, before history dims are set.
# ===========================================================================

"""
    set_bcs(bc_in)

Set the (time-invariant) boundary conditions not yet supplied by the HLM —
currently only the constant soil salinity. Mirrors Fortran `set_bcs`. Must run
once at initialization (not after).
"""
function set_bcs(bc_in::bc_in_type)
    if do_fates_salinity
        fill!(bc_in.salinity_sl, ed_params().bgc_soil_salinity)
    end
    return nothing
end

# ===========================================================================
# SetFatesGlobalElements1 — first FATES routine called: read the param file,
# and (for SP / nocomp) size maxpatch_total + fates_maxPatchesPerSite.
# ===========================================================================

"""
    SetFatesGlobalElements1(use_fates, surf_numpft, surf_numcft, param_reader)

First FATES global-dimension routine. Reads the FATES parameter file, then sets
`maxpatch_total` / `maxpatches_by_landuse` / `fates_maxPatchesPerSite[]` per the
SP / nocomp / fixed-biogeog flags. Mirrors Fortran `SetFatesGlobalElements1`.
The sp/biogeog/nocomp mode flags must already be set (via set_fates_ctrlparms).
"""
function SetFatesGlobalElements1(use_fates::Bool, surf_numpft::Integer,
                                 surf_numcft::Integer, param_reader)
    if use_fates
        # Self explanatory, read the fates parameter file
        FatesReadParameters(param_reader)

        fates_numpft = length(prt_params.wood_density)
        p = ed_params()

        if hlm_use_sp[] == itrue
            # For an SP run we use the primary patches to hold all PFTs — create
            # the same number of patches as the number of PFTs.
            p.maxpatches_by_landuse[primaryland] = fates_numpft
            p.maxpatches_by_landuse[secondaryland:n_landuse_cats] .= 0
            p.maxpatch_total = fates_numpft

            # SP runs need enough CLM/ELM-side patches to hold the LAI data — may
            # be larger than what fates requires. surf_numpft includes the bare
            # ground; maxpatch_total does not (so add 1).
            fates_maxPatchesPerSite[] = max(surf_numpft + surf_numcft, p.maxpatch_total + 1)
        else
            if hlm_use_nocomp[] == itrue
                p.maxpatches_by_landuse[primaryland] =
                    max(p.maxpatches_by_landuse[primaryland], fates_numpft)
                p.maxpatch_total = sum(p.maxpatches_by_landuse)
            end

            # maxpatch_total does not include the bare ground (so add 1)
            fates_maxPatchesPerSite[] = p.maxpatch_total + 1
        end
    end
    return nothing
end

# Convenience accessor: maxpatch_total lives on the ed_params() instance.
maxpatch_total() = ed_params().maxpatch_total

# ===========================================================================
# SetFatesGlobalElements2 — second FATES routine: derive numpft, leaf-age
# classes, history bin counts, VAI bins, uptake modes, max competitors, then
# init Hydro + PARTEH + history maps.
# ===========================================================================

"""
    SetFatesGlobalElements2(use_fates)

Second FATES global-dimension routine. Derives `numpft`, `nleafage`,
`nlevsclass`/`nlevage`/`nlevheight`/`nlevcoage`/`nlevdamage`,
`fates_maxElementsPerPatch`/`PerSite`, the N/P uptake modes, `max_comp_per_site`,
the VAI bin widths (dinc_vai/dlower_vai), the dispersal kernel mode, and runs
`InitHydroGlobals` / `InitPARTEHGlobals` / `fates_history_maps`. Mirrors Fortran
`SetFatesGlobalElements2`.
"""
function SetFatesGlobalElements2(use_fates::Bool)
    if use_fates
        p   = ed_params()
        pft = edpftvarcon_inst()

        # Number of FATES PFTs. The Julia prt_params.wood_density is a plain
        # 1-based Vector (no lbound==0 case), so numpft == length directly.
        # (Fortran handles lbound 0 or 1; here it is always 1-based.)
        numpft[] = length(prt_params.wood_density)

        if numpft[] > maxpft
            println(fates_log(), "The number of PFTs dictated by the FATES parameter file")
            println(fates_log(), "is larger than the maximum allowed. Increase maxpft")
            fates_endrun("SetFatesGlobalElements2: numpft > maxpft")
        end

        # Number of leaf age-classes (2nd dim of leaf_long).
        nleafage[] = size(prt_params.leaf_long, 2)

        nlevsclass[] = length(p.ED_val_history_sizeclass_bin_edges)

        # Restart-file / cohort-array memory sizing.
        if hlm_use_sp[] == itrue
            fates_maxElementsPerPatch[] = num_swb
        else
            fates_maxElementsPerPatch[] = max(num_swb, p.max_cohort_per_patch,
                                              ndcmpy * hlm_maxlevsoil[],
                                              ncwd   * hlm_maxlevsoil[])
        end

        fates_maxElementsPerSite[] = max(
            fates_maxPatchesPerSite[] * fates_maxElementsPerPatch[],
            numWaterMem * numpft[], num_vegtemp_mem, num_elements[],
            nlevsclass[] * numpft[] * n_term_mort_types)

        if hlm_use_planthydro[] == itrue
            fates_maxElementsPerSite[] =
                max(fates_maxElementsPerSite[], nshell * nlevsoi_hyd_max)
        end

        # N/P uptake mode (prescribed if any prescribed_*uptake nonzero).
        if any(x -> abs(x) > nearzero, pft.prescribed_nuptake)
            p.n_uptake_mode = prescribed_n_uptake
        else
            p.n_uptake_mode = coupled_n_uptake
        end
        if any(x -> abs(x) > nearzero, pft.prescribed_puptake)
            p.p_uptake_mode = prescribed_p_uptake
        else
            p.p_uptake_mode = coupled_p_uptake
        end

        # Max nutrient-acquisition competitors per site.
        if hlm_parteh_mode[] == prt_cnp_flex_allom_hyp
            if (p.p_uptake_mode == coupled_p_uptake) || (p.n_uptake_mode == coupled_n_uptake)
                max_comp_per_site[]   = fates_maxElementsPerSite[]
                fates_np_comp_scaling[] = coupled_np_comp_scaling
            else
                max_comp_per_site[]   = 1
                fates_np_comp_scaling[] = trivial_np_comp_scaling
            end
        else
            max_comp_per_site[]   = 1
            fates_np_comp_scaling[] = trivial_np_comp_scaling
        end

        # VAI bin widths for radiative transfer.
        for i in 1:nlevleaf
            p.dinc_vai[i] = p.ED_val_vai_top_bin_width *
                            p.ED_val_vai_width_increase_factor^(i - 1)
        end
        # lower edges of VAI bins
        for i in 1:nlevleaf
            p.dlower_vai[i] = sum(@view p.dinc_vai[1:i])
        end

        # History size/age/height/coage/damage class bin counts (1-indexed).
        nlevage[]    = length(p.ED_val_history_ageclass_bin_edges)
        nlevheight[] = length(p.ED_val_history_height_bin_edges)
        nlevcoage[]  = length(p.ED_val_history_coageclass_bin_edges)
        nlevdamage[] = length(p.ED_val_history_damage_bin_edges)

        # Sanity checks on the bin edge arrays (start at zero, monotone increasing).
        if p.ED_val_history_sizeclass_bin_edges[1] != 0.0
            println(fates_log(), "size class bins specified in parameter file must start at zero")
            fates_endrun("SetFatesGlobalElements2: size class bins must start at zero")
        end
        if p.ED_val_history_ageclass_bin_edges[1] != 0.0
            println(fates_log(), "age class bins specified in parameter file must start at zero")
            fates_endrun("SetFatesGlobalElements2: age class bins must start at zero")
        end
        if p.ED_val_history_height_bin_edges[1] != 0.0
            println(fates_log(), "height class bins specified in parameter file must start at zero")
            fates_endrun("SetFatesGlobalElements2: height class bins must start at zero")
        end
        for i in 2:nlevsclass[]
            if (p.ED_val_history_sizeclass_bin_edges[i] -
                p.ED_val_history_sizeclass_bin_edges[i-1]) <= 0.0
                println(fates_log(), "size class bins must be monotonically increasing")
                fates_endrun("SetFatesGlobalElements2: size class bins not monotone")
            end
        end
        for i in 2:nlevage[]
            if (p.ED_val_history_ageclass_bin_edges[i] -
                p.ED_val_history_ageclass_bin_edges[i-1]) <= 0.0
                println(fates_log(), "age class bins must be monotonically increasing")
                fates_endrun("SetFatesGlobalElements2: age class bins not monotone")
            end
        end
        for i in 2:nlevheight[]
            if (p.ED_val_history_height_bin_edges[i] -
                p.ED_val_history_height_bin_edges[i-1]) <= 0.0
                println(fates_log(), "height class bins must be monotonically increasing")
                fates_endrun("SetFatesGlobalElements2: height class bins not monotone")
            end
        end
        for i in 2:nlevcoage[]
            if (p.ED_val_history_coageclass_bin_edges[i] -
                p.ED_val_history_coageclass_bin_edges[i-1]) <= 0.0
                println(fates_log(), "cohort age class bins must be monotonically increasing")
                fates_endrun("SetFatesGlobalElements2: cohort age class bins not monotone")
            end
        end

        # Set the dispersal kernel mode if any seed-dispersal params are set.
        if any(x -> x < fates_check_param_set, pft.seed_dispersal_pdf_scale)
            fates_dispersal_kernel_mode[] = fates_dispersal_kernel_exponential
        end

        # Init Hydro globals (water retention functions) — needs numpft.
        InitHydroGlobals!(numpft[])

        # Init PARTEH globals + element list/mapping.
        InitPARTEHGlobals()

        # Set the history-output mapping arrays.
        fates_history_maps()
    else
        # FATES off: keep the cohort dim minimal.
        fates_maxElementsPerPatch[] = 1
        fates_maxElementsPerSite[]  = 1
    end
    return nothing
end

# ===========================================================================
# InitTimeAveragingGlobals — define the running-mean time-averaging "definition"
# globals. The Fortran allocates pointer instances of the rmean class and calls
# %define on each; in Julia these definitions live in FatesRunningMeanMod /
# FatesPatchMod as `rmean_def_type` instances. We hold them here as the module
# globals the Fortran names (ema_24hr, fixed_24hr, ...) and `define!` them with
# the host stepsize.
# ===========================================================================

# Module-global running-mean *definitions* (Fortran pointer instances). They are
# `define!`d in InitTimeAveragingGlobals once hlm_stepsize is known.
const ema_24hr            = rmean_def_type()
const fixed_24hr          = rmean_def_type()
const ema_lpa             = rmean_def_type()
const ema_sdlng_emerg_h2o = rmean_def_type()
const ema_sdlng_mort_par  = rmean_def_type()
const ema_sdlng2sap_par   = rmean_def_type()
const ema_sdlng_mdd       = rmean_def_type()
const ema_longterm        = rmean_def_type()

"""
    InitTimeAveragingGlobals()

Instantiate (define) the running-mean time-averaging definition globals from the
HLM step-size + the FATES timescale parameters. Mirrors Fortran
`InitTimeAveragingGlobals`. NOTE the HLM step-size is treated as a global
constant (it would not be if HLM timesteps became dynamic).
"""
function InitTimeAveragingGlobals()
    p  = ed_params()
    dt = hlm_stepsize[]

    define!(ema_24hr,   sec_per_day, dt, moving_ema_window)
    define!(fixed_24hr, sec_per_day, dt, fixed_window)
    # ema_lpa: parameter has units of days
    define!(ema_lpa, p.photo_temp_acclim_timescale * sec_per_day, dt, moving_ema_window)
    define!(ema_sdlng_emerg_h2o, p.sdlng_emerg_h2o_timescale * sec_per_day, dt, moving_ema_window)
    define!(ema_sdlng_mort_par,  p.sdlng_mort_par_timescale  * sec_per_day, dt, moving_ema_window)
    define!(ema_sdlng2sap_par,   p.sdlng2sap_par_timescale   * sec_per_day, dt, moving_ema_window)
    define!(ema_sdlng_mdd,       p.sdlng_mdd_timescale       * sec_per_day, dt, moving_ema_window)
    # ema_longterm: parameter has units of years
    define!(ema_longterm, p.photo_temp_acclim_thome_time * days_per_year * sec_per_day,
            dt, moving_ema_window)
    return nothing
end

# ===========================================================================
# InitPARTEHGlobals — associate PARTEH elements w/ the FATES element list and
# init the relevant PRT global mapping tables.
# ===========================================================================

"""
    InitPARTEHGlobals()

Initialize the Plant Allocation and Reactive Transport (PARTEH) globals: set
`num_elements` / `element_list` / `element_pos` per the parteh mode and call the
matching PRT global-init routine. Mirrors Fortran `InitPARTEHGlobals`.
"""
function InitPARTEHGlobals()
    pm = hlm_parteh_mode[]
    if pm == prt_carbon_allom_hyp
        num_elements[] = 1
        resize!(element_list, num_elements[])
        element_list[1] = carbon12_element
        fill!(element_pos, 0)
        element_pos[carbon12_element] = 1
        InitPRTGlobalAllometricCarbon!()
    elseif pm == prt_cnp_flex_allom_hyp
        num_elements[] = 3
        resize!(element_list, num_elements[])
        element_list[1] = carbon12_element
        element_list[2] = nitrogen_element
        element_list[3] = phosphorus_element
        fill!(element_pos, 0)
        element_pos[carbon12_element]   = 1
        element_pos[nitrogen_element]   = 2
        element_pos[phosphorus_element] = 3
        InitPRTGlobalAllometricCNP()
    else
        println(fates_log(), "You specified an unknown PRT module")
        println(fates_log(), "Check your setting for fates_parteh_mode")
        fates_endrun("InitPARTEHGlobals: unknown parteh mode")
    end
    return nothing
end

# ===========================================================================
# fates_history_maps — allocate + populate the fates_hdim_* history dimension
# maps (multiplexed-index -> component-dimension lookups). The fates_hdim_*
# Refs live in FatesInterfaceTypesMod; here we size + fill them.
# ===========================================================================

"""
    fates_history_maps()

Allocate and populate the `fates_hdim_*` history-output dimension mapping arrays
(the maps from multiplexed history indices like "scpf" back to their component
dimensions size-class/pft/etc.). Mirrors Fortran `fates_history_maps`. Operates
on the FatesInterfaceTypesMod module-global Refs + the ed_params() bin edges.
"""
function fates_history_maps()
    p   = ed_params()
    nsc = nlevsclass[]
    npf = numpft[]
    nag = nlevage[]
    nht = nlevheight[]
    ncoage = nlevcoage[]
    ndam = nlevdamage[]
    nel  = num_elements[]

    # Single-dimension definitions.
    fates_hdim_levsclass[]   = zeros(Float64, nsc)
    fates_hdim_pfmap_levscpf[] = zeros(Int, nsc * npf)
    fates_hdim_scmap_levscpf[] = zeros(Int, nsc * npf)
    fates_hdim_levpft[]      = zeros(Int, npf)
    fates_hdim_levlanduse[]  = zeros(Int, n_landuse_cats)
    fates_hdim_levfuel[]     = zeros(Int, num_fuel_classes)
    fates_hdim_levcwdsc[]    = zeros(Int, ncwd)
    fates_hdim_levage[]      = zeros(Float64, nag)
    fates_hdim_levheight[]   = zeros(Float64, nht)
    fates_hdim_levcoage[]    = zeros(Float64, ncoage)
    fates_hdim_pfmap_levcapf[] = zeros(Int, ncoage * npf)
    fates_hdim_camap_levcapf[] = zeros(Int, ncoage * npf)

    fates_hdim_levcan[]   = zeros(Int, nclmax)
    fates_hdim_levelem[]  = zeros(Int, nel)
    fates_hdim_levleaf[]  = zeros(Float64, nlevleaf)
    fates_hdim_canmap_levcnlf[] = zeros(Int, nlevleaf * nclmax)
    fates_hdim_lfmap_levcnlf[]  = zeros(Int, nlevleaf * nclmax)
    fates_hdim_canmap_levcnlfpf[] = zeros(Int, nlevleaf * nclmax * npf)
    fates_hdim_lfmap_levcnlfpf[]  = zeros(Int, nlevleaf * nclmax * npf)
    fates_hdim_pftmap_levcnlfpf[] = zeros(Int, nlevleaf * nclmax * npf)
    fates_hdim_scmap_levscag[] = zeros(Int, nsc * nag)
    fates_hdim_agmap_levscag[] = zeros(Int, nsc * nag)
    fates_hdim_scmap_levscagpft[] = zeros(Int, nsc * nag * npf)
    fates_hdim_agmap_levscagpft[] = zeros(Int, nsc * nag * npf)
    fates_hdim_pftmap_levscagpft[] = zeros(Int, nsc * nag * npf)
    fates_hdim_agmap_levagepft[] = zeros(Int, nag * npf)
    fates_hdim_pftmap_levagepft[] = zeros(Int, nag * npf)
    fates_hdim_agmap_levagefuel[] = zeros(Int, nag * num_fuel_classes)
    fates_hdim_fscmap_levagefuel[] = zeros(Int, nag * num_fuel_classes)

    fates_hdim_elmap_levelpft[] = zeros(Int, nel * npf)
    fates_hdim_elmap_levelcwd[] = zeros(Int, nel * ncwd)
    fates_hdim_elmap_levelage[] = zeros(Int, nel * nag)
    fates_hdim_pftmap_levelpft[] = zeros(Int, nel * npf)
    fates_hdim_cwdmap_levelcwd[] = zeros(Int, nel * ncwd)
    fates_hdim_agemap_levelage[] = zeros(Int, nel * nag)

    fates_hdim_levdamage[]   = zeros(Float64, ndam)
    fates_hdim_scmap_levcdsc[] = zeros(Int, nsc * ndam)
    fates_hdim_cdmap_levcdsc[] = zeros(Int, nsc * ndam)
    fates_hdim_scmap_levcdpf[] = zeros(Int, nsc * ndam * npf)
    fates_hdim_cdmap_levcdpf[] = zeros(Int, nsc * ndam * npf)
    fates_hdim_pftmap_levcdpf[] = zeros(Int, nsc * ndam * npf)

    # Fill the dimension value arrays directly from the bin edges.
    fates_hdim_levsclass[] .= p.ED_val_history_sizeclass_bin_edges
    fates_hdim_levage[]    .= p.ED_val_history_ageclass_bin_edges
    fates_hdim_levheight[] .= p.ED_val_history_height_bin_edges
    fates_hdim_levcoage[]  .= p.ED_val_history_coageclass_bin_edges
    fates_hdim_levleaf[]   .= p.dlower_vai
    fates_hdim_levdamage[] .= p.ED_val_history_damage_bin_edges

    # Simple index arrays.
    for ipft in 1:npf
        fates_hdim_levpft[][ipft] = ipft
    end
    for ifuel in 1:num_fuel_classes
        fates_hdim_levfuel[][ifuel] = ifuel
    end
    for icwd in 1:ncwd
        fates_hdim_levcwdsc[][icwd] = icwd
    end
    for ican in 1:nclmax
        fates_hdim_levcan[][ican] = ican
    end
    for ilu in 1:n_landuse_cats
        fates_hdim_levlanduse[][ilu] = ilu
    end
    for iel in 1:nel
        fates_hdim_levelem[][iel] = element_list[iel]
    end

    # element x pft
    let i = 0
        for iel in 1:nel, ipft in 1:npf
            i += 1
            fates_hdim_elmap_levelpft[][i]  = iel
            fates_hdim_pftmap_levelpft[][i] = ipft
        end
    end
    # element x cwd
    let i = 0
        for iel in 1:nel, icwd in 1:ncwd
            i += 1
            fates_hdim_elmap_levelcwd[][i]  = iel
            fates_hdim_cwdmap_levelcwd[][i] = icwd
        end
    end
    # element x age
    let i = 0
        for iel in 1:nel, iage in 1:nag
            i += 1
            fates_hdim_elmap_levelage[][i]  = iel
            fates_hdim_agemap_levelage[][i] = iage
        end
    end
    # size-class x pft
    let i = 0
        for ipft in 1:npf, isc in 1:nsc
            i += 1
            fates_hdim_pfmap_levscpf[][i] = ipft
            fates_hdim_scmap_levscpf[][i] = isc
        end
    end
    # cohort-age x pft
    let i = 0
        for ipft in 1:npf, icoage in 1:ncoage
            i += 1
            fates_hdim_pfmap_levcapf[][i] = ipft
            fates_hdim_camap_levcapf[][i] = icoage
        end
    end
    # canopy-layer x leaf-layer
    let i = 0
        for ican in 1:nclmax, ileaf in 1:nlevleaf
            i += 1
            fates_hdim_canmap_levcnlf[][i] = ican
            fates_hdim_lfmap_levcnlf[][i]  = ileaf
        end
    end
    # size-class x age
    let i = 0
        for iage in 1:nag, isc in 1:nsc
            i += 1
            fates_hdim_scmap_levscag[][i] = isc
            fates_hdim_agmap_levscag[][i] = iage
        end
    end
    # size-class x crowndamage
    let i = 0
        for icdam in 1:ndam, isc in 1:nsc
            i += 1
            fates_hdim_scmap_levcdsc[][i] = isc
            fates_hdim_cdmap_levcdsc[][i] = icdam
        end
    end
    # size x crowndamage x pft
    let i = 0
        for ipft in 1:npf, icdam in 1:ndam, isc in 1:nsc
            i += 1
            fates_hdim_scmap_levcdpf[][i]  = isc
            fates_hdim_cdmap_levcdpf[][i]  = icdam
            fates_hdim_pftmap_levcdpf[][i] = ipft
        end
    end
    # can-layer x leaf-layer x pft
    let i = 0
        for ipft in 1:npf, ican in 1:nclmax, ileaf in 1:nlevleaf
            i += 1
            fates_hdim_canmap_levcnlfpf[][i] = ican
            fates_hdim_lfmap_levcnlfpf[][i]  = ileaf
            fates_hdim_pftmap_levcnlfpf[][i] = ipft
        end
    end
    # size-class x age x pft
    let i = 0
        for ipft in 1:npf, iage in 1:nag, isc in 1:nsc
            i += 1
            fates_hdim_scmap_levscagpft[][i]  = isc
            fates_hdim_agmap_levscagpft[][i]  = iage
            fates_hdim_pftmap_levscagpft[][i] = ipft
        end
    end
    # age x pft
    let i = 0
        for ipft in 1:npf, iage in 1:nag
            i += 1
            fates_hdim_agmap_levagepft[][i]  = iage
            fates_hdim_pftmap_levagepft[][i] = ipft
        end
    end
    # age x fuel
    let i = 0
        for iage in 1:nag, ifuel in 1:num_fuel_classes
            i += 1
            fates_hdim_agmap_levagefuel[][i]  = iage
            fates_hdim_fscmap_levagefuel[][i] = ifuel
        end
    end

    return nothing
end

# ===========================================================================
# SetFatesTime — set the FATES clock from the HLM. All sites synchronous.
# ===========================================================================

"""
    SetFatesTime(current_year_in, current_month_in, current_day_in, current_tod_in,
                 current_date_in, reference_date_in, model_day_in, day_of_year_in,
                 days_per_year_in, freq_day_in)

Set the FATES clock module-globals from the HLM. Mirrors Fortran `SetFatesTime`.
Called directly from the HLM each step.
"""
function SetFatesTime(current_year_in::Integer, current_month_in::Integer,
                      current_day_in::Integer, current_tod_in::Integer,
                      current_date_in::Integer, reference_date_in::Integer,
                      model_day_in::Real, day_of_year_in::Integer,
                      days_per_year_in::Integer, freq_day_in::Real)
    hlm_current_year[]   = current_year_in
    hlm_current_month[]  = current_month_in
    hlm_current_day[]    = current_day_in
    hlm_current_tod[]    = current_tod_in
    hlm_current_date[]   = current_date_in
    hlm_reference_date[] = reference_date_in
    hlm_model_day[]      = model_day_in
    hlm_day_of_year[]    = day_of_year_in
    hlm_days_per_year[]  = days_per_year_in
    hlm_freq_day[]       = freq_day_in
    return nothing
end

# ===========================================================================
# set_fates_ctrlparms — transfer HLM control params/dims into FATES. String-tag
# dispatch matching the Fortran. 'flush_to_unset' / 'check_allset' are the
# bookends; the value tags set the hlm_* Refs.
# ===========================================================================

# Local sentinels (match the Fortran parameters in set_fates_ctrlparms).
const _ctrl_unset_int    = -999
const _ctrl_unset_double = -999.9

"""
    set_fates_ctrlparms(tag; ival=nothing, rval=nothing, cval=nothing)

Transfer a single HLM control parameter / dimension into the FATES module-global
`hlm_*` Refs (FatesInterfaceTypesMod). Mirrors Fortran `set_fates_ctrlparms`.
`tag` selects which parameter; the value comes via `ival` (Int), `rval`
(Float64), or `cval` (String). The special tags `"flush_to_unset"` and
`"check_allset"` bracket a transfer sequence (flush all to the unset sentinel,
then validate everything was set).
"""
function set_fates_ctrlparms(tag::AbstractString; ival::Union{Integer,Nothing}=nothing,
                             rval::Union{Real,Nothing}=nothing,
                             cval::Union{AbstractString,Nothing}=nothing)
    t = strip(tag)

    if t == "flush_to_unset"
        if fates_global_verbose()
            println(fates_log(), "Flushing FATES control parameters prior to transfer from host")
        end
        hlm_numSWb[]     = _ctrl_unset_int
        hlm_inir[]       = _ctrl_unset_int
        hlm_ivis[]       = _ctrl_unset_int
        hlm_is_restart[] = _ctrl_unset_int
        hlm_maxlevsoil[] = _ctrl_unset_int
        hlm_name[]       = "unset"
        hlm_hio_ignore_val[] = _ctrl_unset_double
        hlm_masterproc[] = _ctrl_unset_int
        hlm_ipedof[]     = _ctrl_unset_int
        hlm_nu_com[]     = "unset"
        hlm_decomp[]     = "unset"
        hlm_nitrogen_spec[]   = _ctrl_unset_int
        hlm_use_tree_damage[] = _ctrl_unset_int
        hlm_phosphorus_spec[] = _ctrl_unset_int
        hlm_use_ch4[]       = _ctrl_unset_int
        hlm_use_vertsoilc[] = _ctrl_unset_int
        hlm_parteh_mode[]   = _ctrl_unset_int
        hlm_spitfire_mode[] = _ctrl_unset_int
        hlm_seeddisp_cadence[] = _ctrl_unset_int
        hlm_sf_nofire_def[] = _ctrl_unset_int
        hlm_sf_scalar_lightning_def[]     = _ctrl_unset_int
        hlm_sf_successful_ignitions_def[] = _ctrl_unset_int
        hlm_sf_anthro_ignitions_def[]     = _ctrl_unset_int
        hlm_use_planthydro[]  = _ctrl_unset_int
        hlm_use_lu_harvest[]  = _ctrl_unset_int
        hlm_num_lu_harvest_cats[]  = _ctrl_unset_int
        hlm_num_luh2_states[]      = _ctrl_unset_int
        hlm_num_luh2_transitions[] = _ctrl_unset_int
        hlm_use_cohort_age_tracking[] = _ctrl_unset_int
        hlm_use_logging[]  = _ctrl_unset_int
        hlm_use_ed_st3[]   = _ctrl_unset_int
        hlm_use_ed_prescribed_phys[] = _ctrl_unset_int
        hlm_use_fixed_biogeog[] = _ctrl_unset_int
        hlm_use_nocomp[] = _ctrl_unset_int
        hlm_use_sp[]     = _ctrl_unset_int
        hlm_use_inventory_init[]  = _ctrl_unset_int
        hlm_inventory_ctrl_file[] = "unset"
        hlm_hist_level_dynam[] = _ctrl_unset_int
        hlm_hist_level_hifrq[] = _ctrl_unset_int

    elseif t == "check_allset"
        _check_allset()

    else
        # value tags — dispatch by which of ival/rval/cval is present.
        if ival !== nothing
            _set_ctrlparm_ival(t, ival)
        end
        if rval !== nothing
            if t == "hio_ignore_val"
                hlm_hio_ignore_val[] = rval
                fates_global_verbose() &&
                    println(fates_log(), "Transfering hio_ignore_val = ", rval, " to FATES")
            else
                println(fates_log(), "fates NL tag not recognized:", t)
            end
        end
        if cval !== nothing
            if t == "hlm_name"
                hlm_name[] = strip(cval)
                fates_global_verbose() &&
                    println(fates_log(), "Transfering the HLM name = ", strip(cval))
            elseif t == "nu_com"
                hlm_nu_com[] = strip(cval)
                fates_global_verbose() &&
                    println(fates_log(), "Transfering the nutrient competition name = ", strip(cval))
            elseif t == "decomp_method"
                hlm_decomp[] = strip(cval)
                fates_global_verbose() &&
                    println(fates_log(), "Transfering the decomp method name = ", strip(cval))
            elseif t == "inventory_ctrl_file"
                hlm_inventory_ctrl_file[] = strip(cval)
                fates_global_verbose() &&
                    println(fates_log(), "Transfering the inventory control file = ", strip(cval))
            else
                println(fates_log(), "fates NL tag not recognized:", t)
            end
        end
    end

    return nothing
end

# Integer-valued control-parameter tag dispatch (Fortran's ival select case).
function _set_ctrlparm_ival(t::AbstractString, ival::Integer)
    if     t == "masterproc";        hlm_masterproc[] = ival
    elseif t == "num_sw_bbands";     hlm_numSWb[] = ival
    elseif t == "vis_sw_index";      hlm_ivis[] = ival
    elseif t == "nir_sw_index";      hlm_inir[] = ival
    elseif t == "is_restart";        hlm_is_restart[] = ival
    elseif t == "num_lev_soil";      hlm_maxlevsoil[] = ival
    elseif t == "soilwater_ipedof";  hlm_ipedof[] = ival
    elseif t == "use_tree_damage";   hlm_use_tree_damage[] = ival
    elseif t == "nitrogen_spec";     hlm_nitrogen_spec[] = ival
    elseif t == "phosphorus_spec";   hlm_phosphorus_spec[] = ival
    elseif t == "use_ch4";           hlm_use_ch4[] = ival
    elseif t == "use_vertsoilc";     hlm_use_vertsoilc[] = ival
    elseif t == "parteh_mode";       hlm_parteh_mode[] = ival
    elseif t == "seeddisp_cadence";  hlm_seeddisp_cadence[] = ival
    elseif t == "spitfire_mode";     hlm_spitfire_mode[] = ival
    elseif t == "sf_nofire_def";             hlm_sf_nofire_def[] = ival
    elseif t == "sf_scalar_lightning_def";   hlm_sf_scalar_lightning_def[] = ival
    elseif t == "sf_successful_ignitions_def"; hlm_sf_successful_ignitions_def[] = ival
    elseif t == "sf_anthro_ignitions_def";   hlm_sf_anthro_ignitions_def[] = ival
    elseif t == "use_fixed_biogeog"; hlm_use_fixed_biogeog[] = ival
    elseif t == "use_nocomp";        hlm_use_nocomp[] = ival
    elseif t == "use_sp";            hlm_use_sp[] = ival
    elseif t == "use_planthydro";    hlm_use_planthydro[] = ival
    elseif t == "use_lu_harvest";    hlm_use_lu_harvest[] = ival
    elseif t == "num_lu_harvest_cats"; hlm_num_lu_harvest_cats[] = ival
    elseif t == "use_luh2";          hlm_use_luh[] = ival
    elseif t == "use_fates_potentialveg"; hlm_use_potentialveg[] = ival
    elseif t == "num_luh2_states";        hlm_num_luh2_states[] = ival
    elseif t == "num_luh2_transitions";   hlm_num_luh2_transitions[] = ival
    elseif t == "use_cohort_age_tracking"; hlm_use_cohort_age_tracking[] = ival
    elseif t == "use_logging";       hlm_use_logging[] = ival
    elseif t == "use_ed_st3";        hlm_use_ed_st3[] = ival
    elseif t == "use_ed_prescribed_phys"; hlm_use_ed_prescribed_phys[] = ival
    elseif t == "use_inventory_init"; hlm_use_inventory_init[] = ival
    elseif t == "hist_hifrq_dimlevel"; hlm_hist_level_hifrq[] = ival
    elseif t == "hist_dynam_dimlevel"; hlm_hist_level_dynam[] = ival
    else
        println(fates_log(), "fates NL tag not recognized:", t)
        return nothing
    end
    if fates_global_verbose()
        println(fates_log(), "Transfering ", t, " = ", ival, " to FATES")
    end
    return nothing
end

# The 'check_allset' validation block. Errors (fates_endrun) on any unset /
# inconsistent parameter, exactly as the Fortran.
function _check_allset()
    pft = edpftvarcon_inst()

    hlm_numSWb[] == _ctrl_unset_int &&
        fates_endrun("FATES dimension/parameter unset: num_sw_rad_bbands")
    hlm_masterproc[] == _ctrl_unset_int &&
        fates_endrun("FATES parameter unset: hlm_masterproc")
    hlm_numSWb[] != num_swb &&
        fates_endrun("FATES requires 2 shortwave bands (VIS/NIR); HLM signaled otherwise")

    if !(hlm_use_planthydro[] == 1 || hlm_use_planthydro[] == 0)
        fates_endrun("The FATES namelist planthydro flag must be 0 or 1")
    elseif hlm_use_planthydro[] == 1
        if fates_global_verbose()
            println(fates_log(), "use_fates_planthydro is an EXPERIMENTAL FEATURE")
        end
    end

    (hlm_use_lu_harvest[] < 0 || hlm_use_lu_harvest[] > 1) &&
        fates_endrun("The FATES lu_harvest flag must be 0 or 1")
    hlm_num_lu_harvest_cats[] < 0 &&
        fates_endrun("The FATES number of hlm harvest cats must be >= 0")
    hlm_num_luh2_states[] < 0 &&
        fates_endrun("The FATES number of hlm luh state cats must be >= 0")
    hlm_num_luh2_transitions[] < 0 &&
        fates_endrun("The FATES number of hlm luh state transition cats must be >= 0")
    !(hlm_use_logging[] == 1 || hlm_use_logging[] == 0) &&
        fates_endrun("The FATES namelist use_logging flag must be 0 or 1")

    if any(x -> x < fates_check_param_set, pft.mort_ip_age_senescence) &&
       hlm_use_cohort_age_tracking[] == 0
        fates_endrun("Age dependent mortality cannot be on if cohort age tracking is off")
    end

    !(hlm_use_ed_st3[] == 1 || hlm_use_ed_st3[] == 0) &&
        fates_endrun("The FATES namelist stand structure flag must be 0 or 1")
    !(hlm_use_ed_prescribed_phys[] == 1 || hlm_use_ed_prescribed_phys[] == 0) &&
        fates_endrun("The FATES namelist prescribed physiology flag must be 0 or 1")
    (hlm_use_ed_prescribed_phys[] == 1 && hlm_use_ed_st3[] == 1) &&
        fates_endrun("FATES ST3 and prescribed physiology cannot both be on")
    (hlm_use_inventory_init[] == 1 && hlm_use_cohort_age_tracking[] == 1) &&
        fates_endrun("Fates inventory init cannot be used with age dependent mortality")
    !(hlm_use_inventory_init[] == 1 || hlm_use_inventory_init[] == 0) &&
        fates_endrun("The FATES NL inventory flag must be 0 or 1")
    strip(hlm_inventory_ctrl_file[]) == "unset" &&
        fates_endrun("namelist entry for fates inventory control file is unset")

    hlm_ivis[] != ivis &&
        fates_endrun("FATES VIS shortwave index differs from the HLM")
    hlm_inir[] != inir &&
        fates_endrun("FATES NIR shortwave index differs from the HLM")
    hlm_is_restart[] == _ctrl_unset_int &&
        fates_endrun("FATES parameter unset: hlm_is_restart")
    hlm_maxlevsoil[] == _ctrl_unset_int &&
        fates_endrun("FATES dimension/parameter unset: hlm_maxlevsoil")
    strip(hlm_name[]) == "unset" &&
        fates_endrun("FATES dimension/parameter unset: hlm_name")

    strip(hlm_decomp[]) == "unset" &&
        fates_endrun("FATES dimension/parameter unset: hlm_decomp (valid: MIMICS, CENTURY, CTC)")
    if !(strip(hlm_decomp[]) in ("MIMICS", "CENTURY", "CTC", "NONE"))
        fates_endrun("FATES hlm_decomp invalid (valid: NONE, MIMICS, CENTURY, CTC): " * hlm_decomp[])
    end

    strip(hlm_nu_com[]) == "unset" &&
        fates_endrun("FATES dimension/parameter unset: hlm_nu_com")

    if hlm_use_tree_damage[] == _ctrl_unset_int
        fates_endrun("FATES dimension/parameter unset: hlm_use_tree_damage")
    else
        (hlm_use_tree_damage[] == itrue && hlm_parteh_mode[] == prt_cnp_flex_allom_hyp) &&
            fates_endrun("FATES tree damage is not (yet) compatible with CNP allocation")
    end

    hlm_nitrogen_spec[] == _ctrl_unset_int &&
        fates_endrun("FATES parameters unset: hlm_nitrogen_spec")
    hlm_phosphorus_spec[] == _ctrl_unset_int &&
        fates_endrun("FATES parameters unset: hlm_phosphorus_spec")
    abs(hlm_hio_ignore_val[] - _ctrl_unset_double) < 1e-10 &&
        fates_endrun("FATES dimension/parameter unset: hio_ignore")
    hlm_ipedof[] == _ctrl_unset_int &&
        fates_endrun("index for the HLMs pedotransfer function unset: hlm_ipedof")
    hlm_parteh_mode[] == _ctrl_unset_int &&
        fates_endrun("switch deciding the plant reactive transport model is unset: hlm_parteh_mode")
    hlm_seeddisp_cadence[] == _ctrl_unset_int &&
        fates_endrun("switch defining seed dispersal cadence is unset: hlm_seeddisp_cadence")
    hlm_hist_level_dynam[] == _ctrl_unset_int &&
        fates_endrun("switch defining dynamics history level is unset: hlm_hist_level_dynam")
    hlm_hist_level_hifrq[] == _ctrl_unset_int &&
        fates_endrun("switch defining high-frequency history level is unset: hlm_hist_level_hifrq")
    hlm_use_ch4[] == _ctrl_unset_int &&
        fates_endrun("switch for the HLMs CH4 module unset: hlm_use_ch4")
    hlm_use_vertsoilc[] == _ctrl_unset_int &&
        fates_endrun("switch for the HLMs soil carbon discretization unset: hlm_use_vertsoilc")
    hlm_spitfire_mode[] == _ctrl_unset_int &&
        fates_endrun("switch for SPITFIRE unset: hlm_spitfire_mode")
    hlm_sf_nofire_def[] == _ctrl_unset_int &&
        fates_endrun("definition of no-fire mode unset: hlm_sf_nofire_def")
    hlm_sf_scalar_lightning_def[] == _ctrl_unset_int &&
        fates_endrun("definition of scalar lightning mode unset: hlm_sf_scalar_lightning_def")
    hlm_sf_successful_ignitions_def[] == _ctrl_unset_int &&
        fates_endrun("definition of successful ignition mode unset: hlm_sf_successful_ignitions_def")
    hlm_sf_anthro_ignitions_def[] == _ctrl_unset_int &&
        fates_endrun("definition of anthro-ignition mode unset: hlm_sf_anthro_ignitions_def")

    if strip(hlm_name[]) == "CLM" && hlm_parteh_mode[] == 2
        if sum(abs, pft.prescribed_puptake) < nearzero &&
           sum(abs, pft.prescribed_nuptake) < nearzero
            fates_endrun("PARTEH hypothesis 2 only viable with forced uptake BCs for CLM")
        end
    end

    hlm_use_fixed_biogeog[] == _ctrl_unset_int &&
        fates_endrun("switch for fixed biogeog unset: hlm_use_fixed_biogeog")
    hlm_use_nocomp[] == _ctrl_unset_int &&
        fates_endrun("switch for no competition mode unset: hlm_use_nocomp")
    hlm_use_sp[] == _ctrl_unset_int &&
        fates_endrun("switch for SP mode unset: hlm_use_sp")
    hlm_use_cohort_age_tracking[] == _ctrl_unset_int &&
        fates_endrun("switch for cohort_age_tracking unset: hlm_use_cohort_age_tracking")
    (hlm_use_sp[] == itrue && hlm_use_nocomp[] == ifalse) &&
        fates_endrun("SP cannot be on if nocomp mode is off")
    (hlm_use_sp[] == itrue && hlm_use_fixed_biogeog[] == ifalse) &&
        fates_endrun("SP cannot be on if fixed biogeog mode is off")

    if fates_global_verbose()
        println(fates_log(), "Checked. All control parameters sent to FATES.")
    end
    return nothing
end

# ===========================================================================
# FatesReportParameters — drive the various parameter report / check / derive
# / transfer routines.
# ===========================================================================

"""
    FatesReportParameters(masterproc)

Report + check + derive + transfer the FATES parameters. Mirrors Fortran
`FatesReportParameters`.
"""
function FatesReportParameters(masterproc::Bool)
    FatesReportPFTParams(masterproc)
    FatesReportParams(masterproc)
    PRTDerivedParams!()             # Update PARTEH derived constants
    FatesCheckParams(masterproc)    # Check general fates parameters
    PRTCheckParams(masterproc)      # Check PARTEH parameters
    SpitFireCheckParams(masterproc)
    TransferRadParams()
    return nothing
end

# ===========================================================================
# UpdateFatesRMeansTStep — per-step running-mean update of the FATES buffers.
# Walks each site's patch/cohort lists, updates the tveg* + seedling running
# means, and accumulates the site NPP EMA.
# ===========================================================================

# Solar radiation in the shortwave band used for seedling PAR (i.e. PAR).
const _rmeans_ipar = 1
# ema_npp time-scale: 10 day (10*48 steps).
const _ema_npp_tscale = 480.0

"""
    UpdateFatesRMeansTStep(sites, bc_in, bc_out)

Update the FATES running-mean buffers on the model timestep. For each site:
walk the (non-bareground) patches, update `tveg24`/`tveg_lpa`/`tveg_longterm`
from `bc_in.t_veg_pa`, update the seedling-layer PAR/SMP/MDD running means (TRS
regeneration only), accumulate the smoothed site NPP, and copy it to
`bc_out.ema_npp`. Mirrors Fortran `UpdateFatesRMeansTStep`.
"""
function UpdateFatesRMeansTStep(sites::AbstractVector, bc_in::AbstractVector, bc_out::AbstractVector)
    pft = edpftvarcon_inst()
    p   = ed_params()

    for s in 1:length(sites)
        ifp = 0
        site_npp = 0.0
        cpatch = sites[s].oldest_patch
        while cpatch !== nothing
            if cpatch.patchno != 0
                ifp += 1
                UpdateRMean!(cpatch.tveg24,       bc_in[s].t_veg_pa[ifp])
                UpdateRMean!(cpatch.tveg_lpa,      bc_in[s].t_veg_pa[ifp])
                UpdateRMean!(cpatch.tveg_longterm, bc_in[s].t_veg_pa[ifp])

                # Update the seedling-layer running means (TRS regeneration only).
                if p.regeneration_model == TRS_regeneration
                    seedling_par_high, par_high_frac, seedling_par_low, par_low_frac =
                        SeedlingParPatch(cpatch,
                            bc_in[s].solad_parb[ifp, _rmeans_ipar] +
                            bc_in[s].solai_parb[ifp, _rmeans_ipar])

                    new_seedling_layer_par = seedling_par_high * par_high_frac +
                                             seedling_par_low * par_low_frac

                    UpdateRMean!(cpatch.seedling_layer_par24, new_seedling_layer_par)
                    UpdateRMean!(cpatch.sdlng_mort_par,       new_seedling_layer_par)
                    UpdateRMean!(cpatch.sdlng2sap_par,        new_seedling_layer_par)

                    for ipft in 1:numpft[]
                        # Soil moisture at the seedling rooting depth for this pft.
                        ilayer_seedling_root = argmin(
                            abs.(bc_in[s].z_sisl .- pft.seedling_root_depth[ipft]))
                        new_seedling_layer_smp = bc_in[s].smp_sl[ilayer_seedling_root]

                        # New moisture-deficit-day value for this pft.
                        new_seedling_mdd = (abs(pft.seedling_psi_crit[ipft]) -
                            abs(new_seedling_layer_smp)) * (-1.0) * p.sdlng_mdd_timescale
                        # If negative, the soil is wetter than smp_crit -> deficit 0.
                        if new_seedling_mdd < 0.0
                            new_seedling_mdd = 0.0
                        end

                        UpdateRMean!(cpatch.sdlng_emerg_smp[ipft].p, new_seedling_layer_smp)
                        UpdateRMean!(cpatch.sdlng_mdd[ipft].p,       new_seedling_mdd)
                    end
                end

                # Accumulate site NPP from each non-new cohort [kgC/plant/yr] -> [gC/m2/s].
                ccohort = cpatch.tallest
                while ccohort !== nothing
                    if !ccohort.isnew
                        site_npp += ccohort.npp_acc_hold * ccohort.n * area_inv *
                                    g_per_kg * hlm_days_per_year[] / sec_per_day
                    end
                    ccohort = ccohort.shorter
                end
            end
            cpatch = cpatch.younger
        end

        # Smoothed site NPP [gc/m2/yr]. (Fortran reads sites(s) once, AFTER the
        # patch loop — i.e. for the LAST site index. Preserved here verbatim:
        # the smoothing/EMA-write below applies to site `s`, the value of the
        # loop variable on each iteration, which for the per-iteration body is
        # this site. The Fortran's trailing block lives OUTSIDE the patch loop
        # but INSIDE the site loop, so it runs once per site — matched here.)
        if sites[s].ema_npp < -9000.0
            sites[s].ema_npp = site_npp
        else
            sites[s].ema_npp = (1.0 - 1.0 / _ema_npp_tscale) * sites[s].ema_npp +
                               (1.0 / _ema_npp_tscale) * site_npp
        end
        bc_out[s].ema_npp = sites[s].ema_npp
    end

    return nothing
end

# ===========================================================================
# SeedlingParPatch — area-weighted PAR intensity for seedlings in a patch.
# Returns (seedling_par_high, par_high_frac, seedling_par_low, par_low_frac).
# ===========================================================================

"""
    SeedlingParPatch(cpatch, atm_par)
        -> (seedling_par_high, par_high_frac, seedling_par_low, par_low_frac)

Compute the high/low intensity PAR for seedlings in `cpatch` from the canopy-top
PAR `atm_par`, area-weighting light penetration through the lowest canopy
layer(s). Mirrors Fortran `SeedlingParPatch` (returns the four out-args as a
tuple instead of mutating).
"""
function SeedlingParPatch(cpatch, atm_par::Real)
    # Start with the assumption of a single canopy layer.
    seedling_par_high = float(atm_par)
    par_high_frac     = 1.0 - cpatch.total_canopy_area
    par_low_frac      = cpatch.total_canopy_area
    seedling_par_low  = 0.0

    # Work up through the canopy layers from the bottom layer.
    for cl in cpatch.ncl_p:-1:max(1, cpatch.ncl_p - 1)
        cl_par  = 0.0
        cl_area = 0.0
        for ipft in 1:numpft[]
            iv = cpatch.nleaf[cl, ipft]
            # Skip when there are no leaf layers for this pft in this canopy layer.
            if iv != 0
                cl_par += cpatch.canopy_area_profile[cl, ipft, 1] *
                          (cpatch.parprof_pft_dir_z[cl, ipft, iv] +
                           cpatch.parprof_pft_dif_z[cl, ipft, iv])
                cl_area += cpatch.canopy_area_profile[cl, ipft, 1]
            end
        end

        # Scale the par by area, or zero it if the area is near zero.
        if cl_area > nearzero
            cl_par = cl_par / cl_area
        else
            cl_par = 0.0
        end

        if cl < cpatch.ncl_p
            # More than one layer: average light on exposed ground under the veg.
            seedling_par_high = seedling_par_low
            par_high_frac     = 1.0 - cl_area
            seedling_par_low  = cl_par
            par_low_frac      = cl_area
        else
            # Single layer: only set the seedling_par_low.
            seedling_par_low = cl_par
        end
    end

    return seedling_par_high, par_high_frac, seedling_par_low, par_low_frac
end

# ===========================================================================
# DetermineGridCellNeighbors — build the per-gridcell neighbor lists (within the
# PFT seed-dispersal max distance), recording distance + per-PFT dispersal
# probability density. Single-process port (the Fortran's MPI Allgather across
# ranks collapses to the local lat/lon arrays passed in directly).
# ===========================================================================

"""
    DetermineGridCellNeighbors(seeds, numg, gclat, gclon)
        -> Vector{neighborhood_type}

Build the seed-dispersal neighbor lists for `numg` land gridcells given their
latitudes `gclat` and longitudes `gclon`. For every pair within any PFT's
`seed_dispersal_max_dist`, links both cells as mutual neighbors (with the
great-circle distance + the per-PFT dispersal probability density), then
populates each gridcell's `neighbor_indices` array. Returns the neighborhood
vector and stores the (trivial single-process) `ncells_array`/`begg_array` on
`seeds`. Mirrors Fortran `DetermineGridCellNeighbors`.

NOTE: the Fortran gathers lat/lon across MPI ranks (procinfo / ldomain /
MPI_Allgatherv) before this loop. In this single-process port the caller passes
the already-assembled global `gclat`/`gclon` directly. The dispersal-buffer
`seeds` book-keeping (`ncells_array`/`begg_array`) is set to the trivial
single-rank values `[numg]` / `[0]`.
"""
function DetermineGridCellNeighbors(seeds::dispersal_type, numg::Integer,
                                    gclat::AbstractVector, gclon::AbstractVector)
    pft = edpftvarcon_inst()

    # Allocate + initialize the neighborhood array.
    neighbors = [neighborhood_type() for _ in 1:numg]

    # Trivial single-process decomposition book-keeping (1 rank).
    seeds.ncells_array = [Int(numg)]
    seeds.begg_array   = [0]

    # Iterate gridcell pairs; link mutual neighbors within max dispersal distance.
    for gi in 1:(numg - 1)
        for gj in (gi + 1):numg
            g2g_dist = GetNeighborDistance(gi, gj, gclat, gclon)

            if any(d -> d > g2g_dist, pft.seed_dispersal_max_dist)
                # Add gj to gi's neighbor list.
                current_neighbor = neighbor_type()
                current_neighbor.gindex  = gj
                current_neighbor.gc_dist = g2g_dist
                current_neighbor.density_prob = zeros(Float64, numpft[])
                for ipft in 1:numpft[]
                    current_neighbor.density_prob[ipft] = ProbabilityDensity(ipft, g2g_dist)
                end

                if neighbors[gi].first_neighbor !== nothing
                    neighbors[gi].last_neighbor.next_neighbor = current_neighbor
                    neighbors[gi].last_neighbor = current_neighbor
                else
                    neighbors[gi].first_neighbor = current_neighbor
                    neighbors[gi].last_neighbor  = current_neighbor
                end
                neighbors[gi].neighbor_count += 1

                # Add gi to gj's neighbor list as well.
                another_neighbor = neighbor_type()
                another_neighbor.gindex  = gi
                another_neighbor.gc_dist = current_neighbor.gc_dist
                another_neighbor.density_prob = copy(current_neighbor.density_prob)

                if neighbors[gj].first_neighbor !== nothing
                    neighbors[gj].last_neighbor.next_neighbor = another_neighbor
                    neighbors[gj].last_neighbor = another_neighbor
                else
                    neighbors[gj].first_neighbor = another_neighbor
                    neighbors[gj].last_neighbor  = another_neighbor
                end
                neighbors[gj].neighbor_count += 1
            end
        end
    end

    # Populate the gridcell-index array for each gridcell from its linked list.
    for gi in 1:numg
        current_neighbor = neighbors[gi].first_neighbor
        neighbors[gi].neighbor_indices = zeros(Int, neighbors[gi].neighbor_count)
        ni = 1
        while current_neighbor !== nothing
            neighbors[gi].neighbor_indices[ni] = current_neighbor.gindex
            ni += 1
            current_neighbor = current_neighbor.next_neighbor
        end
    end

    return neighbors
end

# ===========================================================================
# FatesReadParameters — register + read + receive the FATES parameter file via
# the HLM-provided reader.
# ===========================================================================

"""
    FatesReadParameters(param_reader)

Register the FATES / SpitFire / PRT / Synchronized parameters, read them from
the parameter file via the HLM `param_reader`, and receive (distribute) them
into the live parameter structures. Mirrors Fortran `FatesReadParameters`.
"""
function FatesReadParameters(param_reader)
    if hlm_masterproc[] == itrue
        println(fates_log(), "FatesParametersInterface :: HLM reading ED/FATES parameters")
    end

    fates_params = fates_parameters_type()
    Init!(fates_params)
    FatesRegisterParams!(fates_params)        # EDParamsMod
    SpitFireRegisterParams!(fates_params)     # SFParamsMod
    PRTRegisterParams!(fates_params)          # PRTParamsFATESMod
    RegisterParams!(FatesSynchronizedParamsInst, fates_params)

    Read!(param_reader, fates_params)

    FatesReceiveParams!(fates_params)
    SpitFireReceiveParams!(fates_params)
    PRTReceiveParams!(fates_params)
    ReceiveParams!(FatesSynchronizedParamsInst, fates_params)

    return nothing
end
