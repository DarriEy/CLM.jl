# EDBtranMod.jl
# Julia port of FATES src/fates/biogeophys/EDBtranMod.F90
#
# Computes the FATES soil-moisture-limited transpiration wetness factor (BTRAN)
# and the root-soil water uptake distribution (ROOTR). Three procedures:
#   * check_layer_water     - logical: is a soil layer's liquid water available
#                             for uptake (positive liquid volume AND warmer than
#                             the soil freezing threshold)?
#   * get_active_suction_layers! - per-site, fill bc_out.active_suction_sl with
#                             check_layer_water over the soil column (or all
#                             .false. if the site's btran filter is off).
#   * btran_ed!             - the main routine: per non-bareground patch, per PFT,
#                             accumulate a root-density-weighted suction-limited
#                             resistance over the soil column into btran_ft, then
#                             normalize to a per-layer uptake fraction, weight the
#                             layers/PFTs by the cohort leaf-area-weighted stomatal
#                             conductance (g_sb_laweight) into bc_out.rootr_pasl,
#                             and (when plant hydraulics is off) the diagnostic
#                             bc_out.btran_pa; finally re-normalize rootr to sum 1.
#                             When plant hydraulics is on, btran_pa instead comes
#                             from BTranForHLMDiagnosticsFromCohortHydr.
#
# Translation notes (per project conventions):
#   * Multi-site entry points take vectors (sites/bc_in/bc_out :: AbstractVector),
#     matching the existing FATES coupling modules (FatesSoilBGCFluxMod /
#     FatesPlantHydraulicsMod). The single-site physics is a `_site!` helper that
#     takes the resolved csite/bc_in_s/bc_out_s, mirroring the in-module loop.
#   * `logical function check_layer_water` -> Bool-returning `check_layer_water`.
#   * Fortran filter loops over the oldest->youngest patch linked list -> while
#     loops over `csite.oldest_patch ... cpatch.younger`, skipping bareground.
#   * `EDPftvarcon_inst%smpsc/smpso` accessed via `EDPftvarcon_inst[]`.
#   * `sites(s)%rootfrac_scr` filled by `set_root_fraction` (already ported in
#     FatesAllometryMod). Fortran `numpft`/`maxpft`/`hlm_use_planthydro` are the
#     module Ref globals / consts already merged.
#   * The Fortran `debug` write to the FATES log on the rootr renormalization is
#     preserved behind a module-local `const _edbtran_debug = false`.
#
# Deps (all already merged into the CLM module): FatesConstantsMod
#   (t_water_freeze_k_1atm, itrue, ifalse, nearzero, nocomp_bareground),
#   EDParamsMod (maxpft, soil_tfrz_thresh), FatesInterfaceTypesMod
#   (numpft, hlm_use_planthydro, bc_in_type, bc_out_type), EDTypesMod
#   (ed_site_type), FatesPatchMod (fates_patch_type), FatesCohortMod
#   (fates_cohort_type), EDPftvarcon (EDPftvarcon_inst), FatesAllometryMod
#   (set_root_fraction), FatesPlantHydraulicsMod
#   (BTranForHLMDiagnosticsFromCohortHydr), FatesGlobals (fates_log).
# Standalone — NOT added to CLMInstances or any dual-copied struct.

# Mirrors the Fortran module-local `logical, parameter :: debug = .false.`
const _edbtran_debug = false

# =====================================================================================

"""
    check_layer_water(h2o_liq_vol, tempk) -> Bool

Return `true` when a soil layer's water is available for plant uptake: it must
hold positive liquid water (`h2o_liq_vol > 0`) AND be warmer than the soil
freezing threshold (`tempk > soil_tfrz_thresh + t_water_freeze_k_1atm`).
Mirrors the Fortran `check_layer_water`.
"""
function check_layer_water(h2o_liq_vol::Real, tempk::Real)
    check_layer_water = false
    if h2o_liq_vol > 0.0
        if tempk > soil_tfrz_thresh + t_water_freeze_k_1atm
            check_layer_water = true
        end
    end
    return check_layer_water
end

# =====================================================================================

"""
    get_active_suction_layers!(sites, bc_in, bc_out)

For each site, fill `bc_out[s].active_suction_sl[j]` with whether soil layer `j`
can have plant water uptake ([`check_layer_water`](@ref) of its liquid volume and
temperature). If the site's `bc_in[s].filter_btran` is off, the whole column is
set inactive (`.false.`). Mirrors the Fortran `get_active_suction_layers`.
"""
function get_active_suction_layers!(sites::AbstractVector, bc_in::AbstractVector,
                                    bc_out::AbstractVector)
    for s in eachindex(sites)
        if bc_in[s].filter_btran
            for j in 1:bc_in[s].nlevsoil
                bc_out[s].active_suction_sl[j] =
                    check_layer_water(bc_in[s].h2o_liqvol_sl[j], bc_in[s].tempk_sl[j])
            end
        else
            fill!(bc_out[s].active_suction_sl, false)
        end
    end
    return nothing
end

# =====================================================================================

"""
    btran_ed!(sites, bc_in, bc_out)

Calculate the transpiration wetness function (BTRAN) and the root water-uptake
distribution (ROOTR) for every site. See [`_btran_ed_site!`](@ref) for the
per-site physics. When plant hydraulics is enabled (`hlm_use_planthydro == itrue`)
the patch-level diagnostic `bc_out.btran_pa` is instead computed by
[`BTranForHLMDiagnosticsFromCohortHydr`](@ref). Mirrors the Fortran `btran_ed`.

Boundary conditions in : `bc_in[s].eff_porosity_sl`, `watsat_sl`, `smp_sl`,
`h2o_liqvol_sl`, `tempk_sl`, `max_rooting_depth_index_col`, `nlevsoil`.
Boundary conditions out: `bc_out[s].rootr_pasl` (root uptake distribution),
`bc_out[s].btran_pa` (patch wetness factor).
"""
function btran_ed!(sites::AbstractVector, bc_in::AbstractVector,
                   bc_out::AbstractVector)
    for s in eachindex(sites)
        _btran_ed_site!(sites[s], bc_in[s], bc_out[s])
    end

    if hlm_use_planthydro[] == itrue
        BTranForHLMDiagnosticsFromCohortHydr(sites, bc_out)
    end

    return nothing
end

"""
    _btran_ed_site!(csite, bc_in, bc_out)

Single-site core of [`btran_ed!`](@ref). Traverses the patch list oldest->youngest
(skipping bareground), and for each veg patch:

  1. For each PFT, set the root-fraction profile (`set_root_fraction`), then over
     the soil column accumulate the suction-limited, root-density-weighted
     resistance `root_resis(ft,j)` into `cpatch.btran_ft(ft)` (only where liquid
     water is available). Layers are then normalized by `btran_ft`.
  2. Build the PFT-weighted stomatal conductance `pftgs` by summing each cohort's
     leaf-area-weighted conductance (`g_sb_laweight`) into its PFT slot.
  3. Distribute the layer resistances into `bc_out.rootr_pasl(ifp,j)`, weighting by
     `pftgs(ft)/sum_pftgs` (or uniformly `1/numpft` on the first timestep when
     `sum_pftgs == 0`).
  4. If plant hydraulics is off, form the patch diagnostic `bc_out.btran_pa(ifp)`
     the same weighted way from `btran_ft`.
  5. Re-normalize `rootr_pasl(ifp,:)` so it sums to 1.
"""
function _btran_ed_site!(csite::ed_site_type, bc_in::bc_in_type, bc_out::bc_out_type)
    edpft = EDPftvarcon_inst[]
    smpsc = edpft.smpsc  # soil water potential at full stomatal close [mm]
    smpso = edpft.smpso  # soil water potential at full stomatal open  [mm]

    np = numpft[]
    nlevsoil = bc_in.nlevsoil

    # Root resistance in each pft x layer (Fortran allocates per-site)
    root_resis = zeros(Float64, np, nlevsoil)

    # pft weighted stomatal conductance [m/s]
    pftgs = zeros(Float64, maxpft)

    fill!(bc_out.rootr_pasl, 0.0)

    ifp = 0
    cpatch = csite.oldest_patch
    while cpatch !== nothing
        if cpatch.nocomp_pft_label != nocomp_bareground  # only for veg patches
            ifp += 1

            # THIS SHOULD REALLY BE A COHORT LOOP ONCE WE HAVE rootfr_ft FOR
            # COHORTS (RGK -- carried verbatim from Fortran)
            fill!(root_resis, 0.0)

            for ft in 1:np
                set_root_fraction(csite.rootfrac_scr, ft, csite.zi_soil;
                                  max_nlevroot=bc_in.max_rooting_depth_index_col)

                cpatch.btran_ft[ft] = 0.0
                for j in 1:nlevsoil
                    # Calculations are only relevant where liquid water exists
                    if check_layer_water(bc_in.h2o_liqvol_sl[j], bc_in.tempk_sl[j])
                        smp_node = max(smpsc[ft], bc_in.smp_sl[j])

                        rresis = min((bc_in.eff_porosity_sl[j] / bc_in.watsat_sl[j]) *
                                     (smp_node - smpsc[ft]) / (smpso[ft] - smpsc[ft]), 1.0)

                        root_resis[ft, j] = csite.rootfrac_scr[j] * rresis

                        # root water uptake is not linearly proportional to root
                        # density, to allow proper deep root function.
                        cpatch.btran_ft[ft] = cpatch.btran_ft[ft] + root_resis[ft, j]
                    else
                        root_resis[ft, j] = 0.0
                    end
                end  # j

                # Normalize root resistances to get layer contribution to ET
                for j in 1:nlevsoil
                    if cpatch.btran_ft[ft] > nearzero
                        root_resis[ft, j] = root_resis[ft, j] / cpatch.btran_ft[ft]
                    else
                        root_resis[ft, j] = 0.0
                    end
                end
            end  # PFT

            # PFT-averaged point level root fraction for extraction purposes.
            # The cohort's conductance g_sb_laweight contains a weighting factor
            # based on the cohort's leaf area. units: [m/s] * [m2]
            fill!(pftgs, 0.0)
            ccohort = cpatch.tallest
            while ccohort !== nothing
                pftgs[ccohort.pft] += ccohort.g_sb_laweight
                ccohort = ccohort.shorter
            end

            # Process the boundary output; necessary for the soil-moisture sink
            # term across the soil layers in the driver/host.
            sum_pftgs = sum(@view pftgs[1:np])

            for j in 1:nlevsoil
                bc_out.rootr_pasl[ifp, j] = 0.0
                for ft in 1:np
                    if sum_pftgs > 0.0  # prevent problem on the first timestep
                        bc_out.rootr_pasl[ifp, j] = bc_out.rootr_pasl[ifp, j] +
                            root_resis[ft, j] * pftgs[ft] / sum_pftgs
                    else
                        bc_out.rootr_pasl[ifp, j] = bc_out.rootr_pasl[ifp, j] +
                            root_resis[ft, j] * 1.0 / float(np)
                    end
                end
            end

            # Calculate the BTRAN passed back to the HLM (diagnostic only). If
            # plant hydraulics is on, the patchxpft btran is overwritten later by
            # BTranForHLMDiagnosticsFromCohortHydr.
            if hlm_use_planthydro[] == ifalse
                bc_out.btran_pa[ifp] = 0.0
                for ft in 1:np
                    if sum_pftgs > 0.0
                        bc_out.btran_pa[ifp] = bc_out.btran_pa[ifp] +
                            cpatch.btran_ft[ft] * pftgs[ft] / sum_pftgs
                    else
                        bc_out.btran_pa[ifp] = bc_out.btran_pa[ifp] +
                            cpatch.btran_ft[ft] * 1.0 / np
                    end
                end
            end

            temprootr = sum(@view bc_out.rootr_pasl[ifp, 1:nlevsoil])

            if abs(1.0 - temprootr) > 1.0e-10 && temprootr > 1.0e-10
                if _edbtran_debug
                    println(fates_log(),
                            "error with rootr in canopy fluxes ", temprootr, " ", sum_pftgs)
                end
                for j in 1:nlevsoil
                    bc_out.rootr_pasl[ifp, j] = bc_out.rootr_pasl[ifp, j] / temprootr
                end
            end
        end  # not bare ground

        cpatch = cpatch.younger
    end

    return nothing
end
