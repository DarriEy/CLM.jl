# FatesBstressMod.jl
# Julia port of FATES src/fates/biogeophys/FatesBstressMod.F90
#
# Computes the stress impact on transpiration from salinity (and, in the
# upstream comments, sulphide) in soils. The water-stress factor (btran) itself
# is computed separately in EDBtranMod (or in plant hydraulics); this module
# adds the *salinity* stress factor `bstress_sal_ft`, a root-density-weighted
# sigmoidal limitation that closely parallels `btran_ed!` in EDBtranMod.jl.
#
# Single public procedure:
#   * btran_sal_stress_fates!  - per patch, per PFT, accumulate a
#                             salinity-limited, root-density-weighted resistance
#                             over the soil column into cpatch.bstress_sal_ft(ft),
#                             only where liquid water is available
#                             (check_layer_water, reused from EDBtranMod).
#
# Translation notes (per project conventions):
#   * Multi-site entry point takes vectors (sites/bc_in :: AbstractVector),
#     matching EDBtranMod.btran_ed! and the other FATES coupling modules. The
#     single-site physics is a `_btran_sal_stress_site!` helper that takes the
#     resolved csite/bc_in_s, mirroring the in-module patch loop.
#   * The Fortran patch loop is oldest->youngest over the linked list; it does
#     NOT skip bareground (verbatim — the Fortran routine has no bareground
#     guard, unlike btran_ed!). Carried as-is.
#   * `numpft`/`maxpft` are the module Ref globals / consts already merged.
#   * `sites(s)%rootfrac_scr` filled by `set_root_fraction` (FatesAllometryMod).
#   * `check_layer_water` is reused directly from EDBtranMod.jl (already merged).
#   * The salinity sigmoid coefficients (1.244, 0.186, -0.132) are FATES magic
#     numbers preserved verbatim from the Fortran.
#
# Quirk preserved: the Fortran header comments mention `bstress_sul_ft` /
#   `sulphide_sl` (sulphide stress), but the routine never actually computes
#   them — only the salinity factor is calculated. The Julia types likewise
#   carry only `bstress_sal_ft` / `salinity_sl`, so no sulphide path is ported.
#
# Deps (all already merged into the CLM module): FatesConstantsMod
#   (t_water_freeze_k_1atm), EDParamsMod (maxpft), FatesInterfaceTypesMod
#   (numpft, bc_in_type, salinity_sl), EDTypesMod (ed_site_type), FatesPatchMod
#   (fates_patch_type), FatesAllometryMod (set_root_fraction), EDBtranMod
#   (check_layer_water).
# Standalone — NOT added to CLMInstances or any dual-copied struct.

# =====================================================================================

"""
    btran_sal_stress_fates!(sites, bc_in)

Calculate the salinity transpiration-stress factor `bstress_sal_ft(ft)` for
every site. Mirrors the Fortran `btran_sal_stress_fates`.

Boundary conditions in : `bc_in[s].salinity_sl` (soil salinity [ppt]),
`bc_in[s].h2o_liqvol_sl`, `bc_in[s].tempk_sl`,
`bc_in[s].max_rooting_depth_index_col`, `bc_in[s].nlevsoil`.
Output: `cpatch.bstress_sal_ft(ft)` on each patch of each site.
"""
function btran_sal_stress_fates!(sites::AbstractVector, bc_in::AbstractVector)
    for s in eachindex(sites)
        _btran_sal_stress_site!(sites[s], bc_in[s])
    end
    return nothing
end

"""
    _btran_sal_stress_site!(csite, bc_in)

Single-site core of [`btran_sal_stress_fates!`](@ref). Traverses the patch list
oldest->youngest. For each patch and each PFT:

  1. Reset `cpatch.bstress_sal_ft(ft)` to 0.
  2. Set the root-fraction profile (`set_root_fraction`).
  3. Over the soil column, where liquid water is available
     ([`check_layer_water`](@ref)), accumulate the root-density-weighted
     salinity resistance `rootfrac_scr(j) * rresis` into `bstress_sal_ft(ft)`,
     where `rresis = min(1.244 / (1 + exp((0.186 - salinity)/(-0.132))), 1)`.
"""
function _btran_sal_stress_site!(csite::ed_site_type, bc_in::bc_in_type)
    np = numpft[]
    nlevsoil = bc_in.nlevsoil

    cpatch = csite.oldest_patch
    while cpatch !== nothing
        for ft in 1:np
            cpatch.bstress_sal_ft[ft] = 0.0

            set_root_fraction(csite.rootfrac_scr, ft, csite.zi_soil;
                              max_nlevroot=bc_in.max_rooting_depth_index_col)

            for j in 1:nlevsoil
                # Calculations are only relevant where liquid water exists
                # (see clm_fates%wrap_btran for the CLM/ELM calculation).
                if check_layer_water(bc_in.h2o_liqvol_sl[j], bc_in.tempk_sl[j])
                    salinity_node = bc_in.salinity_sl[j]

                    rresis = min(1.244 / (1 + exp((0.186 - salinity_node) / (-0.132))), 1.0)

                    cpatch.bstress_sal_ft[ft] =
                        cpatch.bstress_sal_ft[ft] + csite.rootfrac_scr[j] * rresis
                end
            end  # j
        end  # PFT

        cpatch = cpatch.younger
    end

    return nothing
end

# =====================================================================================
