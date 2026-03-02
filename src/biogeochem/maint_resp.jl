# ==========================================================================
# Ported from: src/biogeochem/CNMRespMod.F90
# Maintenance respiration fluxes for coupled carbon-nitrogen code
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level parameters (from params_type in CNMRespMod)
# ---------------------------------------------------------------------------

"""
    MaintRespParams

Parameters for maintenance respiration, read from namelist/file.
Corresponds to `params_type` in `CNMRespMod.F90`.
"""
Base.@kwdef mutable struct MaintRespParams
    br      ::Float64 = SPVAL   # base rate for maintenance respiration (gC/gN/s)
    br_root ::Float64 = SPVAL   # base rate for maintenance respiration for roots (gC/gN/s)
end

"""
    maint_resp_read_params!(params; br=2.525e-6, br_root=SPVAL)

Read/set parameters for maintenance respiration.
Combines `CNMRespReadNML` and `readParams` from `CNMRespMod.F90`.

If `br_root` is not provided (remains SPVAL), it defaults to `br`.
"""
function maint_resp_read_params!(params::MaintRespParams;
                                  br::Float64 = 2.525e-6,
                                  br_root::Float64 = SPVAL)
    params.br = br
    if br_root == SPVAL
        params.br_root = br
    else
        params.br_root = br_root
    end
    return nothing
end

# ---------------------------------------------------------------------------
# PFT constants needed by maintenance respiration (from pftcon)
# ---------------------------------------------------------------------------

"""
    PftConMaintResp

PFT-level parameters used by the maintenance respiration routine.
Contains a subset of fields from `pftconMod` referenced in `CNMRespMod.F90`.
"""
Base.@kwdef mutable struct PftConMaintResp
    woody ::Vector{Float64} = Float64[]   # binary woody flag (1=woody, 0=not woody)
end

# ---------------------------------------------------------------------------
# cn_mresp! — Maintenance respiration
# ---------------------------------------------------------------------------

"""
    cn_mresp!(mask_soilc, mask_soilp, bounds_c, bounds_p,
              params, cn_params_share, pftcon,
              patch, canopystate, soilstate, temperature, photosyns,
              cnveg_cf, cnveg_ns;
              nlevgrnd, npcropmin, nrepr, nlevsno)

Calculate maintenance respiration fluxes (gC/m2/s) for leaf, fine root,
live stem, live coarse root, and crop reproductive tissues.

Ported from `CNMResp` in `CNMRespMod.F90`.
"""
function cn_mresp!(mask_soilc::BitVector, mask_soilp::BitVector,
                   bounds_c::UnitRange{Int}, bounds_p::UnitRange{Int},
                   params::MaintRespParams,
                   cn_params_share::CNSharedParamsData,
                   pftcon::PftConMaintResp,
                   patch::PatchData,
                   canopystate::CanopyStateData,
                   soilstate::SoilStateData,
                   temperature::TemperatureData,
                   photosyns::PhotosynthesisData,
                   cnveg_cf::CNVegCarbonFluxData,
                   cnveg_ns::CNVegNitrogenStateData;
                   nlevgrnd::Int = varpar.nlevgrnd,
                   npcropmin::Int = NPCROPMIN,
                   nrepr::Int = NREPR,
                   nlevsno::Int = varpar.nlevsno)

    # --- Aliases (matching Fortran associate block) ---
    ivt            = patch.itype
    column         = patch.column

    woody          = pftcon.woody

    frac_veg_nosno = canopystate.frac_veg_nosno_patch
    laisun         = canopystate.laisun_patch
    laisha         = canopystate.laisha_patch

    crootfr        = soilstate.crootfr_patch

    # t_soisno is column-level 2D: (ncols, nlevsno+nlevmaxurbgrnd)
    # Fortran indices: (-nlevsno+1:nlevmaxurbgrnd), so soil layer j in Julia
    # is at offset nlevsno + j
    t_soisno       = temperature.t_soisno_col
    t_ref2m        = temperature.t_ref2m_patch
    t10            = temperature.t_a10_patch

    lmrsun         = photosyns.lmrsun_patch
    lmrsha         = photosyns.lmrsha_patch
    rootstem_acc   = photosyns.rootstem_acc

    frootn         = cnveg_ns.frootn_patch
    livestemn      = cnveg_ns.livestemn_patch
    livecrootn     = cnveg_ns.livecrootn_patch
    reproductiven  = cnveg_ns.reproductiven_patch

    leaf_mr        = cnveg_cf.leaf_mr_patch
    froot_mr       = cnveg_cf.froot_mr_patch
    livestem_mr    = cnveg_cf.livestem_mr_patch
    livecroot_mr   = cnveg_cf.livecroot_mr_patch
    reproductive_mr = cnveg_cf.reproductive_mr_patch

    # --- Constants ---
    br      = params.br
    br_root = params.br_root
    Q10     = cn_params_share.Q10

    # Offset for t_soisno column indexing (Fortran -nlevsno+1:.. → Julia 1-based)
    joff = nlevsno

    # --- Column loop: temperature correction factors for each soil layer ---
    # Allocate tcsoi as a local array (ncols × nlevgrnd)
    nc = length(bounds_c)
    tcsoi = Matrix{Float64}(undef, last(bounds_c), nlevgrnd)

    for j in 1:nlevgrnd
        for c in bounds_c
            mask_soilc[c] || continue
            tcsoi[c, j] = Q10^((t_soisno[c, j + joff] - TFRZ - 20.0) / 10.0)
        end
    end

    # --- Patch loop: leaf and live wood maintenance respiration ---
    for p in bounds_p
        mask_soilp[p] || continue

        # Temperature correction from 2m air temperature
        tc = Q10^((t_ref2m[p] - TFRZ - 20.0) / 10.0)

        # Local copies of base rates (may be modified by acclimation)
        br_local      = br
        br_root_local = br_root

        # Acclimation of root and stem respiration fluxes
        if rootstem_acc
            br_local      = br_local      * 10.0^(-0.00794 * ((t10[p] - TFRZ) - 25.0))
            br_root_local = br_root_local * 10.0^(-0.00794 * ((t10[p] - TFRZ) - 25.0))
        end

        # Leaf maintenance respiration
        if frac_veg_nosno[p] == 1
            leaf_mr[p] = lmrsun[p] * laisun[p] * 12.011e-6 +
                         lmrsha[p] * laisha[p] * 12.011e-6
        else
            leaf_mr[p] = 0.0
        end

        # Live stem and live coarse root MR
        if woody[ivt[p]] == 1.0
            livestem_mr[p]  = livestemn[p]  * br_local * tc
            livecroot_mr[p] = livecrootn[p] * br_root_local * tc
        elseif ivt[p] >= npcropmin
            livestem_mr[p] = livestemn[p] * br_local * tc
            for k in 1:nrepr
                reproductive_mr[p, k] = reproductiven[p, k] * br_local * tc
            end
        end

        # --- Initialize fine root MR accumulation ---
        froot_mr[p] = 0.0
    end

    # --- Soil and patch loop: fine root maintenance respiration ---
    for j in 1:nlevgrnd
        for p in bounds_p
            mask_soilp[p] || continue
            c = column[p]

            # Local root base rate (may be modified by acclimation)
            br_root_local = br_root
            if rootstem_acc
                br_root_local = br_root_local * 10.0^(-0.00794 * ((t10[p] - TFRZ) - 25.0))
            end

            froot_mr[p] = froot_mr[p] + frootn[p] * br_root_local * tcsoi[c, j] * crootfr[p, j]
        end
    end

    return nothing
end
