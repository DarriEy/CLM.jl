# ==========================================================================
# Ported from: src/biogeochem/CNAllocationMod.F90
# CN allocation: GPP, maintenance respiration, available C, crop allocation
# fractions, and C/N allometry
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level parameter type (replaces params_type in Fortran)
# ---------------------------------------------------------------------------

"""
    AllocationParams

Parameters for the CN allocation module.

Ported from `params_type` in `CNAllocationMod.F90`.
"""
Base.@kwdef mutable struct AllocationParams
    dayscrecover::Float64 = 30.0   # number of days to recover negative cpool
end

# ---------------------------------------------------------------------------
# PFT constants needed by allocation (from pftcon)
# ---------------------------------------------------------------------------

"""
    PftConAllocation

PFT-level parameters used by allocation routines. Contains a subset of
fields from `pftconMod` that are referenced in `CNAllocationMod.F90`.
"""
Base.@kwdef mutable struct PftConAllocation
    woody       ::Vector{Float64} = Float64[]   # binary woody flag (1=woody, 0=not woody)
    froot_leaf  ::Vector{Float64} = Float64[]   # new fine root C per new leaf C (gC/gC)
    croot_stem  ::Vector{Float64} = Float64[]   # new coarse root C per new stem C (gC/gC)
    stem_leaf   ::Vector{Float64} = Float64[]   # new stem C per new leaf C (gC/gC); -1 = dynamic
    flivewd     ::Vector{Float64} = Float64[]   # fraction of new wood that is live
    leafcn      ::Vector{Float64} = Float64[]   # leaf C:N (gC/gN)
    frootcn     ::Vector{Float64} = Float64[]   # fine root C:N (gC/gN)
    livewdcn    ::Vector{Float64} = Float64[]   # live wood C:N (gC/gN)
    deadwdcn    ::Vector{Float64} = Float64[]   # dead wood C:N (gC/gN)
    graincn     ::Vector{Float64} = Float64[]   # grain C:N (gC/gN)
    grperc      ::Vector{Float64} = Float64[]   # growth respiration fraction
    # Crop-specific allocation parameters
    arooti      ::Vector{Float64} = Float64[]   # initial allocation to roots
    arootf      ::Vector{Float64} = Float64[]   # final allocation to roots
    bfact       ::Vector{Float64} = Float64[]   # exp factor for leaf allocation
    fleafi      ::Vector{Float64} = Float64[]   # initial fraction allocated to leaf
    aleaff      ::Vector{Float64} = Float64[]   # final leaf allocation coefficient
    astemf      ::Vector{Float64} = Float64[]   # final stem allocation coefficient
    allconss    ::Vector{Float64} = Float64[]   # stem allocation decline exponent
    allconsl    ::Vector{Float64} = Float64[]   # leaf allocation decline exponent
    declfact    ::Vector{Float64} = Float64[]   # decline factor for allocation shift
end

# ---------------------------------------------------------------------------
# calc_gpp_mr_availc! — Calculate GPP, MR, and available C
# ---------------------------------------------------------------------------

"""
    calc_gpp_mr_availc!(mask_soilp, bounds,
        alloc_params, pftcon, cn_shared_params,
        patch, crop, photosyns, canopystate,
        cnveg_cs, cnveg_cf;
        c13_cnveg_cf=nothing, c14_cnveg_cf=nothing,
        use_c13=false, use_c14=false,
        npcropmin=17, nrepr=NREPR)

Calculate total GPP, various maintenance respiration terms, and total
available C for allocation.

Ported from `calc_gpp_mr_availc` in `CNAllocationMod.F90`.
"""
@kernel function _alloc_gpp_mr_availc_kernel!(
        psnsun_to_cpool, psnshade_to_cpool,
        c13_psnsun_to_cpool, c13_psnshade_to_cpool,
        c14_psnsun_to_cpool, c14_psnshade_to_cpool,
        gpp_before_downreg, availc,
        leaf_curmr, leaf_xsmr, froot_curmr, froot_xsmr,
        livestem_curmr, livestem_xsmr, livecroot_curmr, livecroot_xsmr,
        reproductive_curmr, reproductive_xsmr,
        xsmrpool_recover_out, cpool_to_xsmrpool,
        @Const(mask_soilp), @Const(itype),
        @Const(psnsun_patch), @Const(psnsha_patch),
        @Const(laisun_patch), @Const(laisha_patch),
        @Const(c13_psnsun_patch), @Const(c13_psnsha_patch),
        @Const(c14_psnsun_patch), @Const(c14_psnsha_patch),
        @Const(woody),
        @Const(leaf_mr), @Const(froot_mr),
        @Const(livestem_mr), @Const(livecroot_mr),
        @Const(reproductive_mr), @Const(croplive),
        @Const(xsmrpool),
        do_c13::Bool, do_c14::Bool, use_matrixcn::Bool,
        dayscrecover::Float64, npcropmin::Int, nrepr::Int)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia

        # Convert psn from umol/m2/s -> gC/m2/s
        # The input psn (psnsun and psnsha) are expressed per unit LAI
        # in the sunlit and shaded canopy. Scale by laisun and laisha.
        psnsun_to_cpool[p]   = psnsun_patch[p] * laisun_patch[p] * 12.011e-6
        psnshade_to_cpool[p] = psnsha_patch[p] * laisha_patch[p] * 12.011e-6

        if do_c13
            c13_psnsun_to_cpool[p]   = c13_psnsun_patch[p] * laisun_patch[p] * 12.011e-6
            c13_psnshade_to_cpool[p] = c13_psnsha_patch[p] * laisha_patch[p] * 12.011e-6
        end

        if do_c14
            c14_psnsun_to_cpool[p]   = c14_psnsun_patch[p] * laisun_patch[p] * 12.011e-6
            c14_psnshade_to_cpool[p] = c14_psnsha_patch[p] * laisha_patch[p] * 12.011e-6
        end

        gpp = psnsun_to_cpool[p] + psnshade_to_cpool[p]
        gpp_before_downreg[p] = gpp

        # get the time step total maintenance respiration
        mr = leaf_mr[p] + froot_mr[p]
        if woody[ivt] == 1.0
            mr = mr + livestem_mr[p] + livecroot_mr[p]
        elseif ivt >= npcropmin
            if croplive[p]
                reproductive_mr_tot = 0.0
                for k in 1:nrepr
                    reproductive_mr_tot = reproductive_mr_tot + reproductive_mr[p, k]
                end
                mr = mr + livestem_mr[p] + reproductive_mr_tot
            end
        end

        # For Matrix solution if mr is very small set it to zero
        if mr < -1.0e-15 && use_matrixcn
            mr = 0.0
        end

        # carbon flux available for allocation
        availc_p = gpp - mr
        availc[p] = availc_p

        # If mr > gpp, then some mr comes from gpp, the rest from cpool (xsmr)
        if mr > 0.0 && availc_p < 0.0
            curmr = gpp
            curmr_ratio = curmr / mr
        else
            curmr_ratio = 1.0
        end

        leaf_curmr[p]      = leaf_mr[p] * curmr_ratio
        leaf_xsmr[p]       = leaf_mr[p] - leaf_curmr[p]
        froot_curmr[p]     = froot_mr[p] * curmr_ratio
        froot_xsmr[p]      = froot_mr[p] - froot_curmr[p]
        livestem_curmr[p]  = livestem_mr[p] * curmr_ratio
        livestem_xsmr[p]   = livestem_mr[p] - livestem_curmr[p]
        livecroot_curmr[p] = livecroot_mr[p] * curmr_ratio
        livecroot_xsmr[p]  = livecroot_mr[p] - livecroot_curmr[p]
        for k in 1:nrepr
            reproductive_curmr[p, k] = reproductive_mr[p, k] * curmr_ratio
            reproductive_xsmr[p, k]  = reproductive_mr[p, k] - reproductive_curmr[p, k]
        end

        # no allocation when available c is negative
        availc[p] = max(availc[p], 0.0)

        # test for an xsmrpool deficit
        if xsmrpool[p] < 0.0
            xsmrpool_recover = -xsmrpool[p] / (dayscrecover * SECSPDAY)
            if xsmrpool_recover < availc[p]
                availc[p] = availc[p] - xsmrpool_recover
            else
                xsmrpool_recover = availc[p]
                availc[p] = 0.0
            end
            xsmrpool_recover_out[p] = xsmrpool_recover
            cpool_to_xsmrpool[p] = xsmrpool_recover
        end
    end
end

function calc_gpp_mr_availc!(mask_soilp::AbstractVector{Bool}, bounds::UnitRange{Int},
                              alloc_params::AllocationParams,
                              pftcon::PftConAllocation,
                              cn_shared_params::CNSharedParamsData,
                              patch::PatchData,
                              crop::CropData,
                              photosyns::PhotosynthesisData,
                              canopystate::CanopyStateData,
                              cnveg_cs::CNVegCarbonStateData,
                              cnveg_cf::CNVegCarbonFluxData;
                              c13_cnveg_cf::Union{CNVegCarbonFluxData,Nothing}=nothing,
                              c14_cnveg_cf::Union{CNVegCarbonFluxData,Nothing}=nothing,
                              use_c13::Bool=false,
                              use_c14::Bool=false,
                              npcropmin::Int=17,
                              nrepr::Int=NREPR)

    dayscrecover = alloc_params.dayscrecover
    use_matrixcn = cn_shared_params.use_matrixcn

    # Resolve config branches on the host
    do_c13 = use_c13 && c13_cnveg_cf !== nothing
    do_c14 = use_c14 && c14_cnveg_cf !== nothing

    # c13/c14 output arrays: use real arrays when active, else reuse the
    # main psn arrays as inert placeholders (kernel never writes them when
    # the corresponding flag is false).
    c13_psnsun_to_cpool   = do_c13 ? c13_cnveg_cf.psnsun_to_cpool_patch   : cnveg_cf.psnsun_to_cpool_patch
    c13_psnshade_to_cpool = do_c13 ? c13_cnveg_cf.psnshade_to_cpool_patch : cnveg_cf.psnshade_to_cpool_patch
    c14_psnsun_to_cpool   = do_c14 ? c14_cnveg_cf.psnsun_to_cpool_patch   : cnveg_cf.psnsun_to_cpool_patch
    c14_psnshade_to_cpool = do_c14 ? c14_cnveg_cf.psnshade_to_cpool_patch : cnveg_cf.psnshade_to_cpool_patch

    _launch!(_alloc_gpp_mr_availc_kernel!,
        cnveg_cf.psnsun_to_cpool_patch, cnveg_cf.psnshade_to_cpool_patch,
        c13_psnsun_to_cpool, c13_psnshade_to_cpool,
        c14_psnsun_to_cpool, c14_psnshade_to_cpool,
        cnveg_cf.gpp_before_downreg_patch, cnveg_cf.availc_patch,
        cnveg_cf.leaf_curmr_patch, cnveg_cf.leaf_xsmr_patch,
        cnveg_cf.froot_curmr_patch, cnveg_cf.froot_xsmr_patch,
        cnveg_cf.livestem_curmr_patch, cnveg_cf.livestem_xsmr_patch,
        cnveg_cf.livecroot_curmr_patch, cnveg_cf.livecroot_xsmr_patch,
        cnveg_cf.reproductive_curmr_patch, cnveg_cf.reproductive_xsmr_patch,
        cnveg_cf.xsmrpool_recover_patch, cnveg_cf.cpool_to_xsmrpool_patch,
        mask_soilp, patch.itype,
        photosyns.psnsun_patch, photosyns.psnsha_patch,
        canopystate.laisun_patch, canopystate.laisha_patch,
        photosyns.c13_psnsun_patch, photosyns.c13_psnsha_patch,
        photosyns.c14_psnsun_patch, photosyns.c14_psnsha_patch,
        pftcon.woody,
        cnveg_cf.leaf_mr_patch, cnveg_cf.froot_mr_patch,
        cnveg_cf.livestem_mr_patch, cnveg_cf.livecroot_mr_patch,
        cnveg_cf.reproductive_mr_patch, crop.croplive_patch,
        cnveg_cs.xsmrpool_patch,
        do_c13, do_c14, use_matrixcn,
        dayscrecover, npcropmin, nrepr)

    return nothing
end

# ---------------------------------------------------------------------------
# calc_crop_allocation_fractions! — Crop allocation fractions
# ---------------------------------------------------------------------------

"""
    calc_crop_allocation_fractions!(mask_pcropp, bounds,
        pftcon, patch, crop, cnveg_state;
        nrepr=NREPR)

Calculate crop allocation fractions to leaf, stem, root and repr, following
AgroIBIS subroutine phenocrop.

Sets: aleaf, astem, aroot, arepr (and under some conditions aleafi, astemi).

Ported from `calc_crop_allocation_fractions` in `CNAllocationMod.F90`.
"""
function calc_crop_allocation_fractions!(mask_pcropp::BitVector, bounds::UnitRange{Int},
                                          pftcon::PftConAllocation,
                                          patch::PatchData,
                                          crop::CropData,
                                          cnveg_state::CNVegStateData;
                                          nrepr::Int=NREPR)

    # Compute crop phase for all patches in the mask
    FT = eltype(cnveg_state.c_allometry_patch)
    crop_phase_out = Vector{FT}(undef, length(mask_pcropp))
    crop_phase!(mask_pcropp, crop, cnveg_state, crop_phase_out)

    for p in bounds
        mask_pcropp[p] || continue

        ivt = patch.itype[p] + 1  # 0-based Fortran → 1-based Julia

        if crop.croplive_patch[p]
            # Phase 1 completed: leaf emergence to start of leaf decline
            if crop_phase_out[p] == cphase_leafemerge

                for k in 1:nrepr
                    cnveg_state.arepr_patch[p, k] = 0.0
                end
                if cnveg_state.peaklai_patch[p] == 1  # lai at maximum allowed
                    cnveg_state.aleaf_patch[p] = 1.0e-5
                    cnveg_state.astem_patch[p] = 0.0
                    cnveg_state.aroot_patch[p] = 1.0 - cnveg_state.aleaf_patch[p]
                else
                    cnveg_state.aroot_patch[p] = max(0.0, min(1.0, pftcon.arooti[ivt] -
                        (pftcon.arooti[ivt] - pftcon.arootf[ivt]) *
                        min(1.0, crop.hui_patch[p] / cnveg_state.gddmaturity_patch[p])))
                    fleaf = pftcon.fleafi[ivt] * (exp(-pftcon.bfact[ivt]) -
                        exp(-pftcon.bfact[ivt] * crop.hui_patch[p] / cnveg_state.huigrain_patch[p])) /
                        (exp(-pftcon.bfact[ivt]) - 1.0)
                    cnveg_state.aleaf_patch[p] = max(1.0e-5, (1.0 - cnveg_state.aroot_patch[p]) * fleaf)
                    cnveg_state.astem_patch[p] = 1.0 - cnveg_state.aleaf_patch[p] - cnveg_state.aroot_patch[p]
                end

                cnveg_state.astemi_patch[p] = cnveg_state.astem_patch[p]
                cnveg_state.aleafi_patch[p] = cnveg_state.aleaf_patch[p]

            # Phase 2 completed: grain fill
            elseif crop_phase_out[p] == cphase_grainfill
                cnveg_state.aroot_patch[p] = max(0.0, min(1.0, pftcon.arooti[ivt] -
                    (pftcon.arooti[ivt] - pftcon.arootf[ivt]) *
                    min(1.0, crop.hui_patch[p] / cnveg_state.gddmaturity_patch[p])))
                if cnveg_state.astemi_patch[p] > pftcon.astemf[ivt]
                    cnveg_state.astem_patch[p] = max(0.0, max(pftcon.astemf[ivt],
                        cnveg_state.astem_patch[p] *
                        (1.0 - min((crop.hui_patch[p] -
                        cnveg_state.huigrain_patch[p]) / ((cnveg_state.gddmaturity_patch[p] * pftcon.declfact[ivt]) -
                        cnveg_state.huigrain_patch[p]), 1.0)^pftcon.allconss[ivt])))
                end

                if cnveg_state.peaklai_patch[p] == 1
                    cnveg_state.aleaf_patch[p] = 1.0e-5
                elseif cnveg_state.aleafi_patch[p] > pftcon.aleaff[ivt]
                    cnveg_state.aleaf_patch[p] = max(1.0e-5, max(pftcon.aleaff[ivt],
                        cnveg_state.aleaf_patch[p] *
                        (1.0 - min((crop.hui_patch[p] -
                        cnveg_state.huigrain_patch[p]) / ((cnveg_state.gddmaturity_patch[p] * pftcon.declfact[ivt]) -
                        cnveg_state.huigrain_patch[p]), 1.0)^pftcon.allconsl[ivt])))
                end

                # All repr allocation goes into the last reproductive pool
                for k in 1:nrepr-1
                    cnveg_state.arepr_patch[p, k] = 0.0
                end
                cnveg_state.arepr_patch[p, nrepr] = 1.0 - cnveg_state.aroot_patch[p] -
                    cnveg_state.astem_patch[p] - cnveg_state.aleaf_patch[p]

            elseif crop_phase_out[p] == cphase_planted
                # pre emergence
                cnveg_state.aleaf_patch[p]  = 1.0
                cnveg_state.aleafi_patch[p] = 1.0
                cnveg_state.astem_patch[p]  = 0.0
                cnveg_state.astemi_patch[p] = 0.0
                cnveg_state.aroot_patch[p]  = 0.0
                for k in 1:nrepr
                    cnveg_state.arepr_patch[p, k] = 0.0
                end

            else
                error("Unexpected crop_phase: $(crop_phase_out[p]) at patch $p")
            end

        else  # .not. croplive
            cnveg_state.aleaf_patch[p]  = 1.0
            cnveg_state.aleafi_patch[p] = 1.0
            cnveg_state.astem_patch[p]  = 0.0
            cnveg_state.astemi_patch[p] = 0.0
            cnveg_state.aroot_patch[p]  = 0.0
            for k in 1:nrepr
                cnveg_state.arepr_patch[p, k] = 0.0
            end
        end

    end

    return nothing
end

# ---------------------------------------------------------------------------
# calc_allometry! — C/N allometry based on allocation fractions
# ---------------------------------------------------------------------------

"""
    calc_allometry!(mask_soilp, bounds,
        pftcon, cn_shared_params, patch,
        cnveg_cf, cnveg_state;
        npcropmin=17, nrepr=NREPR)

Calculate c_allometry and n_allometry terms based on allocation fractions.

Ported from `calc_allometry` in `CNAllocationMod.F90`.
"""
@kernel function _alloc_allometry_kernel!(
        c_allometry, n_allometry,
        @Const(mask_soilp), @Const(itype),
        @Const(froot_leaf), @Const(croot_stem), @Const(stem_leaf),
        @Const(annsum_npp), @Const(flivewd), @Const(grperc),
        @Const(leafcn), @Const(frootcn), @Const(livewdcn), @Const(deadwdcn),
        @Const(woody), @Const(graincn),
        @Const(aroot), @Const(aleaf), @Const(astem), @Const(arepr),
        use_fun::Bool, npcropmin::Int, nrepr::Int)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia

        f1 = froot_leaf[ivt]
        f2 = croot_stem[ivt]

        # modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0
        if stem_leaf[ivt] == -1.0
            f3 = (2.7 / (1.0 + exp(-0.004 * (annsum_npp[p] - 300.0)))) - 0.4
        else
            f3 = stem_leaf[ivt]
        end

        f4  = flivewd[ivt]
        if ivt >= npcropmin
            g1 = 0.25
        else
            g1 = grperc[ivt]
        end
        cnl  = leafcn[ivt]
        cnfr = frootcn[ivt]
        cnlw = livewdcn[ivt]
        cndw = deadwdcn[ivt]

        # based on available C, use constant allometric relationships to
        # determine N requirements
        if !use_fun
            g1a = g1
        else
            g1a = 0.0
        end

        if woody[ivt] == 1.0
            c_allometry[p] = (1.0 + g1a) * (1.0 + f1 + f3 * (1.0 + f2))
            n_allometry[p] = 1.0 / cnl + f1 / cnfr +
                (f3 * f4 * (1.0 + f2)) / cnlw +
                (f3 * (1.0 - f4) * (1.0 + f2)) / cndw
        elseif ivt >= npcropmin  # skip generic crops
            cng = graincn[ivt]
            f1 = aroot[p] / aleaf[p]
            f3 = astem[p] / aleaf[p]
            f5_tot   = 0.0
            f5_n_tot = 0.0
            for k in 1:nrepr
                f5k = arepr[p, k] / aleaf[p]
                f5_tot   = f5_tot + f5k
                f5_n_tot = f5_n_tot + f5k / cng
            end
            c_allometry[p] = (1.0 + g1a) * (1.0 + f1 + f5_tot + f3 * (1.0 + f2))
            n_allometry[p] = 1.0 / cnl + f1 / cnfr + f5_n_tot +
                (f3 * f4 * (1.0 + f2)) / cnlw +
                (f3 * (1.0 - f4) * (1.0 + f2)) / cndw
        else
            c_allometry[p] = 1.0 + g1a + f1 + f1 * g1a
            n_allometry[p] = 1.0 / cnl + f1 / cnfr
        end
    end
end

function calc_allometry!(mask_soilp::AbstractVector{Bool}, bounds::UnitRange{Int},
                          pftcon::PftConAllocation,
                          cn_shared_params::CNSharedParamsData,
                          patch::PatchData,
                          cnveg_cf::CNVegCarbonFluxData,
                          cnveg_state::CNVegStateData;
                          npcropmin::Int=17,
                          nrepr::Int=NREPR)

    use_fun = cn_shared_params.use_fun

    _launch!(_alloc_allometry_kernel!,
        cnveg_state.c_allometry_patch, cnveg_state.n_allometry_patch,
        mask_soilp, patch.itype,
        pftcon.froot_leaf, pftcon.croot_stem, pftcon.stem_leaf,
        cnveg_cf.annsum_npp_patch, pftcon.flivewd, pftcon.grperc,
        pftcon.leafcn, pftcon.frootcn, pftcon.livewdcn, pftcon.deadwdcn,
        pftcon.woody, pftcon.graincn,
        cnveg_state.aroot_patch, cnveg_state.aleaf_patch,
        cnveg_state.astem_patch, cnveg_state.arepr_patch,
        use_fun, npcropmin, nrepr)

    return nothing
end
