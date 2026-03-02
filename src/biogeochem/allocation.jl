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
function calc_gpp_mr_availc!(mask_soilp::BitVector, bounds::UnitRange{Int},
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

    for p in bounds
        mask_soilp[p] || continue

        ivt = patch.itype[p]

        # Convert psn from umol/m2/s -> gC/m2/s
        # The input psn (psnsun and psnsha) are expressed per unit LAI
        # in the sunlit and shaded canopy. Scale by laisun and laisha.
        cnveg_cf.psnsun_to_cpool_patch[p]   = photosyns.psnsun_patch[p] * canopystate.laisun_patch[p] * 12.011e-6
        cnveg_cf.psnshade_to_cpool_patch[p] = photosyns.psnsha_patch[p] * canopystate.laisha_patch[p] * 12.011e-6

        if use_c13 && c13_cnveg_cf !== nothing
            c13_cnveg_cf.psnsun_to_cpool_patch[p]   = photosyns.c13_psnsun_patch[p] * canopystate.laisun_patch[p] * 12.011e-6
            c13_cnveg_cf.psnshade_to_cpool_patch[p] = photosyns.c13_psnsha_patch[p] * canopystate.laisha_patch[p] * 12.011e-6
        end

        if use_c14 && c14_cnveg_cf !== nothing
            c14_cnveg_cf.psnsun_to_cpool_patch[p]   = photosyns.c14_psnsun_patch[p] * canopystate.laisun_patch[p] * 12.011e-6
            c14_cnveg_cf.psnshade_to_cpool_patch[p] = photosyns.c14_psnsha_patch[p] * canopystate.laisha_patch[p] * 12.011e-6
        end

        gpp = cnveg_cf.psnsun_to_cpool_patch[p] + cnveg_cf.psnshade_to_cpool_patch[p]
        cnveg_cf.gpp_before_downreg_patch[p] = gpp

        # get the time step total maintenance respiration
        mr = cnveg_cf.leaf_mr_patch[p] + cnveg_cf.froot_mr_patch[p]
        if pftcon.woody[ivt] == 1.0
            mr = mr + cnveg_cf.livestem_mr_patch[p] + cnveg_cf.livecroot_mr_patch[p]
        elseif ivt >= npcropmin
            if crop.croplive_patch[p]
                reproductive_mr_tot = 0.0
                for k in 1:nrepr
                    reproductive_mr_tot = reproductive_mr_tot + cnveg_cf.reproductive_mr_patch[p, k]
                end
                mr = mr + cnveg_cf.livestem_mr_patch[p] + reproductive_mr_tot
            end
        end

        # For Matrix solution if mr is very small set it to zero
        if mr < -1.0e-15 && use_matrixcn
            mr = 0.0
        end

        # carbon flux available for allocation
        availc = gpp - mr
        cnveg_cf.availc_patch[p] = availc

        # If mr > gpp, then some mr comes from gpp, the rest from cpool (xsmr)
        if mr > 0.0 && availc < 0.0
            curmr = gpp
            curmr_ratio = curmr / mr
        else
            curmr_ratio = 1.0
        end

        cnveg_cf.leaf_curmr_patch[p]      = cnveg_cf.leaf_mr_patch[p] * curmr_ratio
        cnveg_cf.leaf_xsmr_patch[p]       = cnveg_cf.leaf_mr_patch[p] - cnveg_cf.leaf_curmr_patch[p]
        cnveg_cf.froot_curmr_patch[p]     = cnveg_cf.froot_mr_patch[p] * curmr_ratio
        cnveg_cf.froot_xsmr_patch[p]      = cnveg_cf.froot_mr_patch[p] - cnveg_cf.froot_curmr_patch[p]
        cnveg_cf.livestem_curmr_patch[p]  = cnveg_cf.livestem_mr_patch[p] * curmr_ratio
        cnveg_cf.livestem_xsmr_patch[p]   = cnveg_cf.livestem_mr_patch[p] - cnveg_cf.livestem_curmr_patch[p]
        cnveg_cf.livecroot_curmr_patch[p] = cnveg_cf.livecroot_mr_patch[p] * curmr_ratio
        cnveg_cf.livecroot_xsmr_patch[p]  = cnveg_cf.livecroot_mr_patch[p] - cnveg_cf.livecroot_curmr_patch[p]
        for k in 1:nrepr
            cnveg_cf.reproductive_curmr_patch[p, k] = cnveg_cf.reproductive_mr_patch[p, k] * curmr_ratio
            cnveg_cf.reproductive_xsmr_patch[p, k]  = cnveg_cf.reproductive_mr_patch[p, k] - cnveg_cf.reproductive_curmr_patch[p, k]
        end

        # no allocation when available c is negative
        cnveg_cf.availc_patch[p] = max(cnveg_cf.availc_patch[p], 0.0)

        # test for an xsmrpool deficit
        if cnveg_cs.xsmrpool_patch[p] < 0.0
            xsmrpool_recover = -cnveg_cs.xsmrpool_patch[p] / (dayscrecover * SECSPDAY)
            if xsmrpool_recover < cnveg_cf.availc_patch[p]
                cnveg_cf.availc_patch[p] = cnveg_cf.availc_patch[p] - xsmrpool_recover
            else
                xsmrpool_recover = cnveg_cf.availc_patch[p]
                cnveg_cf.availc_patch[p] = 0.0
            end
            cnveg_cf.xsmrpool_recover_patch[p] = xsmrpool_recover
            cnveg_cf.cpool_to_xsmrpool_patch[p] = xsmrpool_recover
        end

    end

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
    crop_phase_out = Vector{Float64}(undef, length(mask_pcropp))
    crop_phase!(mask_pcropp, crop, cnveg_state, crop_phase_out)

    for p in bounds
        mask_pcropp[p] || continue

        ivt = patch.itype[p]

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
function calc_allometry!(mask_soilp::BitVector, bounds::UnitRange{Int},
                          pftcon::PftConAllocation,
                          cn_shared_params::CNSharedParamsData,
                          patch::PatchData,
                          cnveg_cf::CNVegCarbonFluxData,
                          cnveg_state::CNVegStateData;
                          npcropmin::Int=17,
                          nrepr::Int=NREPR)

    use_fun = cn_shared_params.use_fun

    for p in bounds
        mask_soilp[p] || continue

        ivt = patch.itype[p]

        f1 = pftcon.froot_leaf[ivt]
        f2 = pftcon.croot_stem[ivt]

        # modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0
        if pftcon.stem_leaf[ivt] == -1.0
            f3 = (2.7 / (1.0 + exp(-0.004 * (cnveg_cf.annsum_npp_patch[p] - 300.0)))) - 0.4
        else
            f3 = pftcon.stem_leaf[ivt]
        end

        f4  = pftcon.flivewd[ivt]
        if ivt >= npcropmin
            g1 = 0.25
        else
            g1 = pftcon.grperc[ivt]
        end
        cnl  = pftcon.leafcn[ivt]
        cnfr = pftcon.frootcn[ivt]
        cnlw = pftcon.livewdcn[ivt]
        cndw = pftcon.deadwdcn[ivt]

        # based on available C, use constant allometric relationships to
        # determine N requirements
        if !use_fun
            g1a = g1
        else
            g1a = 0.0
        end

        if pftcon.woody[ivt] == 1.0
            cnveg_state.c_allometry_patch[p] = (1.0 + g1a) * (1.0 + f1 + f3 * (1.0 + f2))
            cnveg_state.n_allometry_patch[p] = 1.0 / cnl + f1 / cnfr +
                (f3 * f4 * (1.0 + f2)) / cnlw +
                (f3 * (1.0 - f4) * (1.0 + f2)) / cndw
        elseif ivt >= npcropmin  # skip generic crops
            cng = pftcon.graincn[ivt]
            f1 = cnveg_state.aroot_patch[p] / cnveg_state.aleaf_patch[p]
            f3 = cnveg_state.astem_patch[p] / cnveg_state.aleaf_patch[p]
            f5 = Vector{Float64}(undef, nrepr)
            for k in 1:nrepr
                f5[k] = cnveg_state.arepr_patch[p, k] / cnveg_state.aleaf_patch[p]
            end
            f5_tot   = 0.0
            f5_n_tot = 0.0
            for k in 1:nrepr
                f5_tot   = f5_tot + f5[k]
                f5_n_tot = f5_n_tot + f5[k] / cng
            end
            cnveg_state.c_allometry_patch[p] = (1.0 + g1a) * (1.0 + f1 + f5_tot + f3 * (1.0 + f2))
            cnveg_state.n_allometry_patch[p] = 1.0 / cnl + f1 / cnfr + f5_n_tot +
                (f3 * f4 * (1.0 + f2)) / cnlw +
                (f3 * (1.0 - f4) * (1.0 + f2)) / cndw
        else
            cnveg_state.c_allometry_patch[p] = 1.0 + g1a + f1 + f1 * g1a
            cnveg_state.n_allometry_patch[p] = 1.0 / cnl + f1 / cnfr
        end

    end

    return nothing
end
