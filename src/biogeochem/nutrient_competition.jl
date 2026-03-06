# ==========================================================================
# Ported from: src/biogeochem/NutrientCompetitionCLM45defaultMod.F90
# Plant nutrient demand and competition algorithm (CLM4.5 default method)
#
# Original description:
#   Module contains different subroutines to do soil nutrient competition
#   dynamics.
#   Created by Jinyun Tang, Sep 8, 2014
#   Modified by Mariana Vertenstein, Nov 15, 2014
# ==========================================================================

# ---------------------------------------------------------------------------
# PFT constants needed by nutrient competition (from pftcon)
# ---------------------------------------------------------------------------

"""
    PftConNutrientCompetition

PFT-level parameters used by the nutrient competition routines. Contains a
subset of fields from `pftconMod` referenced in
`NutrientCompetitionCLM45defaultMod.F90`.
"""
Base.@kwdef mutable struct PftConNutrientCompetition
    woody       ::Vector{Float64} = Float64[]   # binary woody flag (1=woody, 0=not woody)
    froot_leaf  ::Vector{Float64} = Float64[]   # new fine root C per new leaf C (gC/gC)
    croot_stem  ::Vector{Float64} = Float64[]   # new coarse root C per new stem C (gC/gC)
    stem_leaf   ::Vector{Float64} = Float64[]   # new stem C per new leaf C (gC/gC); -1 = dynamic
    flivewd     ::Vector{Float64} = Float64[]   # fraction of new wood that is live
    leafcn      ::Vector{Float64} = Float64[]   # leaf C:N (gC/gN)
    frootcn     ::Vector{Float64} = Float64[]   # fine root C:N (gC/gN)
    livewdcn    ::Vector{Float64} = Float64[]   # live wood C:N (gC/gN)
    deadwdcn    ::Vector{Float64} = Float64[]   # dead wood C:N (gC/gN)
    fcur        ::Vector{Float64} = Float64[]   # fraction of allocation to current growth
    graincn     ::Vector{Float64} = Float64[]   # grain C:N (gC/gN)
    grperc      ::Vector{Float64} = Float64[]   # growth respiration fraction
    grpnow      ::Vector{Float64} = Float64[]   # growth respiration fraction released immediately
    # Crop grain-fill C:N override parameters
    fleafcn     ::Vector{Float64} = Float64[]   # leaf C:N during organ fill
    ffrootcn    ::Vector{Float64} = Float64[]   # fine root C:N during organ fill
    fstemcn     ::Vector{Float64} = Float64[]   # stem C:N during organ fill
    astemf      ::Vector{Float64} = Float64[]   # final stem allocation coefficient
    # Deciduous flags
    season_decid ::Vector{Float64} = Float64[]  # binary flag for seasonal-deciduous (0 or 1)
    stress_decid ::Vector{Float64} = Float64[]  # binary flag for stress-deciduous (0 or 1)
end

# ---------------------------------------------------------------------------
# calc_plant_nutrient_demand! -- Calculate plant nitrogen demand
# Wraps calc_plant_nitrogen_demand! (matches Fortran's calc_plant_nutrient_demand)
# ---------------------------------------------------------------------------

"""
    calc_plant_nutrient_demand!(mask_p, bounds, call_is_for_pcrop,
        pftcon, cn_shared_params, patch, crop, cnveg_state,
        cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf;
        dt, npcropmin, nrepr, ntmp_soybean, nirrig_tmp_soybean,
        ntrp_soybean, nirrig_trp_soybean)

Calculate plant nitrogen demand. This is called separately for
non-prognostic-crop points and prognostic-crop points.

Ported from `calc_plant_nutrient_demand` in
`NutrientCompetitionCLM45defaultMod.F90`.
"""
function calc_plant_nutrient_demand!(mask_p::BitVector, bounds::UnitRange{Int},
        call_is_for_pcrop::Bool,
        pftcon::PftConNutrientCompetition,
        cn_shared_params::CNSharedParamsData,
        patch::PatchData,
        crop::CropData,
        cnveg_state::CNVegStateData,
        cnveg_cs::CNVegCarbonStateData,
        cnveg_cf::CNVegCarbonFluxData,
        cnveg_ns::CNVegNitrogenStateData,
        cnveg_nf::CNVegNitrogenFluxData;
        dt::Float64,
        npcropmin::Int=NPCROPMIN,
        nrepr::Int=NREPR,
        ntmp_soybean::Int=23,
        nirrig_tmp_soybean::Int=24,
        ntrp_soybean::Int=77,
        nirrig_trp_soybean::Int=78)

    calc_plant_nitrogen_demand!(mask_p, bounds, call_is_for_pcrop,
        pftcon, cn_shared_params, patch, crop,
        cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf;
        dt=dt,
        npcropmin=npcropmin,
        nrepr=nrepr,
        ntmp_soybean=ntmp_soybean,
        nirrig_tmp_soybean=nirrig_tmp_soybean,
        ntrp_soybean=ntrp_soybean,
        nirrig_trp_soybean=nirrig_trp_soybean)

    return nothing
end

# ---------------------------------------------------------------------------
# calc_plant_nutrient_competition! -- Calculate nutrient yield rate from
# competition. Wraps calc_plant_cn_alloc!
# ---------------------------------------------------------------------------

"""
    calc_plant_nutrient_competition!(mask_soilp, bounds,
        pftcon, cn_shared_params, patch, crop, cnveg_state,
        cnveg_cs, cnveg_cf, cnveg_nf;
        fpg_col, c13_cnveg_cf, c14_cnveg_cf,
        use_c13, use_c14, npcropmin, nrepr)

Calculate nutrient yield rate from competition: distribute available N
between competing patches, allocate C and N to new growth and storage.

Ported from `calc_plant_nutrient_competition` /
`calc_plant_cn_alloc` in `NutrientCompetitionCLM45defaultMod.F90`.
"""
function calc_plant_nutrient_competition!(mask_soilp::BitVector, bounds::UnitRange{Int},
        pftcon::PftConNutrientCompetition,
        cn_shared_params::CNSharedParamsData,
        patch::PatchData,
        crop::CropData,
        cnveg_state::CNVegStateData,
        cnveg_cs::CNVegCarbonStateData,
        cnveg_cf::CNVegCarbonFluxData,
        cnveg_nf::CNVegNitrogenFluxData;
        fpg_col::Vector{Float64},
        c13_cnveg_cf::Union{CNVegCarbonFluxData,Nothing}=nothing,
        c14_cnveg_cf::Union{CNVegCarbonFluxData,Nothing}=nothing,
        use_c13::Bool=false,
        use_c14::Bool=false,
        npcropmin::Int=NPCROPMIN,
        nrepr::Int=NREPR)

    calc_plant_cn_alloc!(mask_soilp, bounds,
        pftcon, cn_shared_params, patch, crop,
        cnveg_state, cnveg_cs, cnveg_cf, cnveg_nf;
        fpg_col=fpg_col,
        c13_cnveg_cf=c13_cnveg_cf,
        c14_cnveg_cf=c14_cnveg_cf,
        use_c13=use_c13,
        use_c14=use_c14,
        npcropmin=npcropmin,
        nrepr=nrepr)

    return nothing
end

# ---------------------------------------------------------------------------
# calc_plant_cn_alloc! -- Private-equivalent: C/N allocation after competition
# ---------------------------------------------------------------------------

"""
    calc_plant_cn_alloc!(mask_soilp, bounds,
        pftcon, cn_shared_params, patch, crop,
        cnveg_state, cnveg_cs, cnveg_cf, cnveg_nf;
        fpg_col, c13_cnveg_cf, c14_cnveg_cf,
        use_c13, use_c14, npcropmin, nrepr)

Distribute the available N between competing patches on the basis of
relative demand, and allocate C and N to new growth and storage.

Ported from `calc_plant_cn_alloc` in
`NutrientCompetitionCLM45defaultMod.F90`.
"""
function calc_plant_cn_alloc!(mask_soilp::BitVector, bounds::UnitRange{Int},
        pftcon::PftConNutrientCompetition,
        cn_shared_params::CNSharedParamsData,
        patch::PatchData,
        crop::CropData,
        cnveg_state::CNVegStateData,
        cnveg_cs::CNVegCarbonStateData,
        cnveg_cf::CNVegCarbonFluxData,
        cnveg_nf::CNVegNitrogenFluxData;
        fpg_col::Vector{Float64},
        c13_cnveg_cf::Union{CNVegCarbonFluxData,Nothing}=nothing,
        c14_cnveg_cf::Union{CNVegCarbonFluxData,Nothing}=nothing,
        use_c13::Bool=false,
        use_c14::Bool=false,
        npcropmin::Int=NPCROPMIN,
        nrepr::Int=NREPR)

    use_fun = cn_shared_params.use_fun

    for p in bounds
        mask_soilp[p] || continue

        c = patch.column[p]
        ivt = patch.itype[p] + 1  # 0-based Fortran → 1-based Julia

        # set local allocation variables
        f1 = pftcon.froot_leaf[ivt]
        f2 = pftcon.croot_stem[ivt]

        # modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
        # constrained so that it does not go lower than 0.2 (under negative annsum_npp)
        # -1.0 means dynamic allocation (trees)
        if pftcon.stem_leaf[ivt] == -1.0
            f3 = (2.7 / (1.0 + exp(-0.004 * (cnveg_cf.annsum_npp_patch[p] - 300.0)))) - 0.4
        else
            f3 = pftcon.stem_leaf[ivt]
        end

        f4   = pftcon.flivewd[ivt]
        g1   = pftcon.grperc[ivt]
        g2   = pftcon.grpnow[ivt]
        cnl  = pftcon.leafcn[ivt]
        cnfr = pftcon.frootcn[ivt]
        cnlw = pftcon.livewdcn[ivt]
        cndw = pftcon.deadwdcn[ivt]
        fcur = pftcon.fcur[ivt]

        f5 = zeros(nrepr)

        if ivt >= npcropmin  # skip 2 generic crops
            if crop.croplive_patch[p] && !isnan(cnveg_state.aleaf_patch[p])
                f1 = cnveg_state.aroot_patch[p] / cnveg_state.aleaf_patch[p]
                f3 = cnveg_state.astem_patch[p] / cnveg_state.aleaf_patch[p]
                for k in 1:nrepr
                    f5[k] = cnveg_state.arepr_patch[p, k] / cnveg_state.aleaf_patch[p]
                end
                g1 = 0.25
            else
                f1 = 0.0
                f3 = 0.0
                for k in 1:nrepr
                    f5[k] = 0.0
                end
                g1 = 0.25
            end
        end

        if use_fun  # if we are using FUN, we get the N available from there
            cnveg_nf.sminn_to_npool_patch[p] = cnveg_nf.sminn_to_plant_fun_patch[p]
        else  # no FUN - get N available from the FPG calculation
            cnveg_nf.sminn_to_npool_patch[p] = cnveg_nf.plant_ndemand_patch[p] * fpg_col[c]
        end

        cnveg_nf.plant_nalloc_patch[p] = cnveg_nf.sminn_to_npool_patch[p] +
                                          cnveg_nf.retransn_to_npool_patch[p]

        cnveg_cf.plant_calloc_patch[p] = cnveg_nf.plant_nalloc_patch[p] *
            (cnveg_state.c_allometry_patch[p] / cnveg_state.n_allometry_patch[p])

        if !use_fun  # ORIGINAL CLM(CN) downregulation code
            cnveg_cf.excess_cflux_patch[p] = cnveg_cf.availc_patch[p] -
                                              cnveg_cf.plant_calloc_patch[p]

            # reduce gpp fluxes due to N limitation
            if cnveg_cf.gpp_before_downreg_patch[p] > 0.0
                cnveg_state.downreg_patch[p] = cnveg_cf.excess_cflux_patch[p] /
                    cnveg_cf.gpp_before_downreg_patch[p]

                cnveg_cf.psnsun_to_cpool_patch[p] = cnveg_cf.psnsun_to_cpool_patch[p] *
                    (1.0 - cnveg_state.downreg_patch[p])
                cnveg_cf.psnshade_to_cpool_patch[p] = cnveg_cf.psnshade_to_cpool_patch[p] *
                    (1.0 - cnveg_state.downreg_patch[p])

                if use_c13 && c13_cnveg_cf !== nothing
                    c13_cnveg_cf.psnsun_to_cpool_patch[p] = c13_cnveg_cf.psnsun_to_cpool_patch[p] *
                        (1.0 - cnveg_state.downreg_patch[p])
                    c13_cnveg_cf.psnshade_to_cpool_patch[p] = c13_cnveg_cf.psnshade_to_cpool_patch[p] *
                        (1.0 - cnveg_state.downreg_patch[p])
                end
                if use_c14 && c14_cnveg_cf !== nothing
                    c14_cnveg_cf.psnsun_to_cpool_patch[p] = c14_cnveg_cf.psnsun_to_cpool_patch[p] *
                        (1.0 - cnveg_state.downreg_patch[p])
                    c14_cnveg_cf.psnshade_to_cpool_patch[p] = c14_cnveg_cf.psnshade_to_cpool_patch[p] *
                        (1.0 - cnveg_state.downreg_patch[p])
                end
            end
        end  # use_fun

        # calculate new leaf C and daily fluxes to current growth and storage pools
        nlc = cnveg_cf.plant_calloc_patch[p] / cnveg_state.c_allometry_patch[p]

        cnveg_cf.cpool_to_leafc_patch[p]          = nlc * fcur
        cnveg_cf.cpool_to_leafc_storage_patch[p]  = nlc * (1.0 - fcur)
        cnveg_cf.cpool_to_frootc_patch[p]         = nlc * f1 * fcur
        cnveg_cf.cpool_to_frootc_storage_patch[p] = nlc * f1 * (1.0 - fcur)

        if pftcon.woody[ivt] == 1.0
            cnveg_cf.cpool_to_livestemc_patch[p]          = nlc * f3 * f4 * fcur
            cnveg_cf.cpool_to_livestemc_storage_patch[p]  = nlc * f3 * f4 * (1.0 - fcur)
            cnveg_cf.cpool_to_deadstemc_patch[p]          = nlc * f3 * (1.0 - f4) * fcur
            cnveg_cf.cpool_to_deadstemc_storage_patch[p]  = nlc * f3 * (1.0 - f4) * (1.0 - fcur)
            cnveg_cf.cpool_to_livecrootc_patch[p]         = nlc * f2 * f3 * f4 * fcur
            cnveg_cf.cpool_to_livecrootc_storage_patch[p] = nlc * f2 * f3 * f4 * (1.0 - fcur)
            cnveg_cf.cpool_to_deadcrootc_patch[p]         = nlc * f2 * f3 * (1.0 - f4) * fcur
            cnveg_cf.cpool_to_deadcrootc_storage_patch[p] = nlc * f2 * f3 * (1.0 - f4) * (1.0 - fcur)
        end

        if ivt >= npcropmin  # skip 2 generic crops
            cnveg_cf.cpool_to_livestemc_patch[p]          = nlc * f3 * f4 * fcur
            cnveg_cf.cpool_to_livestemc_storage_patch[p]  = nlc * f3 * f4 * (1.0 - fcur)
            cnveg_cf.cpool_to_deadstemc_patch[p]          = nlc * f3 * (1.0 - f4) * fcur
            cnveg_cf.cpool_to_deadstemc_storage_patch[p]  = nlc * f3 * (1.0 - f4) * (1.0 - fcur)
            cnveg_cf.cpool_to_livecrootc_patch[p]         = nlc * f2 * f3 * f4 * fcur
            cnveg_cf.cpool_to_livecrootc_storage_patch[p] = nlc * f2 * f3 * f4 * (1.0 - fcur)
            cnveg_cf.cpool_to_deadcrootc_patch[p]         = nlc * f2 * f3 * (1.0 - f4) * fcur
            cnveg_cf.cpool_to_deadcrootc_storage_patch[p] = nlc * f2 * f3 * (1.0 - f4) * (1.0 - fcur)
            for k in 1:nrepr
                cnveg_cf.cpool_to_reproductivec_patch[p, k]         = nlc * f5[k] * fcur
                cnveg_cf.cpool_to_reproductivec_storage_patch[p, k] = nlc * f5[k] * (1.0 - fcur)
            end
        end

        # corresponding N fluxes
        cnveg_nf.npool_to_leafn_patch[p]          = (nlc / cnl) * fcur
        cnveg_nf.npool_to_leafn_storage_patch[p]  = (nlc / cnl) * (1.0 - fcur)
        cnveg_nf.npool_to_frootn_patch[p]         = (nlc * f1 / cnfr) * fcur
        cnveg_nf.npool_to_frootn_storage_patch[p] = (nlc * f1 / cnfr) * (1.0 - fcur)

        if pftcon.woody[ivt] == 1.0
            cnveg_nf.npool_to_livestemn_patch[p]          = (nlc * f3 * f4 / cnlw) * fcur
            cnveg_nf.npool_to_livestemn_storage_patch[p]  = (nlc * f3 * f4 / cnlw) * (1.0 - fcur)
            cnveg_nf.npool_to_deadstemn_patch[p]          = (nlc * f3 * (1.0 - f4) / cndw) * fcur
            cnveg_nf.npool_to_deadstemn_storage_patch[p]  = (nlc * f3 * (1.0 - f4) / cndw) * (1.0 - fcur)
            cnveg_nf.npool_to_livecrootn_patch[p]         = (nlc * f2 * f3 * f4 / cnlw) * fcur
            cnveg_nf.npool_to_livecrootn_storage_patch[p] = (nlc * f2 * f3 * f4 / cnlw) * (1.0 - fcur)
            cnveg_nf.npool_to_deadcrootn_patch[p]         = (nlc * f2 * f3 * (1.0 - f4) / cndw) * fcur
            cnveg_nf.npool_to_deadcrootn_storage_patch[p] = (nlc * f2 * f3 * (1.0 - f4) / cndw) * (1.0 - fcur)
        end

        if ivt >= npcropmin  # skip 2 generic crops
            cng = pftcon.graincn[ivt]
            cnveg_nf.npool_to_livestemn_patch[p]          = (nlc * f3 * f4 / cnlw) * fcur
            cnveg_nf.npool_to_livestemn_storage_patch[p]  = (nlc * f3 * f4 / cnlw) * (1.0 - fcur)
            cnveg_nf.npool_to_deadstemn_patch[p]          = (nlc * f3 * (1.0 - f4) / cndw) * fcur
            cnveg_nf.npool_to_deadstemn_storage_patch[p]  = (nlc * f3 * (1.0 - f4) / cndw) * (1.0 - fcur)
            cnveg_nf.npool_to_livecrootn_patch[p]         = (nlc * f2 * f3 * f4 / cnlw) * fcur
            cnveg_nf.npool_to_livecrootn_storage_patch[p] = (nlc * f2 * f3 * f4 / cnlw) * (1.0 - fcur)
            cnveg_nf.npool_to_deadcrootn_patch[p]         = (nlc * f2 * f3 * (1.0 - f4) / cndw) * fcur
            cnveg_nf.npool_to_deadcrootn_storage_patch[p] = (nlc * f2 * f3 * (1.0 - f4) / cndw) * (1.0 - fcur)
            for k in 1:nrepr
                cnveg_nf.npool_to_reproductiven_patch[p, k]         = (nlc * f5[k] / cng) * fcur
                cnveg_nf.npool_to_reproductiven_storage_patch[p, k] = (nlc * f5[k] / cng) * (1.0 - fcur)
            end
        end

        # Growth respiration storage: carbon that needs to go into growth
        # respiration storage to satisfy all of the storage growth demands
        gresp_storage = cnveg_cf.cpool_to_leafc_storage_patch[p] +
                        cnveg_cf.cpool_to_frootc_storage_patch[p]
        if pftcon.woody[ivt] == 1.0
            gresp_storage += cnveg_cf.cpool_to_livestemc_storage_patch[p]
            gresp_storage += cnveg_cf.cpool_to_deadstemc_storage_patch[p]
            gresp_storage += cnveg_cf.cpool_to_livecrootc_storage_patch[p]
            gresp_storage += cnveg_cf.cpool_to_deadcrootc_storage_patch[p]
        end
        if ivt >= npcropmin  # skip 2 generic crops
            gresp_storage += cnveg_cf.cpool_to_livestemc_storage_patch[p]
            for k in 1:nrepr
                gresp_storage += cnveg_cf.cpool_to_reproductivec_storage_patch[p, k]
            end
        end
        cnveg_cf.cpool_to_gresp_storage_patch[p] = gresp_storage * g1 * (1.0 - g2)

    end  # end patch loop

    return nothing
end

# ---------------------------------------------------------------------------
# calc_plant_nitrogen_demand! -- Private-equivalent: N demand calculation
# ---------------------------------------------------------------------------

"""
    calc_plant_nitrogen_demand!(mask_p, bounds, call_is_for_pcrop,
        pftcon, cn_shared_params, patch, crop,
        cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf;
        dt, npcropmin, nrepr, ntmp_soybean, nirrig_tmp_soybean,
        ntrp_soybean, nirrig_trp_soybean)

Calculate plant nitrogen demand, retranslocated N deployment, and
crop grain-fill retranslocation fluxes.

Sets the following output variables used elsewhere:
- `plant_ndemand_patch`
- `retransn_to_npool_patch`
- `leafn_to_retransn_patch`
- `frootn_to_retransn_patch`
- `livestemn_to_retransn_patch`

Ported from `calc_plant_nitrogen_demand` in
`NutrientCompetitionCLM45defaultMod.F90`.
"""
function calc_plant_nitrogen_demand!(mask_p::BitVector, bounds::UnitRange{Int},
        call_is_for_pcrop::Bool,
        pftcon::PftConNutrientCompetition,
        cn_shared_params::CNSharedParamsData,
        patch::PatchData,
        crop::CropData,
        cnveg_state::CNVegStateData,
        cnveg_cs::CNVegCarbonStateData,
        cnveg_cf::CNVegCarbonFluxData,
        cnveg_ns::CNVegNitrogenStateData,
        cnveg_nf::CNVegNitrogenFluxData;
        dt::Float64,
        npcropmin::Int=NPCROPMIN,
        nrepr::Int=NREPR,
        ntmp_soybean::Int=23,
        nirrig_tmp_soybean::Int=24,
        ntrp_soybean::Int=77,
        nirrig_trp_soybean::Int=78)

    use_fun = cn_shared_params.use_fun

    # loop over patches to assess the total plant N demand
    for p in bounds
        mask_p[p] || continue

        cnveg_nf.plant_ndemand_patch[p] = cnveg_cf.availc_patch[p] *
            (cnveg_state.n_allometry_patch[p] / cnveg_state.c_allometry_patch[p])

        # retranslocated N deployment depends on seasonal cycle of potential GPP
        # (requires one year run to accumulate demand)
        cnveg_state.tempsum_potential_gpp_patch[p] =
            cnveg_state.tempsum_potential_gpp_patch[p] +
            cnveg_cf.gpp_before_downreg_patch[p]

        # carry max retransn info to CN Annual Update
        cnveg_state.tempmax_retransn_patch[p] =
            max(cnveg_state.tempmax_retransn_patch[p], cnveg_ns.retransn_patch[p])
    end

    # Crop grain-fill retranslocation
    if call_is_for_pcrop
        # Compute crop phase for all patches in bounds
        crop_phase_vals = zeros(length(mask_p))
        crop_phase!(mask_p, crop, cnveg_state, crop_phase_vals)

        for p in bounds
            mask_p[p] || continue

            ivt = patch.itype[p] + 1  # 0-based Fortran → 1-based Julia

            if crop.croplive_patch[p]
                if crop_phase_vals[p] == cphase_leafemerge
                    cnveg_state.grain_flag_patch[p] = 0.0  # setting to 0 while in phase 2
                elseif crop_phase_vals[p] == cphase_grainfill
                    # Beth's retranslocation of leafn, stemn, rootn to organ
                    # Filter excess plant N to retransn pool for organ N
                    is_soybean = (ivt == ntmp_soybean || ivt == nirrig_tmp_soybean ||
                                  ivt == ntrp_soybean || ivt == nirrig_trp_soybean)

                    if cnveg_state.astem_patch[p] == pftcon.astemf[ivt] || !is_soybean
                        if cnveg_state.grain_flag_patch[p] == 0.0
                            if !use_fun
                                t1 = 1.0 / dt
                                cnveg_nf.leafn_to_retransn_patch[p] = t1 * (
                                    (cnveg_cs.leafc_patch[p] / pftcon.leafcn[ivt]) -
                                    (cnveg_cs.leafc_patch[p] / pftcon.fleafcn[ivt]))
                                cnveg_nf.livestemn_to_retransn_patch[p] = t1 * (
                                    (cnveg_cs.livestemc_patch[p] / pftcon.livewdcn[ivt]) -
                                    (cnveg_cs.livestemc_patch[p] / pftcon.fstemcn[ivt]))
                                cnveg_nf.frootn_to_retransn_patch[p] = 0.0
                                if pftcon.ffrootcn[ivt] > 0.0
                                    cnveg_nf.frootn_to_retransn_patch[p] = t1 * (
                                        (cnveg_cs.frootc_patch[p] / pftcon.frootcn[ivt]) -
                                        (cnveg_cs.frootc_patch[p] / pftcon.ffrootcn[ivt]))
                                end
                            else  # leafn retrans flux is handled in phenology
                                cnveg_nf.frootn_to_retransn_patch[p] = 0.0
                                cnveg_nf.livestemn_to_retransn_patch[p] = 0.0
                            end
                            cnveg_state.grain_flag_patch[p] = 1.0
                        end
                    end
                end
            end
        end
    end

    # Beth's code: crops pull from retransn pool only during grain fill;
    # retransn pool has N from leaves, stems, and roots for retranslocation
    if !use_fun
        if call_is_for_pcrop
            for p in bounds
                mask_p[p] || continue

                if cnveg_state.grain_flag_patch[p] == 1.0
                    cnveg_nf.avail_retransn_patch[p] = cnveg_nf.plant_ndemand_patch[p]
                else
                    cnveg_nf.avail_retransn_patch[p] = 0.0
                end
            end
        else
            for p in bounds
                mask_p[p] || continue

                if cnveg_state.annsum_potential_gpp_patch[p] > 0.0
                    cnveg_nf.avail_retransn_patch[p] =
                        (cnveg_state.annmax_retransn_patch[p] / 2.0) *
                        (cnveg_cf.gpp_before_downreg_patch[p] /
                         cnveg_state.annsum_potential_gpp_patch[p]) / dt
                else
                    cnveg_nf.avail_retransn_patch[p] = 0.0
                end
            end
        end

        for p in bounds
            mask_p[p] || continue

            ivt = patch.itype[p] + 1  # 0-based Fortran → 1-based Julia

            # make sure available retrans N doesn't exceed storage
            cnveg_nf.avail_retransn_patch[p] = min(cnveg_nf.avail_retransn_patch[p],
                                                     cnveg_ns.retransn_patch[p] / dt)

            # modify plant N demand according to the availability of retranslocated N
            # take from retransn pool at most the flux required to meet plant ndemand
            if cnveg_nf.plant_ndemand_patch[p] > cnveg_nf.avail_retransn_patch[p]
                cnveg_nf.retransn_to_npool_patch[p] = cnveg_nf.avail_retransn_patch[p]
            else
                cnveg_nf.retransn_to_npool_patch[p] = cnveg_nf.plant_ndemand_patch[p]
            end

            if !use_fun
                cnveg_nf.plant_ndemand_patch[p] = cnveg_nf.plant_ndemand_patch[p] -
                    cnveg_nf.retransn_to_npool_patch[p]
            else
                if pftcon.season_decid[ivt] == 1.0 || pftcon.stress_decid[ivt] == 1.0
                    cnveg_nf.plant_ndemand_patch[p] = cnveg_nf.plant_ndemand_patch[p] -
                        cnveg_nf.retransn_to_npool_patch[p]
                end
            end
        end
    end  # use_fun

    return nothing
end
