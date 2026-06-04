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
Base.@kwdef mutable struct PftConNutrientCompetition{FT<:Real, V<:AbstractVector{FT}}
    woody       ::V = Float64[]   # binary woody flag (1=woody, 0=not woody)
    froot_leaf  ::V = Float64[]   # new fine root C per new leaf C (gC/gC)
    croot_stem  ::V = Float64[]   # new coarse root C per new stem C (gC/gC)
    stem_leaf   ::V = Float64[]   # new stem C per new leaf C (gC/gC); -1 = dynamic
    flivewd     ::V = Float64[]   # fraction of new wood that is live
    leafcn      ::V = Float64[]   # leaf C:N (gC/gN)
    frootcn     ::V = Float64[]   # fine root C:N (gC/gN)
    livewdcn    ::V = Float64[]   # live wood C:N (gC/gN)
    deadwdcn    ::V = Float64[]   # dead wood C:N (gC/gN)
    fcur        ::V = Float64[]   # fraction of allocation to current growth
    graincn     ::V = Float64[]   # grain C:N (gC/gN)
    grperc      ::V = Float64[]   # growth respiration fraction
    grpnow      ::V = Float64[]   # growth respiration fraction released immediately
    # Crop grain-fill C:N override parameters
    fleafcn     ::V = Float64[]   # leaf C:N during organ fill
    ffrootcn    ::V = Float64[]   # fine root C:N during organ fill
    fstemcn     ::V = Float64[]   # stem C:N during organ fill
    astemf      ::V = Float64[]   # final stem allocation coefficient
    # Deciduous flags
    season_decid ::V = Float64[]  # binary flag for seasonal-deciduous (0 or 1)
    stress_decid ::V = Float64[]  # binary flag for stress-deciduous (0 or 1)
end
PftConNutrientCompetition{FT}(; kwargs...) where {FT<:Real} = PftConNutrientCompetition{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure PftConNutrientCompetition

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
function calc_plant_nutrient_demand!(mask_p::AbstractVector{Bool}, bounds::UnitRange{Int},
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
        dt::Real,
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
function calc_plant_nutrient_competition!(mask_soilp::AbstractVector{Bool}, bounds::UnitRange{Int},
        pftcon::PftConNutrientCompetition,
        cn_shared_params::CNSharedParamsData,
        patch::PatchData,
        crop::CropData,
        cnveg_state::CNVegStateData,
        cnveg_cs::CNVegCarbonStateData,
        cnveg_cf::CNVegCarbonFluxData,
        cnveg_nf::CNVegNitrogenFluxData;
        fpg_col::AbstractVector{<:Real},
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
# --------------------------------------------------------------------------
# Per-patch C/N allocation kernel (everything except the optional C13/C14
# downregulation writes, which are handled by _cnalloc_iso_kernel! below).
# Loop-carried scalars f5[k] are computed inline per-k (no per-patch alloc).
# --------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Device-view bundles for _cnalloc_main_kernel! (Metal ~31-arg limit).
# V = vector field type, M = matrix field type ([p,k] reproductive arrays).
# Each struct is ONE kernel arg and is @adapt_structure'd for device movement.
# ---------------------------------------------------------------------------
Base.@kwdef struct _CnAllocOut{V,M}
    sminn_to_npool::V; plant_nalloc::V; plant_calloc::V; excess_cflux::V; downreg::V
    psnsun_to_cpool::V; psnshade_to_cpool::V
    cpool_to_leafc::V; cpool_to_leafc_storage::V
    cpool_to_frootc::V; cpool_to_frootc_storage::V
    cpool_to_livestemc::V; cpool_to_livestemc_storage::V
    cpool_to_deadstemc::V; cpool_to_deadstemc_storage::V
    cpool_to_livecrootc::V; cpool_to_livecrootc_storage::V
    cpool_to_deadcrootc::V; cpool_to_deadcrootc_storage::V
    cpool_to_reproductivec::M; cpool_to_reproductivec_storage::M
    npool_to_leafn::V; npool_to_leafn_storage::V
    npool_to_frootn::V; npool_to_frootn_storage::V
    npool_to_livestemn::V; npool_to_livestemn_storage::V
    npool_to_deadstemn::V; npool_to_deadstemn_storage::V
    npool_to_livecrootn::V; npool_to_livecrootn_storage::V
    npool_to_deadcrootn::V; npool_to_deadcrootn_storage::V
    npool_to_reproductiven::M; npool_to_reproductiven_storage::M
    cpool_to_gresp_storage::V
end
Adapt.@adapt_structure _CnAllocOut

Base.@kwdef struct _CnAllocIn{V,M,VB}
    froot_leaf::V; croot_stem::V; stem_leaf::V
    annsum_npp::V; flivewd::V; grperc::V; grpnow::V
    leafcn::V; frootcn::V; livewdcn::V; deadwdcn::V
    fcur_arr::V; graincn::V; woody::V
    croplive::VB; aleaf::V; aroot::V; astem::V
    arepr::M
    sminn_to_plant_fun::V; plant_ndemand::V; fpg_col::V
    retransn_to_npool::V; c_allometry::V; n_allometry::V
    availc::V; gpp_before_downreg::V
end
Adapt.@adapt_structure _CnAllocIn

@kernel function _cnalloc_main_kernel!(
        out::_CnAllocOut, in::_CnAllocIn,
        @Const(mask_soilp), @Const(column), @Const(itype),
        use_fun::Bool, npcropmin::Int, nrepr::Int)
    T = eltype(out.plant_calloc)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        # alias every output struct field to its original loose name
        sminn_to_npool                 = out.sminn_to_npool
        plant_nalloc                   = out.plant_nalloc
        plant_calloc                   = out.plant_calloc
        excess_cflux                   = out.excess_cflux
        downreg                        = out.downreg
        psnsun_to_cpool                = out.psnsun_to_cpool
        psnshade_to_cpool              = out.psnshade_to_cpool
        cpool_to_leafc                 = out.cpool_to_leafc
        cpool_to_leafc_storage         = out.cpool_to_leafc_storage
        cpool_to_frootc                = out.cpool_to_frootc
        cpool_to_frootc_storage        = out.cpool_to_frootc_storage
        cpool_to_livestemc             = out.cpool_to_livestemc
        cpool_to_livestemc_storage     = out.cpool_to_livestemc_storage
        cpool_to_deadstemc             = out.cpool_to_deadstemc
        cpool_to_deadstemc_storage     = out.cpool_to_deadstemc_storage
        cpool_to_livecrootc            = out.cpool_to_livecrootc
        cpool_to_livecrootc_storage    = out.cpool_to_livecrootc_storage
        cpool_to_deadcrootc            = out.cpool_to_deadcrootc
        cpool_to_deadcrootc_storage    = out.cpool_to_deadcrootc_storage
        cpool_to_reproductivec         = out.cpool_to_reproductivec
        cpool_to_reproductivec_storage = out.cpool_to_reproductivec_storage
        npool_to_leafn                 = out.npool_to_leafn
        npool_to_leafn_storage         = out.npool_to_leafn_storage
        npool_to_frootn                = out.npool_to_frootn
        npool_to_frootn_storage        = out.npool_to_frootn_storage
        npool_to_livestemn             = out.npool_to_livestemn
        npool_to_livestemn_storage     = out.npool_to_livestemn_storage
        npool_to_deadstemn             = out.npool_to_deadstemn
        npool_to_deadstemn_storage     = out.npool_to_deadstemn_storage
        npool_to_livecrootn            = out.npool_to_livecrootn
        npool_to_livecrootn_storage    = out.npool_to_livecrootn_storage
        npool_to_deadcrootn            = out.npool_to_deadcrootn
        npool_to_deadcrootn_storage    = out.npool_to_deadcrootn_storage
        npool_to_reproductiven         = out.npool_to_reproductiven
        npool_to_reproductiven_storage = out.npool_to_reproductiven_storage
        cpool_to_gresp_storage         = out.cpool_to_gresp_storage
        # alias every input struct field to its original loose name
        froot_leaf         = in.froot_leaf
        croot_stem         = in.croot_stem
        stem_leaf          = in.stem_leaf
        annsum_npp         = in.annsum_npp
        flivewd            = in.flivewd
        grperc             = in.grperc
        grpnow             = in.grpnow
        leafcn             = in.leafcn
        frootcn            = in.frootcn
        livewdcn           = in.livewdcn
        deadwdcn           = in.deadwdcn
        fcur_arr           = in.fcur_arr
        graincn            = in.graincn
        woody              = in.woody
        croplive           = in.croplive
        aleaf              = in.aleaf
        aroot              = in.aroot
        astem              = in.astem
        arepr              = in.arepr
        sminn_to_plant_fun = in.sminn_to_plant_fun
        plant_ndemand      = in.plant_ndemand
        fpg_col            = in.fpg_col
        retransn_to_npool  = in.retransn_to_npool
        c_allometry        = in.c_allometry
        n_allometry        = in.n_allometry
        availc             = in.availc
        gpp_before_downreg = in.gpp_before_downreg

        c = column[p]
        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia

        # set local allocation variables
        f1 = froot_leaf[ivt]
        f2 = croot_stem[ivt]

        # modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
        # constrained so that it does not go lower than 0.2 (under negative annsum_npp)
        # -1.0 means dynamic allocation (trees)
        if stem_leaf[ivt] == T(-1.0)
            f3 = (T(2.7) / (one(T) + exp(T(-0.004) * (annsum_npp[p] - T(300.0))))) - T(0.4)
        else
            f3 = stem_leaf[ivt]
        end

        f4   = flivewd[ivt]
        g1   = grperc[ivt]
        g2   = grpnow[ivt]
        cnl  = leafcn[ivt]
        cnfr = frootcn[ivt]
        cnlw = livewdcn[ivt]
        cndw = deadwdcn[ivt]
        fcur = fcur_arr[ivt]

        is_crop = ivt >= npcropmin
        crop_alive = is_crop && (croplive[p] != zero(eltype(croplive))) && !isnan(aleaf[p])

        if is_crop  # skip 2 generic crops
            if crop_alive
                f1 = aroot[p] / aleaf[p]
                f3 = astem[p] / aleaf[p]
                g1 = T(0.25)
            else
                f1 = zero(T)
                f3 = zero(T)
                g1 = T(0.25)
            end
        end

        if use_fun  # if we are using FUN, we get the N available from there
            sminn_to_npool[p] = sminn_to_plant_fun[p]
        else  # no FUN - get N available from the FPG calculation
            sminn_to_npool[p] = plant_ndemand[p] * fpg_col[c]
        end

        plant_nalloc[p] = sminn_to_npool[p] + retransn_to_npool[p]

        plant_calloc[p] = plant_nalloc[p] * (c_allometry[p] / n_allometry[p])

        if !use_fun  # ORIGINAL CLM(CN) downregulation code
            excess_cflux[p] = availc[p] - plant_calloc[p]

            # reduce gpp fluxes due to N limitation
            if gpp_before_downreg[p] > zero(T)
                downreg[p] = excess_cflux[p] / gpp_before_downreg[p]

                psnsun_to_cpool[p]   = psnsun_to_cpool[p]   * (one(T) - downreg[p])
                psnshade_to_cpool[p] = psnshade_to_cpool[p] * (one(T) - downreg[p])
            end
        end  # use_fun

        # calculate new leaf C and daily fluxes to current growth and storage pools
        nlc = plant_calloc[p] / c_allometry[p]

        cpool_to_leafc[p]          = nlc * fcur
        cpool_to_leafc_storage[p]  = nlc * (one(T) - fcur)
        cpool_to_frootc[p]         = nlc * f1 * fcur
        cpool_to_frootc_storage[p] = nlc * f1 * (one(T) - fcur)

        if woody[ivt] == one(T)
            cpool_to_livestemc[p]          = nlc * f3 * f4 * fcur
            cpool_to_livestemc_storage[p]  = nlc * f3 * f4 * (one(T) - fcur)
            cpool_to_deadstemc[p]          = nlc * f3 * (one(T) - f4) * fcur
            cpool_to_deadstemc_storage[p]  = nlc * f3 * (one(T) - f4) * (one(T) - fcur)
            cpool_to_livecrootc[p]         = nlc * f2 * f3 * f4 * fcur
            cpool_to_livecrootc_storage[p] = nlc * f2 * f3 * f4 * (one(T) - fcur)
            cpool_to_deadcrootc[p]         = nlc * f2 * f3 * (one(T) - f4) * fcur
            cpool_to_deadcrootc_storage[p] = nlc * f2 * f3 * (one(T) - f4) * (one(T) - fcur)
        end

        if is_crop  # skip 2 generic crops
            cpool_to_livestemc[p]          = nlc * f3 * f4 * fcur
            cpool_to_livestemc_storage[p]  = nlc * f3 * f4 * (one(T) - fcur)
            cpool_to_deadstemc[p]          = nlc * f3 * (one(T) - f4) * fcur
            cpool_to_deadstemc_storage[p]  = nlc * f3 * (one(T) - f4) * (one(T) - fcur)
            cpool_to_livecrootc[p]         = nlc * f2 * f3 * f4 * fcur
            cpool_to_livecrootc_storage[p] = nlc * f2 * f3 * f4 * (one(T) - fcur)
            cpool_to_deadcrootc[p]         = nlc * f2 * f3 * (one(T) - f4) * fcur
            cpool_to_deadcrootc_storage[p] = nlc * f2 * f3 * (one(T) - f4) * (one(T) - fcur)
            for k in 1:nrepr
                f5k = crop_alive ? (arepr[p, k] / aleaf[p]) : zero(T)
                cpool_to_reproductivec[p, k]         = nlc * f5k * fcur
                cpool_to_reproductivec_storage[p, k] = nlc * f5k * (one(T) - fcur)
            end
        end

        # corresponding N fluxes
        npool_to_leafn[p]          = (nlc / cnl) * fcur
        npool_to_leafn_storage[p]  = (nlc / cnl) * (one(T) - fcur)
        npool_to_frootn[p]         = (nlc * f1 / cnfr) * fcur
        npool_to_frootn_storage[p] = (nlc * f1 / cnfr) * (one(T) - fcur)

        if woody[ivt] == one(T)
            npool_to_livestemn[p]          = (nlc * f3 * f4 / cnlw) * fcur
            npool_to_livestemn_storage[p]  = (nlc * f3 * f4 / cnlw) * (one(T) - fcur)
            npool_to_deadstemn[p]          = (nlc * f3 * (one(T) - f4) / cndw) * fcur
            npool_to_deadstemn_storage[p]  = (nlc * f3 * (one(T) - f4) / cndw) * (one(T) - fcur)
            npool_to_livecrootn[p]         = (nlc * f2 * f3 * f4 / cnlw) * fcur
            npool_to_livecrootn_storage[p] = (nlc * f2 * f3 * f4 / cnlw) * (one(T) - fcur)
            npool_to_deadcrootn[p]         = (nlc * f2 * f3 * (one(T) - f4) / cndw) * fcur
            npool_to_deadcrootn_storage[p] = (nlc * f2 * f3 * (one(T) - f4) / cndw) * (one(T) - fcur)
        end

        if is_crop  # skip 2 generic crops
            cng = graincn[ivt]
            npool_to_livestemn[p]          = (nlc * f3 * f4 / cnlw) * fcur
            npool_to_livestemn_storage[p]  = (nlc * f3 * f4 / cnlw) * (one(T) - fcur)
            npool_to_deadstemn[p]          = (nlc * f3 * (one(T) - f4) / cndw) * fcur
            npool_to_deadstemn_storage[p]  = (nlc * f3 * (one(T) - f4) / cndw) * (one(T) - fcur)
            npool_to_livecrootn[p]         = (nlc * f2 * f3 * f4 / cnlw) * fcur
            npool_to_livecrootn_storage[p] = (nlc * f2 * f3 * f4 / cnlw) * (one(T) - fcur)
            npool_to_deadcrootn[p]         = (nlc * f2 * f3 * (one(T) - f4) / cndw) * fcur
            npool_to_deadcrootn_storage[p] = (nlc * f2 * f3 * (one(T) - f4) / cndw) * (one(T) - fcur)
            for k in 1:nrepr
                f5k = crop_alive ? (arepr[p, k] / aleaf[p]) : zero(T)
                npool_to_reproductiven[p, k]         = (nlc * f5k / cng) * fcur
                npool_to_reproductiven_storage[p, k] = (nlc * f5k / cng) * (one(T) - fcur)
            end
        end

        # Growth respiration storage: carbon that needs to go into growth
        # respiration storage to satisfy all of the storage growth demands
        gresp_storage = cpool_to_leafc_storage[p] + cpool_to_frootc_storage[p]
        if woody[ivt] == one(T)
            gresp_storage += cpool_to_livestemc_storage[p]
            gresp_storage += cpool_to_deadstemc_storage[p]
            gresp_storage += cpool_to_livecrootc_storage[p]
            gresp_storage += cpool_to_deadcrootc_storage[p]
        end
        if is_crop  # skip 2 generic crops
            gresp_storage += cpool_to_livestemc_storage[p]
            for k in 1:nrepr
                gresp_storage += cpool_to_reproductivec_storage[p, k]
            end
        end
        cpool_to_gresp_storage[p] = gresp_storage * g1 * (one(T) - g2)
    end
end

# Separate kernel for the optional isotope (C13 or C14) downregulation writes.
# Launched only when the corresponding flag is true AND the flux struct exists,
# so a `nothing` struct is never read on the device. Re-reads downreg + the same
# guards (!use_fun, gpp>0) the main kernel used so the multiply matches exactly.
@kernel function _cnalloc_iso_kernel!(
        iso_psnsun_to_cpool, iso_psnshade_to_cpool,
        @Const(mask_soilp), @Const(downreg), @Const(gpp_before_downreg),
        use_fun::Bool)
    T = eltype(iso_psnsun_to_cpool)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        if !use_fun
            if gpp_before_downreg[p] > zero(T)
                iso_psnsun_to_cpool[p]   = iso_psnsun_to_cpool[p]   * (one(T) - downreg[p])
                iso_psnshade_to_cpool[p] = iso_psnshade_to_cpool[p] * (one(T) - downreg[p])
            end
        end
    end
end

function calc_plant_cn_alloc!(mask_soilp::AbstractVector{Bool}, bounds::UnitRange{Int},
        pftcon::PftConNutrientCompetition,
        cn_shared_params::CNSharedParamsData,
        patch::PatchData,
        crop::CropData,
        cnveg_state::CNVegStateData,
        cnveg_cs::CNVegCarbonStateData,
        cnveg_cf::CNVegCarbonFluxData,
        cnveg_nf::CNVegNitrogenFluxData;
        fpg_col::AbstractVector{<:Real},
        c13_cnveg_cf::Union{CNVegCarbonFluxData,Nothing}=nothing,
        c14_cnveg_cf::Union{CNVegCarbonFluxData,Nothing}=nothing,
        use_c13::Bool=false,
        use_c14::Bool=false,
        npcropmin::Int=NPCROPMIN,
        nrepr::Int=NREPR)

    use_fun = cn_shared_params.use_fun

    out = _CnAllocOut(;
        sminn_to_npool                 = cnveg_nf.sminn_to_npool_patch,
        plant_nalloc                   = cnveg_nf.plant_nalloc_patch,
        plant_calloc                   = cnveg_cf.plant_calloc_patch,
        excess_cflux                   = cnveg_cf.excess_cflux_patch,
        downreg                        = cnveg_state.downreg_patch,
        psnsun_to_cpool                = cnveg_cf.psnsun_to_cpool_patch,
        psnshade_to_cpool              = cnveg_cf.psnshade_to_cpool_patch,
        cpool_to_leafc                 = cnveg_cf.cpool_to_leafc_patch,
        cpool_to_leafc_storage         = cnveg_cf.cpool_to_leafc_storage_patch,
        cpool_to_frootc                = cnveg_cf.cpool_to_frootc_patch,
        cpool_to_frootc_storage        = cnveg_cf.cpool_to_frootc_storage_patch,
        cpool_to_livestemc             = cnveg_cf.cpool_to_livestemc_patch,
        cpool_to_livestemc_storage     = cnveg_cf.cpool_to_livestemc_storage_patch,
        cpool_to_deadstemc             = cnveg_cf.cpool_to_deadstemc_patch,
        cpool_to_deadstemc_storage     = cnveg_cf.cpool_to_deadstemc_storage_patch,
        cpool_to_livecrootc            = cnveg_cf.cpool_to_livecrootc_patch,
        cpool_to_livecrootc_storage    = cnveg_cf.cpool_to_livecrootc_storage_patch,
        cpool_to_deadcrootc            = cnveg_cf.cpool_to_deadcrootc_patch,
        cpool_to_deadcrootc_storage    = cnveg_cf.cpool_to_deadcrootc_storage_patch,
        cpool_to_reproductivec         = cnveg_cf.cpool_to_reproductivec_patch,
        cpool_to_reproductivec_storage = cnveg_cf.cpool_to_reproductivec_storage_patch,
        npool_to_leafn                 = cnveg_nf.npool_to_leafn_patch,
        npool_to_leafn_storage         = cnveg_nf.npool_to_leafn_storage_patch,
        npool_to_frootn                = cnveg_nf.npool_to_frootn_patch,
        npool_to_frootn_storage        = cnveg_nf.npool_to_frootn_storage_patch,
        npool_to_livestemn             = cnveg_nf.npool_to_livestemn_patch,
        npool_to_livestemn_storage     = cnveg_nf.npool_to_livestemn_storage_patch,
        npool_to_deadstemn             = cnveg_nf.npool_to_deadstemn_patch,
        npool_to_deadstemn_storage     = cnveg_nf.npool_to_deadstemn_storage_patch,
        npool_to_livecrootn            = cnveg_nf.npool_to_livecrootn_patch,
        npool_to_livecrootn_storage    = cnveg_nf.npool_to_livecrootn_storage_patch,
        npool_to_deadcrootn            = cnveg_nf.npool_to_deadcrootn_patch,
        npool_to_deadcrootn_storage    = cnveg_nf.npool_to_deadcrootn_storage_patch,
        npool_to_reproductiven         = cnveg_nf.npool_to_reproductiven_patch,
        npool_to_reproductiven_storage = cnveg_nf.npool_to_reproductiven_storage_patch,
        cpool_to_gresp_storage         = cnveg_cf.cpool_to_gresp_storage_patch)

    in = _CnAllocIn(;
        froot_leaf         = pftcon.froot_leaf,
        croot_stem         = pftcon.croot_stem,
        stem_leaf          = pftcon.stem_leaf,
        annsum_npp         = cnveg_cf.annsum_npp_patch,
        flivewd            = pftcon.flivewd,
        grperc             = pftcon.grperc,
        grpnow             = pftcon.grpnow,
        leafcn             = pftcon.leafcn,
        frootcn            = pftcon.frootcn,
        livewdcn           = pftcon.livewdcn,
        deadwdcn           = pftcon.deadwdcn,
        fcur_arr           = pftcon.fcur,
        graincn            = pftcon.graincn,
        woody              = pftcon.woody,
        croplive           = crop.croplive_patch,
        aleaf              = cnveg_state.aleaf_patch,
        aroot              = cnveg_state.aroot_patch,
        astem              = cnveg_state.astem_patch,
        arepr              = cnveg_state.arepr_patch,
        sminn_to_plant_fun = cnveg_nf.sminn_to_plant_fun_patch,
        plant_ndemand      = cnveg_nf.plant_ndemand_patch,
        fpg_col            = fpg_col,
        retransn_to_npool  = cnveg_nf.retransn_to_npool_patch,
        c_allometry        = cnveg_state.c_allometry_patch,
        n_allometry        = cnveg_state.n_allometry_patch,
        availc             = cnveg_cf.availc_patch,
        gpp_before_downreg = cnveg_cf.gpp_before_downreg_patch)

    # mask + integer index vectors onto the state backend (BitVector / host
    # non-bitstype on device). Reference array carries the backend + eltype.
    ref = out.plant_calloc
    FT = eltype(ref)
    nd = length(mask_soilp)
    mask_d = similar(ref, Bool, nd); copyto!(mask_d, collect(Bool, mask_soilp))
    column_d = similar(ref, Int, length(patch.column)); copyto!(column_d, patch.column)
    itype_d = similar(ref, Int, length(patch.itype)); copyto!(itype_d, patch.itype)

    # Struct-first kernel: manual backend launch + synchronize (struct args carry
    # no backend, so we take it from a known device array field).
    backend = _kernel_backend(ref)
    _cnalloc_main_kernel!(backend)(out, in, mask_d, column_d, itype_d,
        use_fun, npcropmin, nrepr; ndrange = nd)
    KA.synchronize(backend)

    if use_c13 && c13_cnveg_cf !== nothing
        _launch!(_cnalloc_iso_kernel!,
            c13_cnveg_cf.psnsun_to_cpool_patch, c13_cnveg_cf.psnshade_to_cpool_patch,
            mask_soilp, cnveg_state.downreg_patch, cnveg_cf.gpp_before_downreg_patch,
            use_fun)
    end
    if use_c14 && c14_cnveg_cf !== nothing
        _launch!(_cnalloc_iso_kernel!,
            c14_cnveg_cf.psnsun_to_cpool_patch, c14_cnveg_cf.psnshade_to_cpool_patch,
            mask_soilp, cnveg_state.downreg_patch, cnveg_cf.gpp_before_downreg_patch,
            use_fun)
    end

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
# --- Loop 1: total plant N demand + per-patch tempsum/tempmax accumulation -
# Each patch updates its OWN tempsum/tempmax accumulator (own-index read-
# modify-write) — fully parallel, NOT a cross-patch reduction.
@kernel function _npdemand_demand_kernel!(plant_ndemand, tempsum_potential_gpp,
        tempmax_retransn,
        @Const(mask), @Const(availc), @Const(n_allometry), @Const(c_allometry),
        @Const(gpp_before_downreg), @Const(retransn))
    p = @index(Global)
    @inbounds if mask[p]
        plant_ndemand[p] = availc[p] * (n_allometry[p] / c_allometry[p])

        # retranslocated N deployment depends on seasonal cycle of potential GPP
        tempsum_potential_gpp[p] = tempsum_potential_gpp[p] + gpp_before_downreg[p]

        # carry max retransn info to CN Annual Update
        tempmax_retransn[p] = max(tempmax_retransn[p], retransn[p])
    end
end

# --- Loop 2: crop grain-fill retranslocation state machine ----------------
# Per-patch grain-fill trigger; is_soybean derived per-patch from ivt. The
# crop_phase values are computed on the host (crop_phase!) and read here.
@kernel function _npdemand_grainfill_kernel!(grain_flag,
        leafn_to_retransn, livestemn_to_retransn, frootn_to_retransn,
        @Const(mask), @Const(itype), @Const(croplive), @Const(crop_phase_vals),
        @Const(astem), @Const(astemf), @Const(leafc), @Const(frootc),
        @Const(livestemc), @Const(leafcn), @Const(fleafcn), @Const(livewdcn),
        @Const(fstemcn), @Const(frootcn), @Const(ffrootcn),
        cphase_leafemerge_in, cphase_grainfill_in,
        ntmp_soybean::Int, nirrig_tmp_soybean::Int,
        ntrp_soybean::Int, nirrig_trp_soybean::Int,
        dt, use_fun::Bool)
    p = @index(Global)
    @inbounds if mask[p]
        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia

        if croplive[p]
            if crop_phase_vals[p] == cphase_leafemerge_in
                grain_flag[p] = zero(eltype(grain_flag))  # setting to 0 while in phase 2
            elseif crop_phase_vals[p] == cphase_grainfill_in
                # Beth's retranslocation of leafn, stemn, rootn to organ
                is_soybean = (ivt == ntmp_soybean || ivt == nirrig_tmp_soybean ||
                              ivt == ntrp_soybean || ivt == nirrig_trp_soybean)

                if astem[p] == astemf[ivt] || !is_soybean
                    if grain_flag[p] == zero(eltype(grain_flag))
                        if !use_fun
                            t1 = one(dt) / dt
                            leafn_to_retransn[p] = t1 * (
                                (leafc[p] / leafcn[ivt]) -
                                (leafc[p] / fleafcn[ivt]))
                            livestemn_to_retransn[p] = t1 * (
                                (livestemc[p] / livewdcn[ivt]) -
                                (livestemc[p] / fstemcn[ivt]))
                            frootn_to_retransn[p] = zero(eltype(frootn_to_retransn))
                            if ffrootcn[ivt] > zero(eltype(ffrootcn))
                                frootn_to_retransn[p] = t1 * (
                                    (frootc[p] / frootcn[ivt]) -
                                    (frootc[p] / ffrootcn[ivt]))
                            end
                        else  # leafn retrans flux is handled in phenology
                            frootn_to_retransn[p] = zero(eltype(frootn_to_retransn))
                            livestemn_to_retransn[p] = zero(eltype(livestemn_to_retransn))
                        end
                        grain_flag[p] = one(eltype(grain_flag))
                    end
                end
            end
        end
    end
end

# --- Loop 3a: avail_retransn for prognostic-crop call (grain-flag gated) ---
@kernel function _npdemand_availretransn_pcrop_kernel!(avail_retransn,
        @Const(mask), @Const(grain_flag), @Const(plant_ndemand))
    p = @index(Global)
    @inbounds if mask[p]
        if grain_flag[p] == one(eltype(grain_flag))
            avail_retransn[p] = plant_ndemand[p]
        else
            avail_retransn[p] = zero(eltype(avail_retransn))
        end
    end
end

# --- Loop 3b: avail_retransn for non-pcrop call (seasonal GPP fraction) ----
@kernel function _npdemand_availretransn_nopcrop_kernel!(avail_retransn,
        @Const(mask), @Const(annsum_potential_gpp), @Const(annmax_retransn),
        @Const(gpp_before_downreg), dt)
    T = eltype(avail_retransn)
    p = @index(Global)
    @inbounds if mask[p]
        if annsum_potential_gpp[p] > zero(T)
            avail_retransn[p] =
                (annmax_retransn[p] / T(2.0)) *
                (gpp_before_downreg[p] / annsum_potential_gpp[p]) / dt
        else
            avail_retransn[p] = zero(T)
        end
    end
end

# --- Loop 4: clamp avail_retransn to storage + modify plant N demand -------
@kernel function _npdemand_retransn_modify_kernel!(avail_retransn,
        retransn_to_npool, plant_ndemand,
        @Const(mask), @Const(itype), @Const(retransn),
        @Const(season_decid), @Const(stress_decid), dt, use_fun::Bool)
    p = @index(Global)
    @inbounds if mask[p]
        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia

        # make sure available retrans N doesn't exceed storage
        avail_retransn[p] = min(avail_retransn[p], retransn[p] / dt)

        # take from retransn pool at most the flux required to meet plant ndemand
        if plant_ndemand[p] > avail_retransn[p]
            retransn_to_npool[p] = avail_retransn[p]
        else
            retransn_to_npool[p] = plant_ndemand[p]
        end

        if !use_fun
            plant_ndemand[p] = plant_ndemand[p] - retransn_to_npool[p]
        else
            if season_decid[ivt] == one(eltype(season_decid)) ||
               stress_decid[ivt] == one(eltype(stress_decid))
                plant_ndemand[p] = plant_ndemand[p] - retransn_to_npool[p]
            end
        end
    end
end

function calc_plant_nitrogen_demand!(mask_p::AbstractVector{Bool}, bounds::UnitRange{Int},
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
        dt::Real,
        npcropmin::Int=NPCROPMIN,
        nrepr::Int=NREPR,
        ntmp_soybean::Int=23,
        nirrig_tmp_soybean::Int=24,
        ntrp_soybean::Int=77,
        nirrig_trp_soybean::Int=78)

    use_fun = cn_shared_params.use_fun

    # loop over patches to assess the total plant N demand
    _launch!(_npdemand_demand_kernel!, cnveg_nf.plant_ndemand_patch,
        cnveg_state.tempsum_potential_gpp_patch,
        cnveg_state.tempmax_retransn_patch,
        mask_p, cnveg_cf.availc_patch, cnveg_state.n_allometry_patch,
        cnveg_state.c_allometry_patch, cnveg_cf.gpp_before_downreg_patch,
        cnveg_ns.retransn_patch)

    # Crop grain-fill retranslocation
    if call_is_for_pcrop
        # Compute crop phase for all patches in bounds (host-side helper,
        # already GPU-capable; the kernel below reads its output array).
        FT = eltype(cnveg_nf.plant_ndemand_patch)
        # device-resident scratch via similar() (a bare zeros(...) is host →
        # would force the CPU backend for crop_phase! and the grainfill kernel).
        crop_phase_vals = fill!(similar(cnveg_nf.plant_ndemand_patch, FT, length(mask_p)), zero(FT))
        crop_phase!(mask_p, crop, cnveg_state, crop_phase_vals)

        _launch!(_npdemand_grainfill_kernel!, cnveg_state.grain_flag_patch,
            cnveg_nf.leafn_to_retransn_patch,
            cnveg_nf.livestemn_to_retransn_patch,
            cnveg_nf.frootn_to_retransn_patch,
            mask_p, patch.itype, crop.croplive_patch, crop_phase_vals,
            cnveg_state.astem_patch, pftcon.astemf,
            cnveg_cs.leafc_patch, cnveg_cs.frootc_patch, cnveg_cs.livestemc_patch,
            pftcon.leafcn, pftcon.fleafcn, pftcon.livewdcn, pftcon.fstemcn,
            pftcon.frootcn, pftcon.ffrootcn,
            FT(cphase_leafemerge), FT(cphase_grainfill),
            ntmp_soybean, nirrig_tmp_soybean, ntrp_soybean, nirrig_trp_soybean,
            FT(dt), use_fun)
    end

    # Beth's code: crops pull from retransn pool only during grain fill;
    # retransn pool has N from leaves, stems, and roots for retranslocation
    if !use_fun
        FT = eltype(cnveg_nf.avail_retransn_patch)
        if call_is_for_pcrop
            _launch!(_npdemand_availretransn_pcrop_kernel!,
                cnveg_nf.avail_retransn_patch,
                mask_p, cnveg_state.grain_flag_patch,
                cnveg_nf.plant_ndemand_patch)
        else
            _launch!(_npdemand_availretransn_nopcrop_kernel!,
                cnveg_nf.avail_retransn_patch,
                mask_p, cnveg_state.annsum_potential_gpp_patch,
                cnveg_state.annmax_retransn_patch,
                cnveg_cf.gpp_before_downreg_patch, FT(dt))
        end

        _launch!(_npdemand_retransn_modify_kernel!,
            cnveg_nf.avail_retransn_patch,
            cnveg_nf.retransn_to_npool_patch,
            cnveg_nf.plant_ndemand_patch,
            mask_p, patch.itype, cnveg_ns.retransn_patch,
            pftcon.season_decid, pftcon.stress_decid, FT(dt), use_fun)
    end  # use_fun

    return nothing
end
