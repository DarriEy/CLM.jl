# ==========================================================================
# Ported from: src/biogeochem/NutrientCompetitionFlexibleCNMod.F90
# Plant nutrient demand and competition algorithm (CLM5 FlexibleCN method)
#
# Original description:
#   Module contains different subroutines to do soil nutrient competition
#   dynamics. This module was copied from NutrientCompetitionCLM45default
#   then flexible-CN modifications were added for the clm50 nitrogen science
#   changes (r120).
#   created by Jinyun Tang, Sep 8, 2014
#   modified by Mariana Vertenstein, Nov 15, 2014
#
# This file is the parallel of `nutrient_competition.jl` (the CLM4.5 default
# method). It is dispatched on `use_flexiblecn` from the CN driver. The two
# methods differ in three places:
#   1. calc_plant_nitrogen_demand: FlexibleCN computes plant_ndemand from
#      Michaelis-Menten kinetics on fine-root C, soil mineral N substrate, a
#      soil-temperature scalar and a leaf-C:N nitrogen scalar (`nscalar`) that
#      flexes the actual leaf C:N between leafcn±10 — instead of the default's
#      availc*(n_allometry/c_allometry). The crop grain-fill retranslocation
#      uses max(leafn - leafc/fleafcn, 0) (the pool-relative form) instead of
#      the default's (leafc/leafcn - leafc/fleafcn).
#   2. calc_plant_cn_alloc: FlexibleCN sets plant_calloc = availc (non-FUN)
#      instead of plant_nalloc*(c_allometry/n_allometry), and has NO GPP
#      downregulation. The npool_to_* fluxes are drawn from the *existing*
#      npool by fractional demand (calc_npool_to_components_flexiblecn).
#   3. carbon_resp_opt==1: high actual-storage-C:N tissue turns over extra C
#      to cpool_to_*_resp (the flexible-C:N turnover term).
#
# The CLM4.5 default reproduces exactly when these branches are not taken.
#
# Two further branches, gated OFF by default and byte-identical when off:
#   * use_matrixcn (cn_shared_params.use_matrixcn): assemble the matrix-CN
#     allocation inputs — matrix_Cinput / matrix_Ninput (the C/N input rates)
#     and matrix_alloc / matrix_nalloc (the B-matrix allocation fractions per
#     veg pool), plus the C13/C14 matrix isotope inputs. The retranslocation-
#     pool transfer registrations (matrix_update_phn) are deferred pending the
#     N-matrix solve (see calc_plant_cn_alloc_flexiblecn! and cn_veg_matrix.jl).
#   * use_crop_agsys: for crop patches, draw the npool_to_* N fluxes from the
#     AgSys-supplied N allocation coefficients (aleaf_n/astem_n/aroot_n/arepr_n)
#     instead of the FlexibleCN fractional-demand draw.
# ==========================================================================

# Minimum N value used when forming a C:N ratio (CNPrecisionControlMod n_min)
const NUTRIENT_COMP_N_MIN = 1.0e-9

# ---------------------------------------------------------------------------
# Matrix-CN vegetation pool indices (clm_varpar.F90:89-109). These are the
# column indices into matrix_alloc_patch / matrix_nalloc_patch (the B-matrix
# carbon / nitrogen allocation fractions) that the use_matrixcn branch writes.
# Fortran parameters are 1-based, so they map directly to Julia 1-based.
# ---------------------------------------------------------------------------
const FLEXCN_ILEAF        = 1
const FLEXCN_ILEAF_ST     = 2
const FLEXCN_IFROOT       = 4
const FLEXCN_IFROOT_ST    = 5
const FLEXCN_ILIVESTEM    = 7
const FLEXCN_ILIVESTEM_ST = 8
const FLEXCN_IDEADSTEM    = 10
const FLEXCN_IDEADSTEM_ST = 11
const FLEXCN_ILIVECROOT   = 13
const FLEXCN_ILIVECROOT_ST = 14
const FLEXCN_IDEADCROOT   = 16
const FLEXCN_IDEADCROOT_ST = 17
const FLEXCN_IGRAIN       = 19
const FLEXCN_IGRAIN_ST    = 20

# ---------------------------------------------------------------------------
# _calc_npool_to_components_agsys! -- per-crop-patch N allocation when running
# with AgSys (the agricultural-systems managed-crop dispatch). Instead of the
# FlexibleCN fractional-demand draw, AgSys supplies its own leaf/stem/root/repr
# N allocation coefficients (aleaf_n, astem_n, aroot_n, arepr_n) and the npool
# is distributed directly in proportion to those. No allocation to coarse roots.
# Ported from `calc_npool_to_components_agsys`
# (NutrientCompetitionFlexibleCNMod.F90:1216-1302).
# ---------------------------------------------------------------------------
@inline function _calc_npool_to_components_agsys!(p::Int,
        npool::Real, fcur::Real, f4::Real,
        aleaf_n::Real, astem_n::Real, aroot_n::Real, arepr_n,
        nf::CNVegNitrogenFluxData, dt::Real, nrepr::Int)
    T = typeof(npool)
    npool_dt = npool / T(dt)

    nf.npool_to_leafn_patch[p]          = aleaf_n * fcur           * npool_dt
    nf.npool_to_leafn_storage_patch[p]  = aleaf_n * (one(T) - fcur) * npool_dt

    nf.npool_to_frootn_patch[p]         = aroot_n * fcur           * npool_dt
    nf.npool_to_frootn_storage_patch[p] = aroot_n * (one(T) - fcur) * npool_dt

    nf.npool_to_livestemn_patch[p]          = astem_n * f4           * fcur           * npool_dt
    nf.npool_to_livestemn_storage_patch[p]  = astem_n * f4           * (one(T) - fcur) * npool_dt
    nf.npool_to_deadstemn_patch[p]          = astem_n * (one(T) - f4) * fcur           * npool_dt
    nf.npool_to_deadstemn_storage_patch[p]  = astem_n * (one(T) - f4) * (one(T) - fcur) * npool_dt

    # Assume no allocation to coarse roots for crops (matches Fortran comment).
    nf.npool_to_livecrootn_patch[p]         = zero(T)
    nf.npool_to_livecrootn_storage_patch[p] = zero(T)
    nf.npool_to_deadcrootn_patch[p]         = zero(T)
    nf.npool_to_deadcrootn_storage_patch[p] = zero(T)

    for k in 1:nrepr
        nf.npool_to_reproductiven_patch[p, k]         = arepr_n[p, k] * fcur           * npool_dt
        nf.npool_to_reproductiven_storage_patch[p, k] = arepr_n[p, k] * (one(T) - fcur) * npool_dt
    end
    return nothing
end

# ---------------------------------------------------------------------------
# calc_plant_nutrient_competition_flexiblecn! -- public competition entry.
# Wraps calc_plant_cn_alloc_flexiblecn!.
# ---------------------------------------------------------------------------

"""
    calc_plant_nutrient_competition_flexiblecn!(mask_soilp, bounds,
        pftcon, cn_shared_params, patch, crop, canopystate,
        cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf,
        soilbgc_ns;
        fpg_col, c13_cnveg_cf, c14_cnveg_cf, use_c13, use_c14,
        carbon_resp_opt, dt, npcropmin, nrepr)

FlexibleCN nutrient-competition: distribute the available N between competing
patches and allocate C and N to new growth and storage using the flexible
leaf-C:N stoichiometry. Parallel of `calc_plant_nutrient_competition!`.

Ported from `calc_plant_nutrient_competition` /
`calc_plant_cn_alloc` in `NutrientCompetitionFlexibleCNMod.F90`.
"""
function calc_plant_nutrient_competition_flexiblecn!(mask_soilp::AbstractVector{Bool},
        bounds::UnitRange{Int},
        pftcon::PftConNutrientCompetition,
        cn_shared_params::CNSharedParamsData,
        patch::PatchData,
        crop::CropData,
        canopystate::CanopyStateData,
        cnveg_state::CNVegStateData,
        cnveg_cs::CNVegCarbonStateData,
        cnveg_cf::CNVegCarbonFluxData,
        cnveg_ns::CNVegNitrogenStateData,
        cnveg_nf::CNVegNitrogenFluxData;
        fpg_col::AbstractVector{<:Real},
        dt::Real=1800.0,
        c13_cnveg_cf::Union{CNVegCarbonFluxData,Nothing}=nothing,
        c14_cnveg_cf::Union{CNVegCarbonFluxData,Nothing}=nothing,
        use_c13::Bool=false,
        use_c14::Bool=false,
        carbon_resp_opt::Int=0,
        use_crop_agsys::Bool=false,
        npcropmin::Int=NPCROPMIN,
        nrepr::Int=NREPR)

    calc_plant_cn_alloc_flexiblecn!(mask_soilp, bounds,
        pftcon, cn_shared_params, patch, crop, canopystate,
        cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf;
        fpg_col=fpg_col, dt=dt,
        c13_cnveg_cf=c13_cnveg_cf, c14_cnveg_cf=c14_cnveg_cf,
        use_c13=use_c13, use_c14=use_c14,
        carbon_resp_opt=carbon_resp_opt,
        use_crop_agsys=use_crop_agsys,
        npcropmin=npcropmin, nrepr=nrepr)

    return nothing
end

# ---------------------------------------------------------------------------
# calc_npool_to_components_flexiblecn! -- per-patch N allocation from npool.
# Distributes a fraction of the EXISTING npool to each tissue in proportion to
# fractional N demand (demand_X / total_demand) * npool / dt.
# Ported from calc_npool_to_components_flexiblecn (Fortran 932-1213).
# ---------------------------------------------------------------------------
@inline function _calc_npool_to_components_flexiblecn!(p::Int, ivt::Int,
        npool::Real, nlc::Real, fcur::Real,
        f1::Real, f2::Real, f3::Real, f4::Real, arepr, aleaf_p,
        crop_alive::Bool, is_crop::Bool, is_woody::Bool,
        leafcn, frootcn, livewdcn, deadwdcn, graincn,
        nf::CNVegNitrogenFluxData, dt::Real, nrepr::Int)
    T = typeof(nlc)
    cnl  = leafcn[ivt]
    cnfr = frootcn[ivt]
    cnlw = livewdcn[ivt]
    cndw = deadwdcn[ivt]

    # 1) tissue N demands from C allocated / C:N ratio
    d_leafn          = (nlc / cnl) * fcur
    d_leafn_storage  = (nlc / cnl) * (one(T) - fcur)
    d_frootn         = (nlc * f1 / cnfr) * fcur
    d_frootn_storage = (nlc * f1 / cnfr) * (one(T) - fcur)

    d_livestemn  = zero(T); d_livestemn_storage  = zero(T)
    d_deadstemn  = zero(T); d_deadstemn_storage  = zero(T)
    d_livecrootn = zero(T); d_livecrootn_storage = zero(T)
    d_deadcrootn = zero(T); d_deadcrootn_storage = zero(T)
    if is_woody || is_crop
        d_livestemn          = (nlc * f3 * f4 / cnlw) * fcur
        d_livestemn_storage  = (nlc * f3 * f4 / cnlw) * (one(T) - fcur)
        d_deadstemn          = (nlc * f3 * (one(T) - f4) / cndw) * fcur
        d_deadstemn_storage  = (nlc * f3 * (one(T) - f4) / cndw) * (one(T) - fcur)
        d_livecrootn         = (nlc * f2 * f3 * f4 / cnlw) * fcur
        d_livecrootn_storage = (nlc * f2 * f3 * f4 / cnlw) * (one(T) - fcur)
        d_deadcrootn         = (nlc * f2 * f3 * (one(T) - f4) / cndw) * fcur
        d_deadcrootn_storage = (nlc * f2 * f3 * (one(T) - f4) / cndw) * (one(T) - fcur)
    end

    # 2) total demand
    total = d_leafn + d_leafn_storage + d_frootn + d_frootn_storage
    if is_woody || is_crop
        total += d_livestemn + d_livestemn_storage + d_deadstemn + d_deadstemn_storage +
                 d_livecrootn + d_livecrootn_storage + d_deadcrootn + d_deadcrootn_storage
    end
    cng = is_crop ? graincn[ivt] : one(T)
    if is_crop
        for k in 1:nrepr
            f5k = crop_alive ? (arepr[p, k] / aleaf_p[p]) : zero(T)
            total += (nlc * f5k / cng) * fcur + (nlc * f5k / cng) * (one(T) - fcur)
        end
    end

    # 3) draw fractional demand * npool / dt
    inv = total == zero(T) ? zero(T) : (npool / dt) / total
    nf.npool_to_leafn_patch[p]          = d_leafn          * inv
    nf.npool_to_leafn_storage_patch[p]  = d_leafn_storage  * inv
    nf.npool_to_frootn_patch[p]         = d_frootn         * inv
    nf.npool_to_frootn_storage_patch[p] = d_frootn_storage * inv
    if is_woody || is_crop
        nf.npool_to_livestemn_patch[p]          = d_livestemn          * inv
        nf.npool_to_livestemn_storage_patch[p]  = d_livestemn_storage  * inv
        nf.npool_to_deadstemn_patch[p]          = d_deadstemn          * inv
        nf.npool_to_deadstemn_storage_patch[p]  = d_deadstemn_storage  * inv
        nf.npool_to_livecrootn_patch[p]         = d_livecrootn         * inv
        nf.npool_to_livecrootn_storage_patch[p] = d_livecrootn_storage * inv
        nf.npool_to_deadcrootn_patch[p]         = d_deadcrootn         * inv
        nf.npool_to_deadcrootn_storage_patch[p] = d_deadcrootn_storage * inv
    end
    if is_crop
        for k in 1:nrepr
            f5k = crop_alive ? (arepr[p, k] / aleaf_p[p]) : zero(T)
            d_repr         = (nlc * f5k / cng) * fcur
            d_repr_storage = (nlc * f5k / cng) * (one(T) - fcur)
            nf.npool_to_reproductiven_patch[p, k]         = d_repr         * inv
            nf.npool_to_reproductiven_storage_patch[p, k] = d_repr_storage * inv
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# calc_plant_cn_alloc_flexiblecn! -- C/N allocation after competition.
# Ported from calc_plant_cn_alloc (Fortran 193-929).
# ---------------------------------------------------------------------------
function calc_plant_cn_alloc_flexiblecn!(mask_soilp::AbstractVector{Bool},
        bounds::UnitRange{Int},
        pftcon::PftConNutrientCompetition,
        cn_shared_params::CNSharedParamsData,
        patch::PatchData,
        crop::CropData,
        canopystate::CanopyStateData,
        cnveg_state::CNVegStateData,
        cnveg_cs::CNVegCarbonStateData,
        cnveg_cf::CNVegCarbonFluxData,
        cnveg_ns::CNVegNitrogenStateData,
        cnveg_nf::CNVegNitrogenFluxData;
        fpg_col::AbstractVector{<:Real},
        dt::Real=1800.0,
        c13_cnveg_cf::Union{CNVegCarbonFluxData,Nothing}=nothing,
        c14_cnveg_cf::Union{CNVegCarbonFluxData,Nothing}=nothing,
        use_c13::Bool=false,
        use_c14::Bool=false,
        carbon_resp_opt::Int=0,
        use_crop_agsys::Bool=false,
        npcropmin::Int=NPCROPMIN,
        nrepr::Int=NREPR)

    use_fun = cn_shared_params.use_fun
    use_matrixcn = cn_shared_params.use_matrixcn
    cf = cnveg_cf; nf = cnveg_nf; cs = cnveg_cs; ns = cnveg_ns; vs = cnveg_state
    SPVAL = 1.0e36

    @inbounds for p in bounds
        mask_soilp[p] || continue
        T = eltype(cf.plant_calloc_patch)
        c = patch.column[p]
        ivt = patch.itype[p] + 1   # 0-based Fortran → 1-based Julia

        f1 = pftcon.froot_leaf[ivt]
        f2 = pftcon.croot_stem[ivt]

        if pftcon.stem_leaf[ivt] == T(-1.0)
            f3 = (T(2.7) / (one(T) + exp(T(-0.004) * (cf.annsum_npp_patch[p] - T(300.0))))) - T(0.4)
        else
            f3 = pftcon.stem_leaf[ivt]
        end

        f4   = pftcon.flivewd[ivt]
        g1   = pftcon.grperc[ivt]
        g2   = pftcon.grpnow[ivt]
        fcur = pftcon.fcur[ivt]
        if pftcon.evergreen[ivt] == one(T)
            fcur = zero(T)
        end

        is_woody = pftcon.woody[ivt] == one(T)
        is_crop  = ivt >= npcropmin
        crop_alive = is_crop && (crop.croplive_patch[p] != zero(eltype(crop.croplive_patch))) &&
                     !isnan(vs.aleaf_patch[p])

        if is_crop
            if crop_alive
                f1 = vs.aroot_patch[p] / vs.aleaf_patch[p]
                f3 = vs.astem_patch[p] / vs.aleaf_patch[p]
                g1 = T(0.25)
            else
                f1 = zero(T)
                f3 = zero(T)
                g1 = T(0.25)
            end
        end

        # N available: from FUN, or from the FPG calc
        if use_fun
            nf.sminn_to_npool_patch[p] = nf.sminn_to_plant_fun_patch[p]
        else
            nf.sminn_to_npool_patch[p] = nf.plant_ndemand_patch[p] * fpg_col[c]
        end
        nf.plant_nalloc_patch[p] = nf.sminn_to_npool_patch[p] + nf.retransn_to_npool_patch[p]
        if use_matrixcn
            # matrix N input (gN/m2/s) = soil mineral N draw (Fortran 446-452)
            nf.matrix_Ninput_patch[p] = nf.sminn_to_npool_patch[p]
        end

        # C available for new growth: FlexibleCN uses npp_growth (FUN) or availc.
        if use_fun
            cf.plant_calloc_patch[p] = cf.npp_growth_patch[p]
            if use_matrixcn
                cf.matrix_Cinput_patch[p] = cf.npp_growth_patch[p]
            end
        else
            cf.plant_calloc_patch[p] = cf.availc_patch[p]
            if use_matrixcn
                cf.matrix_Cinput_patch[p] = cf.availc_patch[p]
            end
        end

        # new leaf C and C fluxes to current growth + storage
        nlc = cf.plant_calloc_patch[p] / vs.c_allometry_patch[p]
        cf.cpool_to_leafc_patch[p]          = nlc * fcur
        cf.cpool_to_leafc_storage_patch[p]  = nlc * (one(T) - fcur)
        cf.cpool_to_frootc_patch[p]         = nlc * f1 * fcur
        cf.cpool_to_frootc_storage_patch[p] = nlc * f1 * (one(T) - fcur)

        if is_woody
            cf.cpool_to_livestemc_patch[p]          = nlc * f3 * f4 * fcur
            cf.cpool_to_livestemc_storage_patch[p]  = nlc * f3 * f4 * (one(T) - fcur)
            cf.cpool_to_deadstemc_patch[p]          = nlc * f3 * (one(T) - f4) * fcur
            cf.cpool_to_deadstemc_storage_patch[p]  = nlc * f3 * (one(T) - f4) * (one(T) - fcur)
            cf.cpool_to_livecrootc_patch[p]         = nlc * f2 * f3 * f4 * fcur
            cf.cpool_to_livecrootc_storage_patch[p] = nlc * f2 * f3 * f4 * (one(T) - fcur)
            cf.cpool_to_deadcrootc_patch[p]         = nlc * f2 * f3 * (one(T) - f4) * fcur
            cf.cpool_to_deadcrootc_storage_patch[p] = nlc * f2 * f3 * (one(T) - f4) * (one(T) - fcur)
        end
        if is_crop
            cf.cpool_to_livestemc_patch[p]          = nlc * f3 * f4 * fcur
            cf.cpool_to_livestemc_storage_patch[p]  = nlc * f3 * f4 * (one(T) - fcur)
            cf.cpool_to_deadstemc_patch[p]          = nlc * f3 * (one(T) - f4) * fcur
            cf.cpool_to_deadstemc_storage_patch[p]  = nlc * f3 * (one(T) - f4) * (one(T) - fcur)
            cf.cpool_to_livecrootc_patch[p]         = nlc * f2 * f3 * f4 * fcur
            cf.cpool_to_livecrootc_storage_patch[p] = nlc * f2 * f3 * f4 * (one(T) - fcur)
            cf.cpool_to_deadcrootc_patch[p]         = nlc * f2 * f3 * (one(T) - f4) * fcur
            cf.cpool_to_deadcrootc_storage_patch[p] = nlc * f2 * f3 * (one(T) - f4) * (one(T) - fcur)
            for k in 1:nrepr
                f5k = crop_alive ? (vs.arepr_patch[p, k] / vs.aleaf_patch[p]) : zero(T)
                cf.cpool_to_reproductivec_patch[p, k]         = nlc * f5k * fcur
                cf.cpool_to_reproductivec_storage_patch[p, k] = nlc * f5k * (one(T) - fcur)
            end
        end

        # ---- matrixcn C allocation B-matrix (Fortran 479-570) ----
        # cpool_to_veg = total C allocated to all veg pools this step. The
        # matrix_alloc(p, pool) is the fraction of that C input going to each
        # pool (the B-matrix column of the matrix-CN allocation input).
        if use_matrixcn
            cpool_to_veg = cf.cpool_to_leafc_patch[p] + cf.cpool_to_leafc_storage_patch[p] +
                           cf.cpool_to_frootc_patch[p] + cf.cpool_to_frootc_storage_patch[p]
            if is_woody || is_crop
                cpool_to_veg += cf.cpool_to_livestemc_patch[p]  + cf.cpool_to_livestemc_storage_patch[p] +
                                cf.cpool_to_deadstemc_patch[p]  + cf.cpool_to_deadstemc_storage_patch[p] +
                                cf.cpool_to_livecrootc_patch[p] + cf.cpool_to_livecrootc_storage_patch[p] +
                                cf.cpool_to_deadcrootc_patch[p] + cf.cpool_to_deadcrootc_storage_patch[p]
            end
            if is_crop
                for k in 1:nrepr
                    cpool_to_veg += cf.cpool_to_reproductivec_patch[p, k] +
                                    cf.cpool_to_reproductivec_storage_patch[p, k]
                end
            end

            cf.matrix_Cinput_patch[p] = cpool_to_veg
            if cpool_to_veg != zero(T)
                inv = one(T) / cpool_to_veg
                cf.matrix_alloc_patch[p, FLEXCN_ILEAF]     = cf.cpool_to_leafc_patch[p]          * inv
                cf.matrix_alloc_patch[p, FLEXCN_ILEAF_ST]  = cf.cpool_to_leafc_storage_patch[p]  * inv
                cf.matrix_alloc_patch[p, FLEXCN_IFROOT]    = cf.cpool_to_frootc_patch[p]         * inv
                cf.matrix_alloc_patch[p, FLEXCN_IFROOT_ST] = cf.cpool_to_frootc_storage_patch[p] * inv
            end
            if (is_woody || is_crop) && cpool_to_veg != zero(T)
                inv = one(T) / cpool_to_veg
                cf.matrix_alloc_patch[p, FLEXCN_ILIVESTEM]     = cf.cpool_to_livestemc_patch[p]          * inv
                cf.matrix_alloc_patch[p, FLEXCN_ILIVESTEM_ST]  = cf.cpool_to_livestemc_storage_patch[p]  * inv
                cf.matrix_alloc_patch[p, FLEXCN_IDEADSTEM]     = cf.cpool_to_deadstemc_patch[p]          * inv
                cf.matrix_alloc_patch[p, FLEXCN_IDEADSTEM_ST]  = cf.cpool_to_deadstemc_storage_patch[p]  * inv
                cf.matrix_alloc_patch[p, FLEXCN_ILIVECROOT]    = cf.cpool_to_livecrootc_patch[p]         * inv
                cf.matrix_alloc_patch[p, FLEXCN_ILIVECROOT_ST] = cf.cpool_to_livecrootc_storage_patch[p] * inv
                cf.matrix_alloc_patch[p, FLEXCN_IDEADCROOT]    = cf.cpool_to_deadcrootc_patch[p]         * inv
                cf.matrix_alloc_patch[p, FLEXCN_IDEADCROOT_ST] = cf.cpool_to_deadcrootc_storage_patch[p] * inv
                if is_crop
                    cf.matrix_alloc_patch[p, FLEXCN_IGRAIN]    = zero(T)
                    cf.matrix_alloc_patch[p, FLEXCN_IGRAIN_ST] = zero(T)
                    for k in 1:nrepr
                        cf.matrix_alloc_patch[p, FLEXCN_IGRAIN] +=
                            cf.cpool_to_reproductivec_patch[p, k] * inv
                        cf.matrix_alloc_patch[p, FLEXCN_IGRAIN_ST] +=
                            cf.cpool_to_reproductivec_storage_patch[p, k] * inv
                    end
                end
            end
        end

        # growth respiration storage
        gresp_storage = cf.cpool_to_leafc_storage_patch[p] + cf.cpool_to_frootc_storage_patch[p]
        if is_woody
            gresp_storage += cf.cpool_to_livestemc_storage_patch[p]
            gresp_storage += cf.cpool_to_deadstemc_storage_patch[p]
            gresp_storage += cf.cpool_to_livecrootc_storage_patch[p]
            gresp_storage += cf.cpool_to_deadcrootc_storage_patch[p]
        end
        if is_crop
            gresp_storage += cf.cpool_to_livestemc_storage_patch[p]
            for k in 1:nrepr
                gresp_storage += cf.cpool_to_reproductivec_storage_patch[p, k]
            end
        end
        cf.cpool_to_gresp_storage_patch[p] = gresp_storage * g1 * (one(T) - g2)

        # N allocation: AgSys (managed-crop) for crops when use_crop_agsys, else
        # the FlexibleCN fractional-demand draw from the existing npool.
        if use_crop_agsys && is_crop
            # AgSys supplies its own N allocation coefficients (Fortran 597-625)
            _calc_npool_to_components_agsys!(p,
                ns.npool_patch[p], fcur, f4,
                vs.aleaf_n_patch[p], vs.astem_n_patch[p], vs.aroot_n_patch[p],
                vs.arepr_n_patch, nf, dt, nrepr)
        else
            _calc_npool_to_components_flexiblecn!(p, ivt,
                ns.npool_patch[p], nlc, fcur, f1, f2, f3, f4,
                vs.arepr_patch, vs.aleaf_patch, crop_alive, is_crop, is_woody,
                pftcon.leafcn, pftcon.frootcn, pftcon.livewdcn, pftcon.deadwdcn, pftcon.graincn,
                nf, dt, nrepr)
        end

        # ---- carbon_resp_opt: flexible-C:N turnover of high-C:N tissue ----
        cf.cpool_to_resp_patch[p]                  = zero(T)
        cf.cpool_to_leafc_resp_patch[p]            = zero(T)
        cf.cpool_to_leafc_storage_resp_patch[p]    = zero(T)
        cf.cpool_to_frootc_resp_patch[p]           = zero(T)
        cf.cpool_to_frootc_storage_resp_patch[p]   = zero(T)
        cf.cpool_to_livecrootc_resp_patch[p]       = zero(T)
        cf.cpool_to_livecrootc_storage_resp_patch[p] = zero(T)
        cf.cpool_to_livestemc_resp_patch[p]        = zero(T)
        cf.cpool_to_livestemc_storage_resp_patch[p] = zero(T)

        lai_tot = canopystate.laisun_patch[p] + canopystate.laisha_patch[p]

        # storage leaf C:N (history diagnostic, also gates carbon_resp_opt leaf turnover)
        actual_storage_leafcn = T(SPVAL)
        if lai_tot > zero(T)
            if ns.leafn_storage_patch[p] == zero(T)
                actual_storage_leafcn = T(SPVAL)
            else
                actual_storage_leafcn = cs.leafc_storage_patch[p] / ns.leafn_storage_patch[p]
            end
        end

        if carbon_resp_opt == 1 && lai_tot > zero(T)
            nminT = T(NUTRIENT_COMP_N_MIN)
            # fine-root storage C:N
            if ns.frootn_storage_patch[p] == zero(T)
                frootcn_actual = cs.frootc_storage_patch[p] / nminT
            else
                frootcn_actual = cs.frootc_storage_patch[p] / ns.frootn_storage_patch[p]
            end

            livestemcn_actual  = zero(T)
            livecrootcn_actual = zero(T)
            if is_woody || is_crop
                if ns.livestemn_storage_patch[p] == zero(T)
                    livestemcn_actual = cs.livestemc_storage_patch[p] / nminT
                else
                    livestemcn_actual = cs.livestemc_storage_patch[p] / ns.livestemn_storage_patch[p]
                end
                if ns.livecrootn_storage_patch[p] == zero(T)
                    livecrootcn_actual = cs.livecrootc_storage_patch[p] / nminT
                else
                    livecrootcn_actual = cs.livecrootc_storage_patch[p] / ns.livecrootn_storage_patch[p]
                end
            end

            leafcn_max  = pftcon.leafcn[ivt]  + T(15.0)
            frootcn_max = pftcon.frootcn[ivt] + T(15.0)

            # leaf: high storage C:N → turn over extra C
            if actual_storage_leafcn > leafcn_max
                frac_resp = (actual_storage_leafcn - leafcn_max) / T(10.0)
                frac_resp = min(one(T), max(zero(T), frac_resp))
                cf.cpool_to_leafc_resp_patch[p]         = frac_resp * cf.cpool_to_leafc_patch[p]
                cf.cpool_to_leafc_storage_resp_patch[p] = frac_resp * cf.cpool_to_leafc_storage_patch[p]
            end
            # fine root
            if frootcn_actual > frootcn_max
                frac_resp = (frootcn_actual - frootcn_max) / T(10.0)
                frac_resp = min(one(T), max(zero(T), frac_resp))
                cf.cpool_to_frootc_resp_patch[p]         = frac_resp * cf.cpool_to_frootc_patch[p]
                cf.cpool_to_frootc_storage_resp_patch[p] = frac_resp * cf.cpool_to_frootc_storage_patch[p]
            end
            if is_woody || is_crop
                livewdcn_max = pftcon.livewdcn[ivt] + T(15.0)
                if livecrootcn_actual > livewdcn_max
                    frac_resp = (livecrootcn_actual - livewdcn_max) / T(10.0)
                    frac_resp = min(one(T), max(zero(T), frac_resp))
                    cf.cpool_to_livecrootc_resp_patch[p]         = frac_resp * cf.cpool_to_livecrootc_patch[p]
                    cf.cpool_to_livecrootc_storage_resp_patch[p] = frac_resp * cf.cpool_to_livecrootc_storage_patch[p]
                end
                if livestemcn_actual > livewdcn_max
                    frac_resp = (livestemcn_actual - livewdcn_max) / T(10.0)
                    frac_resp = min(one(T), max(zero(T), frac_resp))
                    cf.cpool_to_livestemc_resp_patch[p]         = frac_resp * cf.cpool_to_livestemc_patch[p]
                    cf.cpool_to_livestemc_storage_resp_patch[p] = frac_resp * cf.cpool_to_livestemc_storage_patch[p]
                end
            end

            cf.cpool_to_resp_patch[p] =
                cf.cpool_to_leafc_resp_patch[p] + cf.cpool_to_leafc_storage_resp_patch[p] +
                cf.cpool_to_frootc_resp_patch[p] + cf.cpool_to_frootc_storage_resp_patch[p] +
                cf.cpool_to_livecrootc_resp_patch[p] + cf.cpool_to_livecrootc_storage_resp_patch[p] +
                cf.cpool_to_livestemc_resp_patch[p] + cf.cpool_to_livestemc_storage_resp_patch[p]

            if use_matrixcn
                # the C turned over to respiration is removed from the matrix C
                # input (Fortran 835-837)
                cf.matrix_Cinput_patch[p] -= cf.cpool_to_resp_patch[p]
            end
        end

        # ---- matrixcn N allocation B-matrix + isotope C inputs (Fortran 847-924) ----
        # FlexibleCN has NO C13/C14 downregulation, so under the sequential path
        # the isotope psn fluxes are unchanged. Under use_matrixcn the matrix
        # isotope C inputs and the N allocation B-matrix (matrix_nalloc) are
        # written here from the npool_to_* fluxes just computed.
        if use_matrixcn
            psn = cf.psnsun_to_cpool_patch[p] + cf.psnshade_to_cpool_patch[p]
            if use_c13 && c13_cnveg_cf !== nothing && psn != zero(T)
                cf.matrix_C13input_patch[p] = cf.plant_calloc_patch[p] *
                    ((c13_cnveg_cf.psnsun_to_cpool_patch[p] +
                      c13_cnveg_cf.psnshade_to_cpool_patch[p]) / psn)
            end
            if use_c14 && c14_cnveg_cf !== nothing && psn != zero(T)
                cf.matrix_C14input_patch[p] = cf.plant_calloc_patch[p] *
                    ((c14_cnveg_cf.psnsun_to_cpool_patch[p] +
                      c14_cnveg_cf.psnshade_to_cpool_patch[p]) / psn)
            end

            npool_to_veg = nf.npool_to_leafn_patch[p]      + nf.npool_to_leafn_storage_patch[p] +
                           nf.npool_to_frootn_patch[p]     + nf.npool_to_frootn_storage_patch[p] +
                           nf.npool_to_livestemn_patch[p]  + nf.npool_to_livestemn_storage_patch[p] +
                           nf.npool_to_deadstemn_patch[p]  + nf.npool_to_deadstemn_storage_patch[p] +
                           nf.npool_to_livecrootn_patch[p] + nf.npool_to_livecrootn_storage_patch[p] +
                           nf.npool_to_deadcrootn_patch[p] + nf.npool_to_deadcrootn_storage_patch[p]
            if is_crop
                npool_to_veg += nf.npool_to_reproductiven_patch[p, 1] +
                                nf.npool_to_reproductiven_storage_patch[p, 1]
            end

            if npool_to_veg != zero(T)
                invn = one(T) / npool_to_veg
                nf.matrix_nalloc_patch[p, FLEXCN_ILEAF]        = nf.npool_to_leafn_patch[p]          * invn
                nf.matrix_nalloc_patch[p, FLEXCN_ILEAF_ST]     = nf.npool_to_leafn_storage_patch[p]  * invn
                nf.matrix_nalloc_patch[p, FLEXCN_IFROOT]       = nf.npool_to_frootn_patch[p]         * invn
                nf.matrix_nalloc_patch[p, FLEXCN_IFROOT_ST]    = nf.npool_to_frootn_storage_patch[p] * invn
                nf.matrix_nalloc_patch[p, FLEXCN_ILIVESTEM]    = nf.npool_to_livestemn_patch[p]          * invn
                nf.matrix_nalloc_patch[p, FLEXCN_ILIVESTEM_ST] = nf.npool_to_livestemn_storage_patch[p]  * invn
                nf.matrix_nalloc_patch[p, FLEXCN_IDEADSTEM]    = nf.npool_to_deadstemn_patch[p]          * invn
                nf.matrix_nalloc_patch[p, FLEXCN_IDEADSTEM_ST] = nf.npool_to_deadstemn_storage_patch[p]  * invn
                nf.matrix_nalloc_patch[p, FLEXCN_ILIVECROOT]    = nf.npool_to_livecrootn_patch[p]         * invn
                nf.matrix_nalloc_patch[p, FLEXCN_ILIVECROOT_ST] = nf.npool_to_livecrootn_storage_patch[p] * invn
                nf.matrix_nalloc_patch[p, FLEXCN_IDEADCROOT]    = nf.npool_to_deadcrootn_patch[p]         * invn
                nf.matrix_nalloc_patch[p, FLEXCN_IDEADCROOT_ST] = nf.npool_to_deadcrootn_storage_patch[p] * invn
                if is_crop
                    nf.matrix_nalloc_patch[p, FLEXCN_IGRAIN]    = nf.npool_to_reproductiven_patch[p, 1]         * invn
                    nf.matrix_nalloc_patch[p, FLEXCN_IGRAIN_ST] = nf.npool_to_reproductiven_storage_patch[p, 1] * invn
                end
                nf.matrix_Ninput_patch[p] = npool_to_veg - nf.retransn_to_npool_patch[p]
            end

            # DEFERRED: the retranslocation-pool transfer registrations
            # (Fortran 900-922, matrix_update_phn into iretransn_to_* transfer
            # indices) require the N-matrix turnover/transfer machinery
            # (matrix_update_phn! + populated iretransn_to_*_ph indices), which
            # is NOT yet ported (see cn_veg_matrix.jl: the explicit N-matrix
            # solve is deferred). The matrix_nalloc / matrix_Ninput B-matrix
            # writes above ARE ported; the retransn transfer accounting is the
            # only piece that depends on that missing infrastructure.
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# calc_plant_nutrient_demand_flexiblecn! -- public N-demand entry.
# Wraps calc_plant_nitrogen_demand_flexiblecn!.
# ---------------------------------------------------------------------------

"""
    calc_plant_nutrient_demand_flexiblecn!(mask_p, bounds, call_is_for_pcrop,
        pftcon, cn_shared_params, patch, crop, canopystate,
        cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf,
        soilbgc_ns, soilbgc_cf;
        dt, dzsoi_decomp, nlevdecomp, npcropmin, nrepr,
        ntmp_soybean, nirrig_tmp_soybean, ntrp_soybean, nirrig_trp_soybean)

FlexibleCN plant nitrogen demand. Parallel of `calc_plant_nutrient_demand!`.
Plant N demand is computed from Michaelis-Menten kinetics on fine-root C, a
soil mineral-N substrate term, a soil-temperature scalar and a flexible
leaf-C:N nitrogen scalar.

Ported from `calc_plant_nitrogen_demand` in
`NutrientCompetitionFlexibleCNMod.F90`.
"""
function calc_plant_nutrient_demand_flexiblecn!(mask_p::AbstractVector{Bool},
        bounds::UnitRange{Int}, call_is_for_pcrop::Bool,
        pftcon::PftConNutrientCompetition,
        cn_shared_params::CNSharedParamsData,
        patch::PatchData,
        crop::CropData,
        canopystate::CanopyStateData,
        cnveg_state::CNVegStateData,
        cnveg_cs::CNVegCarbonStateData,
        cnveg_cf::CNVegCarbonFluxData,
        cnveg_ns::CNVegNitrogenStateData,
        cnveg_nf::CNVegNitrogenFluxData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData,
        soilbgc_cf::SoilBiogeochemCarbonFluxData;
        dt::Real,
        dzsoi_decomp::AbstractVector{<:Real},
        nlevdecomp::Int,
        npcropmin::Int=NPCROPMIN,
        nrepr::Int=NREPR,
        ntmp_soybean::Int=23,
        nirrig_tmp_soybean::Int=24,
        ntrp_soybean::Int=77,
        nirrig_trp_soybean::Int=78)

    use_fun = cn_shared_params.use_fun
    cf = cnveg_cf; nf = cnveg_nf; cs = cnveg_cs; ns = cnveg_ns; vs = cnveg_state
    SPVAL = 1.0e36
    Vmax_N = 2.7e-8
    Kmin   = 1.0
    T = eltype(nf.plant_ndemand_patch)

    @inbounds for p in bounds
        mask_p[p] || continue
        c = patch.column[p]
        ivt = patch.itype[p] + 1

        # actual leaf C:N, flexed between leafcn±10
        if ns.leafn_patch[p] < T(NUTRIENT_COMP_N_MIN)
            actual_leafcn = T(SPVAL)
        else
            actual_leafcn = cs.leafc_patch[p] / ns.leafn_patch[p]
        end
        leafcn_min = pftcon.leafcn[ivt] - T(10.0)
        leafcn_max = pftcon.leafcn[ivt] + T(10.0)
        actual_leafcn = max(actual_leafcn, leafcn_min - T(0.0001))
        actual_leafcn = min(actual_leafcn, leafcn_max)

        nscalar = (actual_leafcn - leafcn_min) / (leafcn_max - leafcn_min)
        nscalar = min(max(zero(T), nscalar), one(T))

        # soil mineral-N substrate term (column, vertically integrated)
        sminn_total = zero(T)
        for j in 1:nlevdecomp
            sminn_total += soilbgc_ns.sminn_vr_col[c, j] * dzsoi_decomp[j]
        end
        substrate_term = sminn_total / (sminn_total + T(Kmin))

        temp_scalar = soilbgc_cf.t_scalar_col[c, 1]
        temp_scalar = min(max(zero(T), temp_scalar), one(T))

        if use_fun
            # FUN: plant_ndemand is just a max draw on soil N pools
            nf.plant_ndemand_patch[p] = cf.availc_patch[p] *
                (vs.n_allometry_patch[p] / vs.c_allometry_patch[p])
        else
            lai_tot = canopystate.laisun_patch[p] + canopystate.laisha_patch[p]
            if lai_tot > zero(T)
                nf.plant_ndemand_patch[p] = T(Vmax_N) * cs.frootc_patch[p] *
                    substrate_term * temp_scalar * nscalar
            else
                nf.plant_ndemand_patch[p] = zero(T)
            end
            if actual_leafcn < leafcn_min
                nf.plant_ndemand_patch[p] = zero(T)
            end
        end

        # retranslocated N deployment depends on seasonal cycle of potential GPP
        vs.tempsum_potential_gpp_patch[p] += cf.gpp_before_downreg_patch[p]
        vs.tempmax_retransn_patch[p] = max(vs.tempmax_retransn_patch[p], ns.retransn_patch[p])
    end

    # crop grain-fill retranslocation (pool-relative form)
    if call_is_for_pcrop
        crop_phase_vals = fill!(similar(nf.plant_ndemand_patch, T, length(mask_p)), zero(T))
        crop_phase!(mask_p, crop, vs, crop_phase_vals)

        @inbounds for p in bounds
            mask_p[p] || continue
            (crop.croplive_patch[p] != zero(eltype(crop.croplive_patch))) || continue
            ivt = patch.itype[p] + 1
            cphase = crop_phase_vals[p]
            if cphase == T(cphase_leafemerge)
                vs.grain_flag_patch[p] = zero(T)
            elseif cphase == T(cphase_grainfill)
                is_soybean = (ivt == ntmp_soybean || ivt == nirrig_tmp_soybean ||
                              ivt == ntrp_soybean || ivt == nirrig_trp_soybean)
                if vs.astem_patch[p] == pftcon.astemf[ivt] || !is_soybean
                    if vs.grain_flag_patch[p] == zero(T)
                        t1 = one(T) / T(dt)
                        nf.leafn_to_retransn_patch[p] = t1 *
                            max(ns.leafn_patch[p] - (cs.leafc_patch[p] / pftcon.fleafcn[ivt]), zero(T))
                        nf.livestemn_to_retransn_patch[p] = t1 *
                            max(ns.livestemn_patch[p] - (cs.livestemc_patch[p] / pftcon.fstemcn[ivt]), zero(T))
                        nf.frootn_to_retransn_patch[p] = zero(T)
                        if pftcon.ffrootcn[ivt] > zero(T)
                            nf.frootn_to_retransn_patch[p] = t1 *
                                max(ns.frootn_patch[p] - (cs.frootc_patch[p] / pftcon.ffrootcn[ivt]), zero(T))
                        end
                        vs.grain_flag_patch[p] = one(T)
                    end
                end
            end
        end
    end

    # crops pull from retransn only during grain fill; others get seasonal fraction
    if call_is_for_pcrop
        @inbounds for p in bounds
            mask_p[p] || continue
            if vs.grain_flag_patch[p] == one(T)
                nf.avail_retransn_patch[p] = nf.plant_ndemand_patch[p]
            else
                nf.avail_retransn_patch[p] = zero(T)
            end
        end
    else
        @inbounds for p in bounds
            mask_p[p] || continue
            if vs.annsum_potential_gpp_patch[p] > zero(T)
                nf.avail_retransn_patch[p] = (vs.annmax_retransn_patch[p] / T(2.0)) *
                    (cf.gpp_before_downreg_patch[p] / vs.annsum_potential_gpp_patch[p]) / T(dt)
            else
                nf.avail_retransn_patch[p] = zero(T)
            end
        end
    end

    # clamp to storage + modify plant N demand
    @inbounds for p in bounds
        mask_p[p] || continue
        nf.avail_retransn_patch[p] = min(nf.avail_retransn_patch[p], ns.retransn_patch[p] / T(dt))
        if nf.plant_ndemand_patch[p] > nf.avail_retransn_patch[p]
            nf.retransn_to_npool_patch[p] = nf.avail_retransn_patch[p]
        else
            nf.retransn_to_npool_patch[p] = nf.plant_ndemand_patch[p]
        end
        if !use_fun
            nf.plant_ndemand_patch[p] -= nf.retransn_to_npool_patch[p]
        end
    end

    return nothing
end
