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
#
# GPU kernelization (mirrors nutrient_competition.jl): both public entry points
# run as KernelAbstractions kernels (one thread per patch — every write is to the
# patch's own index, so it is race-free and byte-identical to the sequential host
# loop on the KA CPU backend). The many touched carbon/nitrogen state & flux
# arrays are grouped into `Adapt.@adapt_structure` device-view bundles to stay
# under Metal's ~31-arg launch limit; pftcon params (a host-global) travel inside
# the input bundle. The two per-patch npool-allocation helpers are inlined into the
# allocation kernel. All literals are carried at the working element type `T`, so
# the kernel compiles on Metal (no Float64) while the Float64 host path is unchanged.
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

# ===========================================================================
# Allocation: device-view bundles + main kernel
# ===========================================================================

# Written outputs (carbon-flux + nitrogen-flux patch arrays). V = float vectors,
# M = float matrices (the [patch, k]/[patch, pool] arrays).
Base.@kwdef struct _FlexAllocOut{V,M}
    sminn_to_npool::V; plant_nalloc::V; plant_calloc::V
    cpool_to_leafc::V; cpool_to_leafc_storage::V
    cpool_to_frootc::V; cpool_to_frootc_storage::V
    cpool_to_livestemc::V; cpool_to_livestemc_storage::V
    cpool_to_deadstemc::V; cpool_to_deadstemc_storage::V
    cpool_to_livecrootc::V; cpool_to_livecrootc_storage::V
    cpool_to_deadcrootc::V; cpool_to_deadcrootc_storage::V
    cpool_to_reproductivec::M; cpool_to_reproductivec_storage::M
    cpool_to_gresp_storage::V
    npool_to_leafn::V; npool_to_leafn_storage::V
    npool_to_frootn::V; npool_to_frootn_storage::V
    npool_to_livestemn::V; npool_to_livestemn_storage::V
    npool_to_deadstemn::V; npool_to_deadstemn_storage::V
    npool_to_livecrootn::V; npool_to_livecrootn_storage::V
    npool_to_deadcrootn::V; npool_to_deadcrootn_storage::V
    npool_to_reproductiven::M; npool_to_reproductiven_storage::M
    cpool_to_resp::V
    cpool_to_leafc_resp::V; cpool_to_leafc_storage_resp::V
    cpool_to_frootc_resp::V; cpool_to_frootc_storage_resp::V
    cpool_to_livecrootc_resp::V; cpool_to_livecrootc_storage_resp::V
    cpool_to_livestemc_resp::V; cpool_to_livestemc_storage_resp::V
    matrix_Cinput::V; matrix_Ninput::V; matrix_C13input::V; matrix_C14input::V
    matrix_alloc::M; matrix_nalloc::M
end
Adapt.@adapt_structure _FlexAllocOut

# Read-only inputs (state reads + pftcon params). VB = the (Bool-ish) croplive vec.
Base.@kwdef struct _FlexAllocIn{V,M,VB}
    # pftcon parameter vectors (indexed by ivt)
    froot_leaf::V; croot_stem::V; stem_leaf::V; flivewd::V; grperc::V; grpnow::V
    fcur_arr::V; evergreen::V; woody::V
    leafcn::V; frootcn::V; livewdcn::V; deadwdcn::V; graincn::V
    # cnveg carbon-flux reads
    annsum_npp::V; availc::V; npp_growth::V; psnsun_to_cpool::V; psnshade_to_cpool::V
    # cnveg nitrogen-flux reads
    sminn_to_plant_fun::V; plant_ndemand::V; retransn_to_npool::V
    fpg_col::V
    # cnveg carbon-state storage reads
    leafc_storage::V; frootc_storage::V; livestemc_storage::V; livecrootc_storage::V
    # cnveg nitrogen-state reads
    npool::V; leafn_storage::V; frootn_storage::V; livestemn_storage::V; livecrootn_storage::V
    # cnveg (vegetation) state reads
    c_allometry::V; aleaf::V; aroot::V; astem::V; arepr::M
    aleaf_n::V; astem_n::V; aroot_n::V; arepr_n::M
    # canopy state reads
    laisun::V; laisha::V
    # crop state read
    croplive::VB
end
Adapt.@adapt_structure _FlexAllocIn

# --- Allocation kernel: one thread per patch --------------------------------
@kernel function _flexalloc_main_kernel!(out::_FlexAllocOut, in::_FlexAllocIn,
        @Const(mask_soilp), @Const(column), @Const(itype),
        @Const(c13_psnsun), @Const(c13_psnshade), @Const(c14_psnsun), @Const(c14_psnshade),
        use_fun::Bool, use_matrixcn::Bool, use_crop_agsys::Bool, carbon_resp_opt::Int,
        use_c13::Bool, use_c14::Bool, npcropmin::Int, nrepr::Int, dt)
    T = eltype(out.plant_calloc)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        # alias output struct fields to their loose (Fortran) names
        sminn_to_npool = out.sminn_to_npool; plant_nalloc = out.plant_nalloc
        plant_calloc = out.plant_calloc
        cpool_to_leafc = out.cpool_to_leafc; cpool_to_leafc_storage = out.cpool_to_leafc_storage
        cpool_to_frootc = out.cpool_to_frootc; cpool_to_frootc_storage = out.cpool_to_frootc_storage
        cpool_to_livestemc = out.cpool_to_livestemc; cpool_to_livestemc_storage = out.cpool_to_livestemc_storage
        cpool_to_deadstemc = out.cpool_to_deadstemc; cpool_to_deadstemc_storage = out.cpool_to_deadstemc_storage
        cpool_to_livecrootc = out.cpool_to_livecrootc; cpool_to_livecrootc_storage = out.cpool_to_livecrootc_storage
        cpool_to_deadcrootc = out.cpool_to_deadcrootc; cpool_to_deadcrootc_storage = out.cpool_to_deadcrootc_storage
        cpool_to_reproductivec = out.cpool_to_reproductivec; cpool_to_reproductivec_storage = out.cpool_to_reproductivec_storage
        cpool_to_gresp_storage = out.cpool_to_gresp_storage
        npool_to_leafn = out.npool_to_leafn; npool_to_leafn_storage = out.npool_to_leafn_storage
        npool_to_frootn = out.npool_to_frootn; npool_to_frootn_storage = out.npool_to_frootn_storage
        npool_to_livestemn = out.npool_to_livestemn; npool_to_livestemn_storage = out.npool_to_livestemn_storage
        npool_to_deadstemn = out.npool_to_deadstemn; npool_to_deadstemn_storage = out.npool_to_deadstemn_storage
        npool_to_livecrootn = out.npool_to_livecrootn; npool_to_livecrootn_storage = out.npool_to_livecrootn_storage
        npool_to_deadcrootn = out.npool_to_deadcrootn; npool_to_deadcrootn_storage = out.npool_to_deadcrootn_storage
        npool_to_reproductiven = out.npool_to_reproductiven; npool_to_reproductiven_storage = out.npool_to_reproductiven_storage
        cpool_to_resp = out.cpool_to_resp
        cpool_to_leafc_resp = out.cpool_to_leafc_resp; cpool_to_leafc_storage_resp = out.cpool_to_leafc_storage_resp
        cpool_to_frootc_resp = out.cpool_to_frootc_resp; cpool_to_frootc_storage_resp = out.cpool_to_frootc_storage_resp
        cpool_to_livecrootc_resp = out.cpool_to_livecrootc_resp; cpool_to_livecrootc_storage_resp = out.cpool_to_livecrootc_storage_resp
        cpool_to_livestemc_resp = out.cpool_to_livestemc_resp; cpool_to_livestemc_storage_resp = out.cpool_to_livestemc_storage_resp
        matrix_Cinput = out.matrix_Cinput; matrix_Ninput = out.matrix_Ninput
        matrix_C13input = out.matrix_C13input; matrix_C14input = out.matrix_C14input
        matrix_alloc = out.matrix_alloc; matrix_nalloc = out.matrix_nalloc
        # alias input struct fields
        froot_leaf = in.froot_leaf; croot_stem = in.croot_stem; stem_leaf = in.stem_leaf
        flivewd = in.flivewd; grperc = in.grperc; grpnow = in.grpnow
        fcur_arr = in.fcur_arr; evergreen = in.evergreen; woody = in.woody
        leafcn = in.leafcn; frootcn = in.frootcn; livewdcn = in.livewdcn
        deadwdcn = in.deadwdcn; graincn = in.graincn
        annsum_npp = in.annsum_npp; availc = in.availc; npp_growth = in.npp_growth
        psnsun_to_cpool = in.psnsun_to_cpool; psnshade_to_cpool = in.psnshade_to_cpool
        sminn_to_plant_fun = in.sminn_to_plant_fun; plant_ndemand = in.plant_ndemand
        retransn_to_npool = in.retransn_to_npool
        leafc_storage = in.leafc_storage; frootc_storage = in.frootc_storage
        livestemc_storage = in.livestemc_storage; livecrootc_storage = in.livecrootc_storage
        npool = in.npool; leafn_storage = in.leafn_storage; frootn_storage = in.frootn_storage
        livestemn_storage = in.livestemn_storage; livecrootn_storage = in.livecrootn_storage
        c_allometry = in.c_allometry; aleaf = in.aleaf; aroot = in.aroot; astem = in.astem
        arepr = in.arepr; aleaf_n = in.aleaf_n; astem_n = in.astem_n; aroot_n = in.aroot_n
        arepr_n = in.arepr_n; laisun = in.laisun; laisha = in.laisha; croplive = in.croplive

        SPVAL = T(1.0e36)
        c = column[p]
        ivt = itype[p] + 1   # 0-based Fortran → 1-based Julia

        f1 = froot_leaf[ivt]
        f2 = croot_stem[ivt]

        if stem_leaf[ivt] == T(-1.0)
            f3 = (T(2.7) / (one(T) + exp(T(-0.004) * (annsum_npp[p] - T(300.0))))) - T(0.4)
        else
            f3 = stem_leaf[ivt]
        end

        f4   = flivewd[ivt]
        g1   = grperc[ivt]
        g2   = grpnow[ivt]
        fcur = fcur_arr[ivt]
        if evergreen[ivt] == one(T)
            fcur = zero(T)
        end

        is_woody = woody[ivt] == one(T)
        is_crop  = ivt >= npcropmin
        crop_alive = is_crop && (croplive[p] != zero(eltype(croplive))) &&
                     !isnan(aleaf[p])

        if is_crop
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

        # N available: from FUN, or from the FPG calc
        if use_fun
            sminn_to_npool[p] = sminn_to_plant_fun[p]
        else
            sminn_to_npool[p] = plant_ndemand[p] * in.fpg_col[c]
        end
        plant_nalloc[p] = sminn_to_npool[p] + retransn_to_npool[p]
        if use_matrixcn
            matrix_Ninput[p] = sminn_to_npool[p]
        end

        # C available for new growth: FlexibleCN uses npp_growth (FUN) or availc.
        if use_fun
            plant_calloc[p] = npp_growth[p]
            if use_matrixcn
                matrix_Cinput[p] = npp_growth[p]
            end
        else
            plant_calloc[p] = availc[p]
            if use_matrixcn
                matrix_Cinput[p] = availc[p]
            end
        end

        # new leaf C and C fluxes to current growth + storage
        nlc = plant_calloc[p] / c_allometry[p]
        cpool_to_leafc[p]          = nlc * fcur
        cpool_to_leafc_storage[p]  = nlc * (one(T) - fcur)
        cpool_to_frootc[p]         = nlc * f1 * fcur
        cpool_to_frootc_storage[p] = nlc * f1 * (one(T) - fcur)

        if is_woody
            cpool_to_livestemc[p]          = nlc * f3 * f4 * fcur
            cpool_to_livestemc_storage[p]  = nlc * f3 * f4 * (one(T) - fcur)
            cpool_to_deadstemc[p]          = nlc * f3 * (one(T) - f4) * fcur
            cpool_to_deadstemc_storage[p]  = nlc * f3 * (one(T) - f4) * (one(T) - fcur)
            cpool_to_livecrootc[p]         = nlc * f2 * f3 * f4 * fcur
            cpool_to_livecrootc_storage[p] = nlc * f2 * f3 * f4 * (one(T) - fcur)
            cpool_to_deadcrootc[p]         = nlc * f2 * f3 * (one(T) - f4) * fcur
            cpool_to_deadcrootc_storage[p] = nlc * f2 * f3 * (one(T) - f4) * (one(T) - fcur)
        end
        if is_crop
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

        # ---- matrixcn C allocation B-matrix ----
        if use_matrixcn
            cpool_to_veg = cpool_to_leafc[p] + cpool_to_leafc_storage[p] +
                           cpool_to_frootc[p] + cpool_to_frootc_storage[p]
            if is_woody || is_crop
                cpool_to_veg += cpool_to_livestemc[p]  + cpool_to_livestemc_storage[p] +
                                cpool_to_deadstemc[p]  + cpool_to_deadstemc_storage[p] +
                                cpool_to_livecrootc[p] + cpool_to_livecrootc_storage[p] +
                                cpool_to_deadcrootc[p] + cpool_to_deadcrootc_storage[p]
            end
            if is_crop
                for k in 1:nrepr
                    cpool_to_veg += cpool_to_reproductivec[p, k] +
                                    cpool_to_reproductivec_storage[p, k]
                end
            end

            matrix_Cinput[p] = cpool_to_veg
            if cpool_to_veg != zero(T)
                inv = one(T) / cpool_to_veg
                matrix_alloc[p, FLEXCN_ILEAF]     = cpool_to_leafc[p]          * inv
                matrix_alloc[p, FLEXCN_ILEAF_ST]  = cpool_to_leafc_storage[p]  * inv
                matrix_alloc[p, FLEXCN_IFROOT]    = cpool_to_frootc[p]         * inv
                matrix_alloc[p, FLEXCN_IFROOT_ST] = cpool_to_frootc_storage[p] * inv
            end
            if (is_woody || is_crop) && cpool_to_veg != zero(T)
                inv = one(T) / cpool_to_veg
                matrix_alloc[p, FLEXCN_ILIVESTEM]     = cpool_to_livestemc[p]          * inv
                matrix_alloc[p, FLEXCN_ILIVESTEM_ST]  = cpool_to_livestemc_storage[p]  * inv
                matrix_alloc[p, FLEXCN_IDEADSTEM]     = cpool_to_deadstemc[p]          * inv
                matrix_alloc[p, FLEXCN_IDEADSTEM_ST]  = cpool_to_deadstemc_storage[p]  * inv
                matrix_alloc[p, FLEXCN_ILIVECROOT]    = cpool_to_livecrootc[p]         * inv
                matrix_alloc[p, FLEXCN_ILIVECROOT_ST] = cpool_to_livecrootc_storage[p] * inv
                matrix_alloc[p, FLEXCN_IDEADCROOT]    = cpool_to_deadcrootc[p]         * inv
                matrix_alloc[p, FLEXCN_IDEADCROOT_ST] = cpool_to_deadcrootc_storage[p] * inv
                if is_crop
                    matrix_alloc[p, FLEXCN_IGRAIN]    = zero(T)
                    matrix_alloc[p, FLEXCN_IGRAIN_ST] = zero(T)
                    ag = zero(T); ags = zero(T)
                    for k in 1:nrepr
                        ag  += cpool_to_reproductivec[p, k]         * inv
                        ags += cpool_to_reproductivec_storage[p, k] * inv
                    end
                    matrix_alloc[p, FLEXCN_IGRAIN]    = ag
                    matrix_alloc[p, FLEXCN_IGRAIN_ST] = ags
                end
            end
        end

        # growth respiration storage
        gresp_storage = cpool_to_leafc_storage[p] + cpool_to_frootc_storage[p]
        if is_woody
            gresp_storage += cpool_to_livestemc_storage[p]
            gresp_storage += cpool_to_deadstemc_storage[p]
            gresp_storage += cpool_to_livecrootc_storage[p]
            gresp_storage += cpool_to_deadcrootc_storage[p]
        end
        if is_crop
            gresp_storage += cpool_to_livestemc_storage[p]
            for k in 1:nrepr
                gresp_storage += cpool_to_reproductivec_storage[p, k]
            end
        end
        cpool_to_gresp_storage[p] = gresp_storage * g1 * (one(T) - g2)

        # ---- N allocation from npool: AgSys (managed crop) or FlexibleCN draw ----
        if use_crop_agsys && is_crop
            # AgSys: distribute npool by supplied coefficients (aleaf_n/astem_n/aroot_n/arepr_n)
            npool_dt = npool[p] / dt
            npool_to_leafn[p]          = aleaf_n[p] * fcur           * npool_dt
            npool_to_leafn_storage[p]  = aleaf_n[p] * (one(T) - fcur) * npool_dt
            npool_to_frootn[p]         = aroot_n[p] * fcur           * npool_dt
            npool_to_frootn_storage[p] = aroot_n[p] * (one(T) - fcur) * npool_dt
            npool_to_livestemn[p]          = astem_n[p] * f4           * fcur           * npool_dt
            npool_to_livestemn_storage[p]  = astem_n[p] * f4           * (one(T) - fcur) * npool_dt
            npool_to_deadstemn[p]          = astem_n[p] * (one(T) - f4) * fcur           * npool_dt
            npool_to_deadstemn_storage[p]  = astem_n[p] * (one(T) - f4) * (one(T) - fcur) * npool_dt
            npool_to_livecrootn[p]         = zero(T)
            npool_to_livecrootn_storage[p] = zero(T)
            npool_to_deadcrootn[p]         = zero(T)
            npool_to_deadcrootn_storage[p] = zero(T)
            for k in 1:nrepr
                npool_to_reproductiven[p, k]         = arepr_n[p, k] * fcur           * npool_dt
                npool_to_reproductiven_storage[p, k] = arepr_n[p, k] * (one(T) - fcur) * npool_dt
            end
        else
            # FlexibleCN fractional-demand draw from the existing npool.
            cnl  = leafcn[ivt]; cnfr = frootcn[ivt]; cnlw = livewdcn[ivt]; cndw = deadwdcn[ivt]
            d_leafn          = (nlc / cnl) * fcur
            d_leafn_storage  = (nlc / cnl) * (one(T) - fcur)
            d_frootn         = (nlc * f1 / cnfr) * fcur
            d_frootn_storage = (nlc * f1 / cnfr) * (one(T) - fcur)
            d_livestemn = zero(T); d_livestemn_storage = zero(T)
            d_deadstemn = zero(T); d_deadstemn_storage = zero(T)
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
            total = d_leafn + d_leafn_storage + d_frootn + d_frootn_storage
            if is_woody || is_crop
                total += d_livestemn + d_livestemn_storage + d_deadstemn + d_deadstemn_storage +
                         d_livecrootn + d_livecrootn_storage + d_deadcrootn + d_deadcrootn_storage
            end
            cng = is_crop ? graincn[ivt] : one(T)
            if is_crop
                for k in 1:nrepr
                    f5k = crop_alive ? (arepr[p, k] / aleaf[p]) : zero(T)
                    total += (nlc * f5k / cng) * fcur + (nlc * f5k / cng) * (one(T) - fcur)
                end
            end
            invn = total == zero(T) ? zero(T) : (npool[p] / dt) / total
            npool_to_leafn[p]          = d_leafn          * invn
            npool_to_leafn_storage[p]  = d_leafn_storage  * invn
            npool_to_frootn[p]         = d_frootn         * invn
            npool_to_frootn_storage[p] = d_frootn_storage * invn
            if is_woody || is_crop
                npool_to_livestemn[p]          = d_livestemn          * invn
                npool_to_livestemn_storage[p]  = d_livestemn_storage  * invn
                npool_to_deadstemn[p]          = d_deadstemn          * invn
                npool_to_deadstemn_storage[p]  = d_deadstemn_storage  * invn
                npool_to_livecrootn[p]         = d_livecrootn         * invn
                npool_to_livecrootn_storage[p] = d_livecrootn_storage * invn
                npool_to_deadcrootn[p]         = d_deadcrootn         * invn
                npool_to_deadcrootn_storage[p] = d_deadcrootn_storage * invn
            end
            if is_crop
                for k in 1:nrepr
                    f5k = crop_alive ? (arepr[p, k] / aleaf[p]) : zero(T)
                    d_repr         = (nlc * f5k / cng) * fcur
                    d_repr_storage = (nlc * f5k / cng) * (one(T) - fcur)
                    npool_to_reproductiven[p, k]         = d_repr         * invn
                    npool_to_reproductiven_storage[p, k] = d_repr_storage * invn
                end
            end
        end

        # ---- carbon_resp_opt: flexible-C:N turnover of high-C:N tissue ----
        cpool_to_resp[p]                  = zero(T)
        cpool_to_leafc_resp[p]            = zero(T)
        cpool_to_leafc_storage_resp[p]    = zero(T)
        cpool_to_frootc_resp[p]           = zero(T)
        cpool_to_frootc_storage_resp[p]   = zero(T)
        cpool_to_livecrootc_resp[p]       = zero(T)
        cpool_to_livecrootc_storage_resp[p] = zero(T)
        cpool_to_livestemc_resp[p]        = zero(T)
        cpool_to_livestemc_storage_resp[p] = zero(T)

        lai_tot = laisun[p] + laisha[p]

        # storage leaf C:N (history diagnostic, also gates carbon_resp_opt leaf turnover)
        actual_storage_leafcn = SPVAL
        if lai_tot > zero(T)
            if leafn_storage[p] == zero(T)
                actual_storage_leafcn = SPVAL
            else
                actual_storage_leafcn = leafc_storage[p] / leafn_storage[p]
            end
        end

        if carbon_resp_opt == 1 && lai_tot > zero(T)
            nminT = T(NUTRIENT_COMP_N_MIN)
            if frootn_storage[p] == zero(T)
                frootcn_actual = frootc_storage[p] / nminT
            else
                frootcn_actual = frootc_storage[p] / frootn_storage[p]
            end

            livestemcn_actual  = zero(T)
            livecrootcn_actual = zero(T)
            if is_woody || is_crop
                if livestemn_storage[p] == zero(T)
                    livestemcn_actual = livestemc_storage[p] / nminT
                else
                    livestemcn_actual = livestemc_storage[p] / livestemn_storage[p]
                end
                if livecrootn_storage[p] == zero(T)
                    livecrootcn_actual = livecrootc_storage[p] / nminT
                else
                    livecrootcn_actual = livecrootc_storage[p] / livecrootn_storage[p]
                end
            end

            leafcn_max  = leafcn[ivt]  + T(15.0)
            frootcn_max = frootcn[ivt] + T(15.0)

            if actual_storage_leafcn > leafcn_max
                frac_resp = (actual_storage_leafcn - leafcn_max) / T(10.0)
                frac_resp = min(one(T), max(zero(T), frac_resp))
                cpool_to_leafc_resp[p]         = frac_resp * cpool_to_leafc[p]
                cpool_to_leafc_storage_resp[p] = frac_resp * cpool_to_leafc_storage[p]
            end
            if frootcn_actual > frootcn_max
                frac_resp = (frootcn_actual - frootcn_max) / T(10.0)
                frac_resp = min(one(T), max(zero(T), frac_resp))
                cpool_to_frootc_resp[p]         = frac_resp * cpool_to_frootc[p]
                cpool_to_frootc_storage_resp[p] = frac_resp * cpool_to_frootc_storage[p]
            end
            if is_woody || is_crop
                livewdcn_max = livewdcn[ivt] + T(15.0)
                if livecrootcn_actual > livewdcn_max
                    frac_resp = (livecrootcn_actual - livewdcn_max) / T(10.0)
                    frac_resp = min(one(T), max(zero(T), frac_resp))
                    cpool_to_livecrootc_resp[p]         = frac_resp * cpool_to_livecrootc[p]
                    cpool_to_livecrootc_storage_resp[p] = frac_resp * cpool_to_livecrootc_storage[p]
                end
                if livestemcn_actual > livewdcn_max
                    frac_resp = (livestemcn_actual - livewdcn_max) / T(10.0)
                    frac_resp = min(one(T), max(zero(T), frac_resp))
                    cpool_to_livestemc_resp[p]         = frac_resp * cpool_to_livestemc[p]
                    cpool_to_livestemc_storage_resp[p] = frac_resp * cpool_to_livestemc_storage[p]
                end
            end

            cpool_to_resp[p] =
                cpool_to_leafc_resp[p] + cpool_to_leafc_storage_resp[p] +
                cpool_to_frootc_resp[p] + cpool_to_frootc_storage_resp[p] +
                cpool_to_livecrootc_resp[p] + cpool_to_livecrootc_storage_resp[p] +
                cpool_to_livestemc_resp[p] + cpool_to_livestemc_storage_resp[p]

            if use_matrixcn
                matrix_Cinput[p] -= cpool_to_resp[p]
            end
        end

        # ---- matrixcn N allocation B-matrix + isotope C inputs ----
        if use_matrixcn
            psn = psnsun_to_cpool[p] + psnshade_to_cpool[p]
            if use_c13 && psn != zero(T)
                matrix_C13input[p] = plant_calloc[p] *
                    ((c13_psnsun[p] + c13_psnshade[p]) / psn)
            end
            if use_c14 && psn != zero(T)
                matrix_C14input[p] = plant_calloc[p] *
                    ((c14_psnsun[p] + c14_psnshade[p]) / psn)
            end

            npool_to_veg = npool_to_leafn[p]      + npool_to_leafn_storage[p] +
                           npool_to_frootn[p]     + npool_to_frootn_storage[p] +
                           npool_to_livestemn[p]  + npool_to_livestemn_storage[p] +
                           npool_to_deadstemn[p]  + npool_to_deadstemn_storage[p] +
                           npool_to_livecrootn[p] + npool_to_livecrootn_storage[p] +
                           npool_to_deadcrootn[p] + npool_to_deadcrootn_storage[p]
            if is_crop
                npool_to_veg += npool_to_reproductiven[p, 1] +
                                npool_to_reproductiven_storage[p, 1]
            end

            if npool_to_veg != zero(T)
                invnn = one(T) / npool_to_veg
                matrix_nalloc[p, FLEXCN_ILEAF]        = npool_to_leafn[p]          * invnn
                matrix_nalloc[p, FLEXCN_ILEAF_ST]     = npool_to_leafn_storage[p]  * invnn
                matrix_nalloc[p, FLEXCN_IFROOT]       = npool_to_frootn[p]         * invnn
                matrix_nalloc[p, FLEXCN_IFROOT_ST]    = npool_to_frootn_storage[p] * invnn
                matrix_nalloc[p, FLEXCN_ILIVESTEM]    = npool_to_livestemn[p]          * invnn
                matrix_nalloc[p, FLEXCN_ILIVESTEM_ST] = npool_to_livestemn_storage[p]  * invnn
                matrix_nalloc[p, FLEXCN_IDEADSTEM]    = npool_to_deadstemn[p]          * invnn
                matrix_nalloc[p, FLEXCN_IDEADSTEM_ST] = npool_to_deadstemn_storage[p]  * invnn
                matrix_nalloc[p, FLEXCN_ILIVECROOT]    = npool_to_livecrootn[p]         * invnn
                matrix_nalloc[p, FLEXCN_ILIVECROOT_ST] = npool_to_livecrootn_storage[p] * invnn
                matrix_nalloc[p, FLEXCN_IDEADCROOT]    = npool_to_deadcrootn[p]         * invnn
                matrix_nalloc[p, FLEXCN_IDEADCROOT_ST] = npool_to_deadcrootn_storage[p] * invnn
                if is_crop
                    matrix_nalloc[p, FLEXCN_IGRAIN]    = npool_to_reproductiven[p, 1]         * invnn
                    matrix_nalloc[p, FLEXCN_IGRAIN_ST] = npool_to_reproductiven_storage[p, 1] * invnn
                end
                matrix_Ninput[p] = npool_to_veg - retransn_to_npool[p]
            end
        end
    end
end

# ---------------------------------------------------------------------------
# calc_plant_nutrient_competition_flexiblecn! -- public competition entry.
# Wraps calc_plant_cn_alloc_flexiblecn!.
# ---------------------------------------------------------------------------

"""
    calc_plant_nutrient_competition_flexiblecn!(mask_soilp, bounds, ...)

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
# calc_plant_cn_alloc_flexiblecn! -- C/N allocation after competition.
# Ported from calc_plant_cn_alloc (Fortran 193-929). Runs as one per-patch KA
# kernel (byte-identical on the KA CPU backend).
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

    out = _FlexAllocOut(;
        sminn_to_npool = nf.sminn_to_npool_patch, plant_nalloc = nf.plant_nalloc_patch,
        plant_calloc = cf.plant_calloc_patch,
        cpool_to_leafc = cf.cpool_to_leafc_patch, cpool_to_leafc_storage = cf.cpool_to_leafc_storage_patch,
        cpool_to_frootc = cf.cpool_to_frootc_patch, cpool_to_frootc_storage = cf.cpool_to_frootc_storage_patch,
        cpool_to_livestemc = cf.cpool_to_livestemc_patch, cpool_to_livestemc_storage = cf.cpool_to_livestemc_storage_patch,
        cpool_to_deadstemc = cf.cpool_to_deadstemc_patch, cpool_to_deadstemc_storage = cf.cpool_to_deadstemc_storage_patch,
        cpool_to_livecrootc = cf.cpool_to_livecrootc_patch, cpool_to_livecrootc_storage = cf.cpool_to_livecrootc_storage_patch,
        cpool_to_deadcrootc = cf.cpool_to_deadcrootc_patch, cpool_to_deadcrootc_storage = cf.cpool_to_deadcrootc_storage_patch,
        cpool_to_reproductivec = cf.cpool_to_reproductivec_patch, cpool_to_reproductivec_storage = cf.cpool_to_reproductivec_storage_patch,
        cpool_to_gresp_storage = cf.cpool_to_gresp_storage_patch,
        npool_to_leafn = nf.npool_to_leafn_patch, npool_to_leafn_storage = nf.npool_to_leafn_storage_patch,
        npool_to_frootn = nf.npool_to_frootn_patch, npool_to_frootn_storage = nf.npool_to_frootn_storage_patch,
        npool_to_livestemn = nf.npool_to_livestemn_patch, npool_to_livestemn_storage = nf.npool_to_livestemn_storage_patch,
        npool_to_deadstemn = nf.npool_to_deadstemn_patch, npool_to_deadstemn_storage = nf.npool_to_deadstemn_storage_patch,
        npool_to_livecrootn = nf.npool_to_livecrootn_patch, npool_to_livecrootn_storage = nf.npool_to_livecrootn_storage_patch,
        npool_to_deadcrootn = nf.npool_to_deadcrootn_patch, npool_to_deadcrootn_storage = nf.npool_to_deadcrootn_storage_patch,
        npool_to_reproductiven = nf.npool_to_reproductiven_patch, npool_to_reproductiven_storage = nf.npool_to_reproductiven_storage_patch,
        cpool_to_resp = cf.cpool_to_resp_patch,
        cpool_to_leafc_resp = cf.cpool_to_leafc_resp_patch, cpool_to_leafc_storage_resp = cf.cpool_to_leafc_storage_resp_patch,
        cpool_to_frootc_resp = cf.cpool_to_frootc_resp_patch, cpool_to_frootc_storage_resp = cf.cpool_to_frootc_storage_resp_patch,
        cpool_to_livecrootc_resp = cf.cpool_to_livecrootc_resp_patch, cpool_to_livecrootc_storage_resp = cf.cpool_to_livecrootc_storage_resp_patch,
        cpool_to_livestemc_resp = cf.cpool_to_livestemc_resp_patch, cpool_to_livestemc_storage_resp = cf.cpool_to_livestemc_storage_resp_patch,
        matrix_Cinput = cf.matrix_Cinput_patch, matrix_Ninput = nf.matrix_Ninput_patch,
        matrix_C13input = cf.matrix_C13input_patch, matrix_C14input = cf.matrix_C14input_patch,
        matrix_alloc = cf.matrix_alloc_patch, matrix_nalloc = nf.matrix_nalloc_patch)

    in = _FlexAllocIn(;
        froot_leaf = pftcon.froot_leaf, croot_stem = pftcon.croot_stem, stem_leaf = pftcon.stem_leaf,
        flivewd = pftcon.flivewd, grperc = pftcon.grperc, grpnow = pftcon.grpnow,
        fcur_arr = pftcon.fcur, evergreen = pftcon.evergreen, woody = pftcon.woody,
        leafcn = pftcon.leafcn, frootcn = pftcon.frootcn, livewdcn = pftcon.livewdcn,
        deadwdcn = pftcon.deadwdcn, graincn = pftcon.graincn,
        annsum_npp = cf.annsum_npp_patch, availc = cf.availc_patch, npp_growth = cf.npp_growth_patch,
        psnsun_to_cpool = cf.psnsun_to_cpool_patch, psnshade_to_cpool = cf.psnshade_to_cpool_patch,
        sminn_to_plant_fun = nf.sminn_to_plant_fun_patch, plant_ndemand = nf.plant_ndemand_patch,
        retransn_to_npool = nf.retransn_to_npool_patch, fpg_col = fpg_col,
        leafc_storage = cs.leafc_storage_patch, frootc_storage = cs.frootc_storage_patch,
        livestemc_storage = cs.livestemc_storage_patch, livecrootc_storage = cs.livecrootc_storage_patch,
        npool = ns.npool_patch, leafn_storage = ns.leafn_storage_patch, frootn_storage = ns.frootn_storage_patch,
        livestemn_storage = ns.livestemn_storage_patch, livecrootn_storage = ns.livecrootn_storage_patch,
        c_allometry = vs.c_allometry_patch, aleaf = vs.aleaf_patch, aroot = vs.aroot_patch, astem = vs.astem_patch,
        arepr = vs.arepr_patch, aleaf_n = vs.aleaf_n_patch, astem_n = vs.astem_n_patch, aroot_n = vs.aroot_n_patch,
        arepr_n = vs.arepr_n_patch, laisun = canopystate.laisun_patch, laisha = canopystate.laisha_patch,
        croplive = crop.croplive_patch)

    # mask + integer index vectors onto the state backend; reference array carries
    # the backend + eltype. The optional isotope psn arrays fall back to the base
    # cf psn arrays when the isotope struct is absent (guarded by use_c13/use_c14,
    # so the placeholder is never read).
    ref = out.plant_calloc
    FT = eltype(ref)
    nd = length(mask_soilp)
    mask_d   = similar(ref, Bool, nd); copyto!(mask_d, collect(Bool, mask_soilp))
    column_d = similar(ref, Int, length(patch.column)); copyto!(column_d, patch.column)
    itype_d  = similar(ref, Int, length(patch.itype));  copyto!(itype_d, patch.itype)
    c13_sun = (c13_cnveg_cf !== nothing) ? c13_cnveg_cf.psnsun_to_cpool_patch   : cf.psnsun_to_cpool_patch
    c13_sha = (c13_cnveg_cf !== nothing) ? c13_cnveg_cf.psnshade_to_cpool_patch : cf.psnshade_to_cpool_patch
    c14_sun = (c14_cnveg_cf !== nothing) ? c14_cnveg_cf.psnsun_to_cpool_patch   : cf.psnsun_to_cpool_patch
    c14_sha = (c14_cnveg_cf !== nothing) ? c14_cnveg_cf.psnshade_to_cpool_patch : cf.psnshade_to_cpool_patch

    backend = _kernel_backend(ref)
    _flexalloc_main_kernel!(backend)(out, in, mask_d, column_d, itype_d,
        c13_sun, c13_sha, c14_sun, c14_sha,
        use_fun, use_matrixcn, use_crop_agsys, carbon_resp_opt,
        use_c13, use_c14, npcropmin, nrepr, FT(dt); ndrange = nd)
    KA.synchronize(backend)

    return nothing
end

# ===========================================================================
# N demand: per-patch kernels (parallel the nutrient_competition.jl demand set)
# ===========================================================================

# --- Loop 1: FlexibleCN plant N demand (Michaelis-Menten) + accumulation ----
@kernel function _flexdemand_main_kernel!(plant_ndemand, tempsum_potential_gpp, tempmax_retransn,
        @Const(mask), @Const(column), @Const(itype),
        @Const(leafn), @Const(leafc), @Const(leafcn_pft), @Const(frootc),
        @Const(availc), @Const(n_allometry), @Const(c_allometry),
        @Const(laisun), @Const(laisha), @Const(sminn_vr), @Const(dzsoi_decomp),
        @Const(t_scalar), @Const(gpp_before_downreg), @Const(retransn),
        use_fun::Bool, nlevdecomp::Int, Vmax_N, Kmin, spval, dt)
    T = eltype(plant_ndemand)
    p = @index(Global)
    @inbounds if mask[p]
        c = column[p]
        ivt = itype[p] + 1

        if leafn[p] < T(NUTRIENT_COMP_N_MIN)
            actual_leafcn = spval
        else
            actual_leafcn = leafc[p] / leafn[p]
        end
        leafcn_min = leafcn_pft[ivt] - T(10.0)
        leafcn_max = leafcn_pft[ivt] + T(10.0)
        actual_leafcn = max(actual_leafcn, leafcn_min - T(0.0001))
        actual_leafcn = min(actual_leafcn, leafcn_max)

        nscalar = (actual_leafcn - leafcn_min) / (leafcn_max - leafcn_min)
        nscalar = min(max(zero(T), nscalar), one(T))

        sminn_total = zero(T)
        for j in 1:nlevdecomp
            sminn_total += sminn_vr[c, j] * dzsoi_decomp[j]
        end
        substrate_term = sminn_total / (sminn_total + Kmin)

        temp_scalar = t_scalar[c, 1]
        temp_scalar = min(max(zero(T), temp_scalar), one(T))

        if use_fun
            plant_ndemand[p] = availc[p] * (n_allometry[p] / c_allometry[p])
        else
            lai_tot = laisun[p] + laisha[p]
            if lai_tot > zero(T)
                plant_ndemand[p] = Vmax_N * frootc[p] * substrate_term * temp_scalar * nscalar
            else
                plant_ndemand[p] = zero(T)
            end
            if actual_leafcn < leafcn_min
                plant_ndemand[p] = zero(T)
            end
        end

        tempsum_potential_gpp[p] += gpp_before_downreg[p]
        tempmax_retransn[p] = max(tempmax_retransn[p], retransn[p])
    end
end

# --- Loop 2: crop grain-fill retranslocation (pool-relative FlexibleCN form) -
@kernel function _flexdemand_grainfill_kernel!(grain_flag,
        leafn_to_retransn, livestemn_to_retransn, frootn_to_retransn,
        @Const(mask), @Const(itype), @Const(croplive), @Const(crop_phase_vals),
        @Const(astem), @Const(astemf),
        @Const(leafn), @Const(leafc), @Const(livestemn), @Const(livestemc),
        @Const(frootn), @Const(frootc), @Const(fleafcn), @Const(fstemcn), @Const(ffrootcn),
        cphase_leafemerge_in, cphase_grainfill_in,
        ntmp_soybean::Int, nirrig_tmp_soybean::Int, ntrp_soybean::Int, nirrig_trp_soybean::Int, dt)
    T = eltype(grain_flag)
    p = @index(Global)
    @inbounds if mask[p]
        if croplive[p] != zero(eltype(croplive))
            ivt = itype[p] + 1
            cphase = crop_phase_vals[p]
            if cphase == cphase_leafemerge_in
                grain_flag[p] = zero(T)
            elseif cphase == cphase_grainfill_in
                is_soybean = (ivt == ntmp_soybean || ivt == nirrig_tmp_soybean ||
                              ivt == ntrp_soybean || ivt == nirrig_trp_soybean)
                if astem[p] == astemf[ivt] || !is_soybean
                    if grain_flag[p] == zero(T)
                        t1 = one(T) / dt
                        leafn_to_retransn[p] = t1 *
                            max(leafn[p] - (leafc[p] / fleafcn[ivt]), zero(T))
                        livestemn_to_retransn[p] = t1 *
                            max(livestemn[p] - (livestemc[p] / fstemcn[ivt]), zero(T))
                        frootn_to_retransn[p] = zero(T)
                        if ffrootcn[ivt] > zero(T)
                            frootn_to_retransn[p] = t1 *
                                max(frootn[p] - (frootc[p] / ffrootcn[ivt]), zero(T))
                        end
                        grain_flag[p] = one(T)
                    end
                end
            end
        end
    end
end

# --- Loop 3a: avail_retransn for prognostic-crop call (grain-flag gated) -----
@kernel function _flexdemand_availretransn_pcrop_kernel!(avail_retransn,
        @Const(mask), @Const(grain_flag), @Const(plant_ndemand))
    T = eltype(avail_retransn)
    p = @index(Global)
    @inbounds if mask[p]
        if grain_flag[p] == one(T)
            avail_retransn[p] = plant_ndemand[p]
        else
            avail_retransn[p] = zero(T)
        end
    end
end

# --- Loop 3b: avail_retransn for non-pcrop call (seasonal GPP fraction) ------
@kernel function _flexdemand_availretransn_nopcrop_kernel!(avail_retransn,
        @Const(mask), @Const(annsum_potential_gpp), @Const(annmax_retransn),
        @Const(gpp_before_downreg), dt)
    T = eltype(avail_retransn)
    p = @index(Global)
    @inbounds if mask[p]
        if annsum_potential_gpp[p] > zero(T)
            avail_retransn[p] = (annmax_retransn[p] / T(2.0)) *
                (gpp_before_downreg[p] / annsum_potential_gpp[p]) / dt
        else
            avail_retransn[p] = zero(T)
        end
    end
end

# --- Loop 4: clamp avail_retransn to storage + modify plant N demand ---------
@kernel function _flexdemand_retransn_modify_kernel!(avail_retransn, retransn_to_npool, plant_ndemand,
        @Const(mask), @Const(retransn), use_fun::Bool, dt)
    T = eltype(avail_retransn)
    p = @index(Global)
    @inbounds if mask[p]
        avail_retransn[p] = min(avail_retransn[p], retransn[p] / dt)
        if plant_ndemand[p] > avail_retransn[p]
            retransn_to_npool[p] = avail_retransn[p]
        else
            retransn_to_npool[p] = plant_ndemand[p]
        end
        if !use_fun
            plant_ndemand[p] -= retransn_to_npool[p]
        end
    end
end

# ---------------------------------------------------------------------------
# calc_plant_nutrient_demand_flexiblecn! -- public N-demand entry.
# ---------------------------------------------------------------------------

"""
    calc_plant_nutrient_demand_flexiblecn!(mask_p, bounds, call_is_for_pcrop, ...)

FlexibleCN plant nitrogen demand. Parallel of `calc_plant_nutrient_demand!`.
Plant N demand is computed from Michaelis-Menten kinetics on fine-root C, a
soil mineral-N substrate term, a soil-temperature scalar and a flexible
leaf-C:N nitrogen scalar. Runs as per-patch KA kernels.

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
    FT = eltype(nf.plant_ndemand_patch)
    Vmax_N = FT(2.7e-8); Kmin = FT(1.0); SPVAL = FT(1.0e36)

    # Loop 1: FlexibleCN plant N demand (Michaelis-Menten) + accumulation.
    # dzsoi_decomp is a host physical-constant vector → move to the state backend.
    dz_k = _to_backend_like(nf.plant_ndemand_patch, FT, dzsoi_decomp)
    _launch!(_flexdemand_main_kernel!, nf.plant_ndemand_patch,
        vs.tempsum_potential_gpp_patch, vs.tempmax_retransn_patch,
        mask_p, patch.column, patch.itype,
        ns.leafn_patch, cs.leafc_patch, pftcon.leafcn, cs.frootc_patch,
        cf.availc_patch, vs.n_allometry_patch, vs.c_allometry_patch,
        canopystate.laisun_patch, canopystate.laisha_patch,
        soilbgc_ns.sminn_vr_col, dz_k, soilbgc_cf.t_scalar_col,
        cf.gpp_before_downreg_patch, ns.retransn_patch,
        use_fun, nlevdecomp, Vmax_N, Kmin, SPVAL, FT(dt))

    # Loop 2: crop grain-fill retranslocation
    if call_is_for_pcrop
        crop_phase_vals = fill!(similar(nf.plant_ndemand_patch, FT, length(mask_p)), zero(FT))
        crop_phase!(mask_p, crop, vs, crop_phase_vals)

        _launch!(_flexdemand_grainfill_kernel!, vs.grain_flag_patch,
            nf.leafn_to_retransn_patch, nf.livestemn_to_retransn_patch, nf.frootn_to_retransn_patch,
            mask_p, patch.itype, crop.croplive_patch, crop_phase_vals,
            vs.astem_patch, pftcon.astemf,
            ns.leafn_patch, cs.leafc_patch, ns.livestemn_patch, cs.livestemc_patch,
            ns.frootn_patch, cs.frootc_patch, pftcon.fleafcn, pftcon.fstemcn, pftcon.ffrootcn,
            FT(cphase_leafemerge), FT(cphase_grainfill),
            ntmp_soybean, nirrig_tmp_soybean, ntrp_soybean, nirrig_trp_soybean, FT(dt))
    end

    # Loop 3: crops pull from retransn only during grain fill; others seasonal fraction
    if call_is_for_pcrop
        _launch!(_flexdemand_availretransn_pcrop_kernel!, nf.avail_retransn_patch,
            mask_p, vs.grain_flag_patch, nf.plant_ndemand_patch)
    else
        _launch!(_flexdemand_availretransn_nopcrop_kernel!, nf.avail_retransn_patch,
            mask_p, vs.annsum_potential_gpp_patch, vs.annmax_retransn_patch,
            cf.gpp_before_downreg_patch, FT(dt))
    end

    # Loop 4: clamp to storage + modify plant N demand
    _launch!(_flexdemand_retransn_modify_kernel!,
        nf.avail_retransn_patch, nf.retransn_to_npool_patch, nf.plant_ndemand_patch,
        mask_p, ns.retransn_patch, use_fun, FT(dt))

    return nothing
end
