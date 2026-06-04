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
Base.@kwdef mutable struct PftConAllocation{FT<:Real, V<:AbstractVector{FT}}
    woody       ::V = Float64[]   # binary woody flag (1=woody, 0=not woody)
    froot_leaf  ::V = Float64[]   # new fine root C per new leaf C (gC/gC)
    croot_stem  ::V = Float64[]   # new coarse root C per new stem C (gC/gC)
    stem_leaf   ::V = Float64[]   # new stem C per new leaf C (gC/gC); -1 = dynamic
    flivewd     ::V = Float64[]   # fraction of new wood that is live
    leafcn      ::V = Float64[]   # leaf C:N (gC/gN)
    frootcn     ::V = Float64[]   # fine root C:N (gC/gN)
    livewdcn    ::V = Float64[]   # live wood C:N (gC/gN)
    deadwdcn    ::V = Float64[]   # dead wood C:N (gC/gN)
    graincn     ::V = Float64[]   # grain C:N (gC/gN)
    grperc      ::V = Float64[]   # growth respiration fraction
    # Crop-specific allocation parameters
    arooti      ::V = Float64[]   # initial allocation to roots
    arootf      ::V = Float64[]   # final allocation to roots
    bfact       ::V = Float64[]   # exp factor for leaf allocation
    fleafi      ::V = Float64[]   # initial fraction allocated to leaf
    aleaff      ::V = Float64[]   # final leaf allocation coefficient
    astemf      ::V = Float64[]   # final stem allocation coefficient
    allconss    ::V = Float64[]   # stem allocation decline exponent
    allconsl    ::V = Float64[]   # leaf allocation decline exponent
    declfact    ::V = Float64[]   # decline factor for allocation shift
end
PftConAllocation{FT}(; kwargs...) where {FT<:Real} = PftConAllocation{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure PftConAllocation

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
# --- Device-view bundles for _alloc_gpp_mr_availc_kernel! (Metal arg-count safe) ---

# Output arrays (mixes patch-vectors + reproductive matrices) -> {V,M}
Base.@kwdef struct _AllocGppOut{V,M}
    psnsun_to_cpool::V
    psnshade_to_cpool::V
    c13_psnsun_to_cpool::V
    c13_psnshade_to_cpool::V
    c14_psnsun_to_cpool::V
    c14_psnshade_to_cpool::V
    gpp_before_downreg::V
    availc::V
    leaf_curmr::V
    leaf_xsmr::V
    froot_curmr::V
    froot_xsmr::V
    livestem_curmr::V
    livestem_xsmr::V
    livecroot_curmr::V
    livecroot_xsmr::V
    reproductive_curmr::M
    reproductive_xsmr::M
    xsmrpool_recover_out::V
    cpool_to_xsmrpool::V
end
Adapt.@adapt_structure _AllocGppOut

# Read-only input arrays (mixes patch-vectors + reproductive_mr matrix) -> {V,M}
Base.@kwdef struct _AllocGppIn{V,M}
    psnsun_patch::V
    psnsha_patch::V
    laisun_patch::V
    laisha_patch::V
    c13_psnsun_patch::V
    c13_psnsha_patch::V
    c14_psnsun_patch::V
    c14_psnsha_patch::V
    woody::V
    leaf_mr::V
    froot_mr::V
    livestem_mr::V
    livecroot_mr::V
    reproductive_mr::M
    xsmrpool::V
end
Adapt.@adapt_structure _AllocGppIn

# isbits scalar params at working precision
Base.@kwdef struct _AllocGppScalars{T}
    dayscrecover::T
    secspday::T
end

@kernel function _alloc_gpp_mr_availc_kernel!(
        out::_AllocGppOut, in::_AllocGppIn,
        @Const(mask_soilp), @Const(itype), @Const(croplive),
        scal::_AllocGppScalars,
        do_c13::Bool, do_c14::Bool, use_matrixcn::Bool,
        npcropmin::Int, nrepr::Int)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        T = eltype(out.availc)

        # alias every field to a Fortran-named local; body stays verbatim except literal -> T
        psnsun_to_cpool      = out.psnsun_to_cpool
        psnshade_to_cpool    = out.psnshade_to_cpool
        c13_psnsun_to_cpool  = out.c13_psnsun_to_cpool
        c13_psnshade_to_cpool = out.c13_psnshade_to_cpool
        c14_psnsun_to_cpool  = out.c14_psnsun_to_cpool
        c14_psnshade_to_cpool = out.c14_psnshade_to_cpool
        gpp_before_downreg   = out.gpp_before_downreg
        availc               = out.availc
        leaf_curmr           = out.leaf_curmr
        leaf_xsmr            = out.leaf_xsmr
        froot_curmr          = out.froot_curmr
        froot_xsmr           = out.froot_xsmr
        livestem_curmr       = out.livestem_curmr
        livestem_xsmr        = out.livestem_xsmr
        livecroot_curmr      = out.livecroot_curmr
        livecroot_xsmr       = out.livecroot_xsmr
        reproductive_curmr   = out.reproductive_curmr
        reproductive_xsmr    = out.reproductive_xsmr
        xsmrpool_recover_out = out.xsmrpool_recover_out
        cpool_to_xsmrpool    = out.cpool_to_xsmrpool

        psnsun_patch    = in.psnsun_patch
        psnsha_patch    = in.psnsha_patch
        laisun_patch    = in.laisun_patch
        laisha_patch    = in.laisha_patch
        c13_psnsun_patch = in.c13_psnsun_patch
        c13_psnsha_patch = in.c13_psnsha_patch
        c14_psnsun_patch = in.c14_psnsun_patch
        c14_psnsha_patch = in.c14_psnsha_patch
        woody           = in.woody
        leaf_mr         = in.leaf_mr
        froot_mr        = in.froot_mr
        livestem_mr     = in.livestem_mr
        livecroot_mr    = in.livecroot_mr
        reproductive_mr = in.reproductive_mr
        xsmrpool        = in.xsmrpool

        dayscrecover = scal.dayscrecover
        secspday     = scal.secspday

        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia

        # Convert psn from umol/m2/s -> gC/m2/s
        # The input psn (psnsun and psnsha) are expressed per unit LAI
        # in the sunlit and shaded canopy. Scale by laisun and laisha.
        psnsun_to_cpool[p]   = psnsun_patch[p] * laisun_patch[p] * T(12.011e-6)
        psnshade_to_cpool[p] = psnsha_patch[p] * laisha_patch[p] * T(12.011e-6)

        if do_c13
            c13_psnsun_to_cpool[p]   = c13_psnsun_patch[p] * laisun_patch[p] * T(12.011e-6)
            c13_psnshade_to_cpool[p] = c13_psnsha_patch[p] * laisha_patch[p] * T(12.011e-6)
        end

        if do_c14
            c14_psnsun_to_cpool[p]   = c14_psnsun_patch[p] * laisun_patch[p] * T(12.011e-6)
            c14_psnshade_to_cpool[p] = c14_psnsha_patch[p] * laisha_patch[p] * T(12.011e-6)
        end

        gpp = psnsun_to_cpool[p] + psnshade_to_cpool[p]
        gpp_before_downreg[p] = gpp

        # get the time step total maintenance respiration
        mr = leaf_mr[p] + froot_mr[p]
        if woody[ivt] == one(T)
            mr = mr + livestem_mr[p] + livecroot_mr[p]
        elseif ivt >= npcropmin
            if croplive[p]
                reproductive_mr_tot = zero(T)
                for k in 1:nrepr
                    reproductive_mr_tot = reproductive_mr_tot + reproductive_mr[p, k]
                end
                mr = mr + livestem_mr[p] + reproductive_mr_tot
            end
        end

        # For Matrix solution if mr is very small set it to zero
        if mr < T(-1.0e-15) && use_matrixcn
            mr = zero(T)
        end

        # carbon flux available for allocation
        availc_p = gpp - mr
        availc[p] = availc_p

        # If mr > gpp, then some mr comes from gpp, the rest from cpool (xsmr)
        if mr > zero(T) && availc_p < zero(T)
            curmr = gpp
            curmr_ratio = curmr / mr
        else
            curmr_ratio = one(T)
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
        availc[p] = max(availc[p], zero(T))

        # test for an xsmrpool deficit
        if xsmrpool[p] < zero(T)
            xsmrpool_recover = -xsmrpool[p] / (dayscrecover * secspday)
            if xsmrpool_recover < availc[p]
                availc[p] = availc[p] - xsmrpool_recover
            else
                xsmrpool_recover = availc[p]
                availc[p] = zero(T)
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

    # Working precision (from an output array element type)
    T = eltype(cnveg_cf.availc_patch)

    out = _AllocGppOut(;
        psnsun_to_cpool       = cnveg_cf.psnsun_to_cpool_patch,
        psnshade_to_cpool     = cnveg_cf.psnshade_to_cpool_patch,
        c13_psnsun_to_cpool   = c13_psnsun_to_cpool,
        c13_psnshade_to_cpool = c13_psnshade_to_cpool,
        c14_psnsun_to_cpool   = c14_psnsun_to_cpool,
        c14_psnshade_to_cpool = c14_psnshade_to_cpool,
        gpp_before_downreg    = cnveg_cf.gpp_before_downreg_patch,
        availc                = cnveg_cf.availc_patch,
        leaf_curmr            = cnveg_cf.leaf_curmr_patch,
        leaf_xsmr             = cnveg_cf.leaf_xsmr_patch,
        froot_curmr           = cnveg_cf.froot_curmr_patch,
        froot_xsmr            = cnveg_cf.froot_xsmr_patch,
        livestem_curmr        = cnveg_cf.livestem_curmr_patch,
        livestem_xsmr         = cnveg_cf.livestem_xsmr_patch,
        livecroot_curmr       = cnveg_cf.livecroot_curmr_patch,
        livecroot_xsmr        = cnveg_cf.livecroot_xsmr_patch,
        reproductive_curmr    = cnveg_cf.reproductive_curmr_patch,
        reproductive_xsmr     = cnveg_cf.reproductive_xsmr_patch,
        xsmrpool_recover_out  = cnveg_cf.xsmrpool_recover_patch,
        cpool_to_xsmrpool     = cnveg_cf.cpool_to_xsmrpool_patch)

    in = _AllocGppIn(;
        psnsun_patch     = photosyns.psnsun_patch,
        psnsha_patch     = photosyns.psnsha_patch,
        laisun_patch     = canopystate.laisun_patch,
        laisha_patch     = canopystate.laisha_patch,
        c13_psnsun_patch = photosyns.c13_psnsun_patch,
        c13_psnsha_patch = photosyns.c13_psnsha_patch,
        c14_psnsun_patch = photosyns.c14_psnsun_patch,
        c14_psnsha_patch = photosyns.c14_psnsha_patch,
        woody            = pftcon.woody,
        leaf_mr          = cnveg_cf.leaf_mr_patch,
        froot_mr         = cnveg_cf.froot_mr_patch,
        livestem_mr      = cnveg_cf.livestem_mr_patch,
        livecroot_mr     = cnveg_cf.livecroot_mr_patch,
        reproductive_mr  = cnveg_cf.reproductive_mr_patch,
        xsmrpool         = cnveg_cs.xsmrpool_patch)

    scal = _AllocGppScalars(; dayscrecover = T(dayscrecover), secspday = T(SECSPDAY))

    backend = _kernel_backend(out.availc)
    _alloc_gpp_mr_availc_kernel!(backend)(
        out, in,
        mask_soilp, patch.itype, crop.croplive_patch,
        scal,
        do_c13, do_c14, use_matrixcn,
        npcropmin, nrepr;
        ndrange = length(mask_soilp))
    KA.synchronize(backend)

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
@kernel function _alloc_crop_fractions_kernel!(
        aleaf, astem, aroot, aleafi, astemi, arepr,
        @Const(mask_pcropp), @Const(itype), @Const(croplive),
        @Const(crop_phase_out), @Const(peaklai),
        @Const(hui), @Const(gddmaturity), @Const(huigrain),
        @Const(arooti), @Const(arootf), @Const(fleafi), @Const(bfact),
        @Const(astemf), @Const(declfact), @Const(allconss),
        @Const(aleaff), @Const(allconsl),
        cphase_planted_in, cphase_leafemerge_in, cphase_grainfill_in,
        nrepr::Int)
    T = eltype(aleaf)
    p = @index(Global)
    @inbounds if mask_pcropp[p]
        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia

        if croplive[p]
            # Phase 1 completed: leaf emergence to start of leaf decline
            if crop_phase_out[p] == cphase_leafemerge

                for k in 1:nrepr
                    arepr[p, k] = zero(T)
                end
                if peaklai[p] == 1  # lai at maximum allowed
                    aleaf[p] = T(1.0e-5)
                    astem[p] = zero(T)
                    aroot[p] = one(T) - aleaf[p]
                else
                    aroot[p] = max(zero(T), min(one(T), arooti[ivt] -
                        (arooti[ivt] - arootf[ivt]) *
                        min(one(T), hui[p] / gddmaturity[p])))
                    fleaf = fleafi[ivt] * (exp(-bfact[ivt]) -
                        exp(-bfact[ivt] * hui[p] / huigrain[p])) /
                        (exp(-bfact[ivt]) - one(T))
                    aleaf[p] = max(T(1.0e-5), (one(T) - aroot[p]) * fleaf)
                    astem[p] = one(T) - aleaf[p] - aroot[p]
                end

                astemi[p] = astem[p]
                aleafi[p] = aleaf[p]

            # Phase 2 completed: grain fill
            elseif crop_phase_out[p] == cphase_grainfill
                aroot[p] = max(zero(T), min(one(T), arooti[ivt] -
                    (arooti[ivt] - arootf[ivt]) *
                    min(one(T), hui[p] / gddmaturity[p])))
                if astemi[p] > astemf[ivt]
                    astem[p] = max(zero(T), max(astemf[ivt],
                        astem[p] *
                        (one(T) - min((hui[p] -
                        huigrain[p]) / ((gddmaturity[p] * declfact[ivt]) -
                        huigrain[p]), one(T))^allconss[ivt])))
                end

                if peaklai[p] == 1
                    aleaf[p] = T(1.0e-5)
                elseif aleafi[p] > aleaff[ivt]
                    aleaf[p] = max(T(1.0e-5), max(aleaff[ivt],
                        aleaf[p] *
                        (one(T) - min((hui[p] -
                        huigrain[p]) / ((gddmaturity[p] * declfact[ivt]) -
                        huigrain[p]), one(T))^allconsl[ivt])))
                end

                # All repr allocation goes into the last reproductive pool
                for k in 1:nrepr-1
                    arepr[p, k] = zero(T)
                end
                arepr[p, nrepr] = one(T) - aroot[p] - astem[p] - aleaf[p]

            else  # cphase_planted: pre emergence
                aleaf[p]  = one(T)
                aleafi[p] = one(T)
                astem[p]  = zero(T)
                astemi[p] = zero(T)
                aroot[p]  = zero(T)
                for k in 1:nrepr
                    arepr[p, k] = zero(T)
                end
            end

        else  # .not. croplive
            aleaf[p]  = one(T)
            aleafi[p] = one(T)
            astem[p]  = zero(T)
            astemi[p] = zero(T)
            aroot[p]  = zero(T)
            for k in 1:nrepr
                arepr[p, k] = zero(T)
            end
        end
    end
end

function calc_crop_allocation_fractions!(mask_pcropp::AbstractVector{Bool}, bounds::UnitRange{Int},
                                          pftcon::PftConAllocation,
                                          patch::PatchData,
                                          crop::CropData,
                                          cnveg_state::CNVegStateData;
                                          nrepr::Int=NREPR)

    # Compute crop phase for all patches in the mask
    FT = eltype(cnveg_state.c_allometry_patch)
    crop_phase_out = similar(cnveg_state.c_allometry_patch, length(mask_pcropp))
    crop_phase!(mask_pcropp, crop, cnveg_state, crop_phase_out)

    # Crop-phase code constants at working precision (host-resolved scalars).
    cphase_planted_in    = FT(cphase_planted)
    cphase_leafemerge_in = FT(cphase_leafemerge)
    cphase_grainfill_in  = FT(cphase_grainfill)

    # Hoisted error guard: the original loop throws on an unexpected crop_phase
    # code for a live crop. The kernel below is error-free, so replicate the
    # check on the host first (never fires for valid data).
    for p in bounds
        mask_pcropp[p] || continue
        crop.croplive_patch[p] || continue
        cp = crop_phase_out[p]
        if cp != cphase_leafemerge_in && cp != cphase_grainfill_in && cp != cphase_planted_in
            error("Unexpected crop_phase: $(cp) at patch $p")
        end
    end

    _launch!(_alloc_crop_fractions_kernel!,
        cnveg_state.aleaf_patch, cnveg_state.astem_patch, cnveg_state.aroot_patch,
        cnveg_state.aleafi_patch, cnveg_state.astemi_patch, cnveg_state.arepr_patch,
        mask_pcropp, patch.itype, crop.croplive_patch,
        crop_phase_out, cnveg_state.peaklai_patch,
        crop.hui_patch, cnveg_state.gddmaturity_patch, cnveg_state.huigrain_patch,
        pftcon.arooti, pftcon.arootf, pftcon.fleafi, pftcon.bfact,
        pftcon.astemf, pftcon.declfact, pftcon.allconss,
        pftcon.aleaff, pftcon.allconsl,
        cphase_planted_in, cphase_leafemerge_in, cphase_grainfill_in,
        nrepr)

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
    T = eltype(c_allometry)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia

        f1 = froot_leaf[ivt]
        f2 = croot_stem[ivt]

        # modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0
        if stem_leaf[ivt] == T(-1.0)
            f3 = (T(2.7) / (one(T) + exp(T(-0.004) * (annsum_npp[p] - T(300.0))))) - T(0.4)
        else
            f3 = stem_leaf[ivt]
        end

        f4  = flivewd[ivt]
        if ivt >= npcropmin
            g1 = T(0.25)
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
            g1a = zero(T)
        end

        if woody[ivt] == one(T)
            c_allometry[p] = (one(T) + g1a) * (one(T) + f1 + f3 * (one(T) + f2))
            n_allometry[p] = one(T) / cnl + f1 / cnfr +
                (f3 * f4 * (one(T) + f2)) / cnlw +
                (f3 * (one(T) - f4) * (one(T) + f2)) / cndw
        elseif ivt >= npcropmin  # skip generic crops
            cng = graincn[ivt]
            f1 = aroot[p] / aleaf[p]
            f3 = astem[p] / aleaf[p]
            f5_tot   = zero(T)
            f5_n_tot = zero(T)
            for k in 1:nrepr
                f5k = arepr[p, k] / aleaf[p]
                f5_tot   = f5_tot + f5k
                f5_n_tot = f5_n_tot + f5k / cng
            end
            c_allometry[p] = (one(T) + g1a) * (one(T) + f1 + f5_tot + f3 * (one(T) + f2))
            n_allometry[p] = one(T) / cnl + f1 / cnfr + f5_n_tot +
                (f3 * f4 * (one(T) + f2)) / cnlw +
                (f3 * (one(T) - f4) * (one(T) + f2)) / cndw
        else
            c_allometry[p] = one(T) + g1a + f1 + f1 * g1a
            n_allometry[p] = one(T) / cnl + f1 / cnfr
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
