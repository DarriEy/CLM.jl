# ==========================================================================
# Matrix-CN veg-C wiring (phase 1) — connect the (already-validated) matrix
# solver to the live CN cycle without touching the GPU-kernelized phenology/
# gap/fire modules.
#
# Approach (see memory matrix-cn-wiring-scope): the per-flux transfer fields
# (leafc_storage_to_xfer_patch, …) are computed by phenology REGARDLESS of
# use_matrixcn (they feed the c_state_update donor decrements). So instead of
# threading ~250 stateful matrix_update calls into the flux kernels, a SEPARATE
# accumulation pass reads each flux field + its donor pool, recovers the transfer
# rate (= flux/donor), and records it via matrix_update_phc! (the proven
# accumulator). Then cn_veg_matrix_solve_c! advances all pools in one solve.
#
# At matrixcheck=false the accumulation is purely additive/order-independent and
# the matrix advance reproduces the sequential c_state_update EXACTLY (proven in
# test_cn_veg_matrix). The turnover cap (matrixcheck=true) is a later refinement.
#
# This file wires the VEG-C PHENOLOGY transfers (17 natveg / 18 crop). Gap-
# mortality + fire accumulation + the N + soil paths follow the same pattern.
# ==========================================================================

# Pool index (varpar constants) -> CNVegCarbonState carbon pool value at patch p.
@inline function _vegc_pool_val(cs::CNVegCarbonStateData, i::Int, p::Int)
    @inbounds if     i == ILEAF          ; return cs.leafc_patch[p]
    elseif i == ILEAF_ST       ; return cs.leafc_storage_patch[p]
    elseif i == ILEAF_XF       ; return cs.leafc_xfer_patch[p]
    elseif i == IFROOT         ; return cs.frootc_patch[p]
    elseif i == IFROOT_ST      ; return cs.frootc_storage_patch[p]
    elseif i == IFROOT_XF      ; return cs.frootc_xfer_patch[p]
    elseif i == ILIVESTEM      ; return cs.livestemc_patch[p]
    elseif i == ILIVESTEM_ST   ; return cs.livestemc_storage_patch[p]
    elseif i == ILIVESTEM_XF   ; return cs.livestemc_xfer_patch[p]
    elseif i == IDEADSTEM      ; return cs.deadstemc_patch[p]
    elseif i == IDEADSTEM_ST   ; return cs.deadstemc_storage_patch[p]
    elseif i == IDEADSTEM_XF   ; return cs.deadstemc_xfer_patch[p]
    elseif i == ILIVECROOT     ; return cs.livecrootc_patch[p]
    elseif i == ILIVECROOT_ST  ; return cs.livecrootc_storage_patch[p]
    elseif i == ILIVECROOT_XF  ; return cs.livecrootc_xfer_patch[p]
    elseif i == IDEADCROOT     ; return cs.deadcrootc_patch[p]
    elseif i == IDEADCROOT_ST  ; return cs.deadcrootc_storage_patch[p]
    elseif i == IDEADCROOT_XF  ; return cs.deadcrootc_xfer_patch[p]
    elseif i == IGRAIN         ; return cs.reproductivec_patch[p, 1]         # crop grain (reproductive)
    elseif i == IGRAIN_ST      ; return cs.reproductivec_storage_patch[p, 1]
    elseif i == IGRAIN_XF      ; return cs.reproductivec_xfer_patch[p, 1]
    else                       ; return 0.0
    end
end

# The 17 natural-veg phenology transfers (14 in-veg storage/xfer/pool +
# livewood-turnover, then 3 trailing pure out-of-veg losses = litterfall). The
# solve treats the last ncphouttrans entries as turnover-only (no A entry), so
# the out-transfers MUST be the trailing indices. Crop adds grain->iout (18/4).
const _PHC_DONER_NATVEG = Int[
    ILEAF_ST, ILEAF_XF, IFROOT_ST, IFROOT_XF,
    ILIVESTEM_ST, ILIVESTEM_XF, IDEADSTEM_ST, IDEADSTEM_XF,
    ILIVECROOT_ST, ILIVECROOT_XF, IDEADCROOT_ST, IDEADCROOT_XF,
    ILIVESTEM, ILIVECROOT,               # livewood turnover (live->dead)
    ILEAF, IFROOT, ILIVESTEM]            # out: litterfall (leaf, froot, livestem)
const _PHC_RECEIVER_NATVEG = Int[
    ILEAF_XF, ILEAF, IFROOT_XF, IFROOT,
    ILIVESTEM_XF, ILIVESTEM, IDEADSTEM_XF, IDEADSTEM,
    ILIVECROOT_XF, ILIVECROOT, IDEADCROOT_XF, IDEADCROOT,
    IDEADSTEM, IDEADCROOT,
    NVEGPOOL_NATVEG + 1, NVEGPOOL_NATVEG + 1, NVEGPOOL_NATVEG + 1]  # iout

"""
    cn_veg_matrix_c_topology!(cf; use_crop=false)

Set up the veg-C phenology transfer topology on the flux instance: assign the
named `_ph` transfer indices (1..ncphtrans) and populate
`matrix_phtransfer_doner_patch` / `matrix_phtransfer_receiver_patch`. Idempotent;
call once after the flux instance is sized. The internal ordering is our own but
consistent (the solve reads the donor/receiver arrays), with the pure out-of-veg
losses as the trailing `ncphouttrans` indices as the solver requires.
"""
function cn_veg_matrix_c_topology!(cf::CNVegCarbonFluxData; use_crop::Bool = false,
                                   nvegcpool::Int = NVEGPOOL_NATVEG)
    doner    = copy(_PHC_DONER_NATVEG)
    receiver = copy(_PHC_RECEIVER_NATVEG)
    if use_crop
        # For crop, "outside veg" is ioutc = nvegcpool+1 (grain occupies 19..21), so
        # remap the natveg out-receivers (NVEGPOOL_NATVEG+1) to ioutc. The solver never
        # reads out-transfer receivers, but this keeps the arrays semantically correct.
        ioutc = nvegcpool + 1
        for k in eachindex(receiver)
            receiver[k] == NVEGPOOL_NATVEG + 1 && (receiver[k] = ioutc)
        end
        # crop adds grain->iout as a 4th out-transfer (index 18).
        push!(doner, IGRAIN); push!(receiver, ioutc)
        cf.igrain_to_iout_ph = length(doner)
    end
    cf.matrix_phtransfer_doner_patch    = doner
    cf.matrix_phtransfer_receiver_patch = receiver

    # Size the per-step transfer/turnover accumulators (the flux init leaves them
    # empty). np from a phenology flux field the init already sized.
    np = length(cf.leafc_to_litter_patch)
    cf.matrix_phtransfer_patch = zeros(np, length(doner))
    cf.matrix_phturnover_patch = zeros(np, nvegcpool)
    cf.matrix_alloc_patch  = zeros(np, nvegcpool)
    cf.matrix_Cinput_patch = zeros(np)

    cf.ileafst_to_ileafxf_ph        = 1
    cf.ileafxf_to_ileaf_ph          = 2
    cf.ifrootst_to_ifrootxf_ph      = 3
    cf.ifrootxf_to_ifroot_ph        = 4
    cf.ilivestemst_to_ilivestemxf_ph = 5
    cf.ilivestemxf_to_ilivestem_ph  = 6
    cf.ideadstemst_to_ideadstemxf_ph = 7
    cf.ideadstemxf_to_ideadstem_ph  = 8
    cf.ilivecrootst_to_ilivecrootxf_ph = 9
    cf.ilivecrootxf_to_ilivecroot_ph = 10
    cf.ideadcrootst_to_ideadcrootxf_ph = 11
    cf.ideadcrootxf_to_ideadcroot_ph = 12
    cf.ilivestem_to_ideadstem_ph    = 13
    cf.ilivecroot_to_ideadcroot_ph  = 14
    cf.ileaf_to_iout_ph             = 15
    cf.ifroot_to_iout_ph            = 16
    cf.ilivestem_to_iout_ph         = 17

    # --- Gap mortality: 18 pure out-of-veg losses (pool i -> iout), index i.
    # (grain does not gap/burn → gm/fi stay the 18 natveg pools even for crop.)
    iout = nvegcpool + 1
    cf.matrix_gmtransfer_doner_patch    = collect(1:NVEGPOOL_NATVEG)   # ILEAF..IDEADCROOT_XF
    cf.matrix_gmtransfer_receiver_patch = fill(iout, NVEGPOOL_NATVEG)
    cf.matrix_gmtransfer_patch = zeros(np, NVEGPOOL_NATVEG)
    cf.matrix_gmturnover_patch = zeros(np, nvegcpool)
    cf.ileaf_to_iout_gm=1;  cf.ileafst_to_iout_gm=2;  cf.ileafxf_to_iout_gm=3
    cf.ifroot_to_iout_gm=4; cf.ifrootst_to_iout_gm=5; cf.ifrootxf_to_iout_gm=6
    cf.ilivestem_to_iout_gm=7;  cf.ilivestemst_to_iout_gm=8;  cf.ilivestemxf_to_iout_gm=9
    cf.ideadstem_to_iout_gm=10; cf.ideadstemst_to_iout_gm=11; cf.ideadstemxf_to_iout_gm=12
    cf.ilivecroot_to_iout_gm=13; cf.ilivecrootst_to_iout_gm=14; cf.ilivecrootxf_to_iout_gm=15
    cf.ideadcroot_to_iout_gm=16; cf.ideadcrootst_to_iout_gm=17; cf.ideadcrootxf_to_iout_gm=18

    # --- Fire: 2 in-veg live->dead conversions (indices 1,2) + 18 pure out-of-veg
    # losses (indices 3..20, pool i -> iout). ncfiouttrans=18 trailing → the 2
    # in-veg transfers MUST be first.
    fi_doner    = vcat([ILIVESTEM, ILIVECROOT], collect(1:NVEGPOOL_NATVEG))
    fi_receiver = vcat([IDEADSTEM, IDEADCROOT], fill(iout, NVEGPOOL_NATVEG))
    cf.matrix_fitransfer_doner_patch    = fi_doner
    cf.matrix_fitransfer_receiver_patch = fi_receiver
    cf.matrix_fitransfer_patch = zeros(np, length(fi_doner))
    cf.matrix_fiturnover_patch = zeros(np, nvegcpool)
    cf.ilivestem_to_ideadstem_fi=1; cf.ilivecroot_to_ideadcroot_fi=2
    cf.ileaf_to_iout_fi=3;  cf.ileafst_to_iout_fi=4;  cf.ileafxf_to_iout_fi=5
    cf.ifroot_to_iout_fi=6; cf.ifrootst_to_iout_fi=7; cf.ifrootxf_to_iout_fi=8
    cf.ilivestem_to_iout_fi=9;  cf.ilivestemst_to_iout_fi=10; cf.ilivestemxf_to_iout_fi=11
    cf.ideadstem_to_iout_fi=12; cf.ideadstemst_to_iout_fi=13; cf.ideadstemxf_to_iout_fi=14
    cf.ilivecroot_to_iout_fi=15; cf.ilivecrootst_to_iout_fi=16; cf.ilivecrootxf_to_iout_fi=17
    cf.ideadcroot_to_iout_fi=18; cf.ideadcrootst_to_iout_fi=19; cf.ideadcrootxf_to_iout_fi=20
    return nothing
end

# The 18 fire out-loss flux fields per pool: (to_fire, to_litter_fire) summed.
const _FIC_TO_FIRE = Symbol[
    :m_leafc_to_fire_patch, :m_leafc_storage_to_fire_patch, :m_leafc_xfer_to_fire_patch,
    :m_frootc_to_fire_patch, :m_frootc_storage_to_fire_patch, :m_frootc_xfer_to_fire_patch,
    :m_livestemc_to_fire_patch, :m_livestemc_storage_to_fire_patch, :m_livestemc_xfer_to_fire_patch,
    :m_deadstemc_to_fire_patch, :m_deadstemc_storage_to_fire_patch, :m_deadstemc_xfer_to_fire_patch,
    :m_livecrootc_to_fire_patch, :m_livecrootc_storage_to_fire_patch, :m_livecrootc_xfer_to_fire_patch,
    :m_deadcrootc_to_fire_patch, :m_deadcrootc_storage_to_fire_patch, :m_deadcrootc_xfer_to_fire_patch]
const _FIC_TO_LITTER = Symbol[
    :m_leafc_to_litter_fire_patch, :m_leafc_storage_to_litter_fire_patch, :m_leafc_xfer_to_litter_fire_patch,
    :m_frootc_to_litter_fire_patch, :m_frootc_storage_to_litter_fire_patch, :m_frootc_xfer_to_litter_fire_patch,
    :m_livestemc_to_litter_fire_patch, :m_livestemc_storage_to_litter_fire_patch, :m_livestemc_xfer_to_litter_fire_patch,
    :m_deadstemc_to_litter_fire_patch, :m_deadstemc_storage_to_litter_fire_patch, :m_deadstemc_xfer_to_litter_fire_patch,
    :m_livecrootc_to_litter_fire_patch, :m_livecrootc_storage_to_litter_fire_patch, :m_livecrootc_xfer_to_litter_fire_patch,
    :m_deadcrootc_to_litter_fire_patch, :m_deadcrootc_storage_to_litter_fire_patch, :m_deadcrootc_xfer_to_litter_fire_patch]

# The 12 allocation targets: (veg pool, cpool_to_* flux field). Allocation fills
# the leaf/froot/live+deadwood pools + their storage (never the xfer pools).
const _ALLOC_TARGET = Tuple{Int,Symbol}[
    (ILEAF, :cpool_to_leafc_patch), (ILEAF_ST, :cpool_to_leafc_storage_patch),
    (IFROOT, :cpool_to_frootc_patch), (IFROOT_ST, :cpool_to_frootc_storage_patch),
    (ILIVESTEM, :cpool_to_livestemc_patch), (ILIVESTEM_ST, :cpool_to_livestemc_storage_patch),
    (IDEADSTEM, :cpool_to_deadstemc_patch), (IDEADSTEM_ST, :cpool_to_deadstemc_storage_patch),
    (ILIVECROOT, :cpool_to_livecrootc_patch), (ILIVECROOT_ST, :cpool_to_livecrootc_storage_patch),
    (IDEADCROOT, :cpool_to_deadcrootc_patch), (IDEADCROOT_ST, :cpool_to_deadcrootc_storage_patch)]

"""
    cn_veg_matrix_alloc_c!(cf, mask_soilp, bounds_patch)

Assemble the veg-C allocation B-input from the cpool_to_* fluxes: the solve adds
`matrix_alloc[p,i] * matrix_Cinput[p] * dt` to pool i, so set matrix_Cinput to the
total allocated C rate and matrix_alloc to the per-pool fraction (the faithful
fraction form, reusable by the isotope B-input). Zero when nothing is allocated.
"""
# The 2 crop grain allocation targets (reproductive C; 2D fields indexed [p,1]).
const _ALLOC_TARGET_CROP = Tuple{Int,Symbol}[
    (IGRAIN, :cpool_to_reproductivec_patch),
    (IGRAIN_ST, :cpool_to_reproductivec_storage_patch)]

# Allocation B-input (12 cpool_to_* targets [+2 crop grain]): matrix_Cinput = total
# allocated rate; matrix_alloc = per-pool fraction. One thread per patch.
Base.@kwdef struct _MatAllocOut{V,M}
    alloc::M; cinput::V
end
Adapt.@adapt_structure _MatAllocOut
Base.@kwdef struct _MatAllocIn{V,M}
    a1::V; a2::V; a3::V; a4::V; a5::V; a6::V; a7::V; a8::V; a9::V; a10::V; a11::V; a12::V
    g1::M; g2::M
end
Adapt.@adapt_structure _MatAllocIn

@kernel function _mat_alloc_c_kernel!(out::_MatAllocOut, in::_MatAllocIn, @Const(mask), use_crop::Bool, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask[p]
        alloc = out.alloc; cin = out.cinput
        tot = in.a1[p] + in.a2[p] + in.a3[p] + in.a4[p] + in.a5[p] + in.a6[p] +
              in.a7[p] + in.a8[p] + in.a9[p] + in.a10[p] + in.a11[p] + in.a12[p]
        if use_crop; tot += in.g1[p, 1] + in.g2[p, 1]; end
        cin[p] = tot
        if tot > 0
            alloc[p, ILEAF] = in.a1[p] / tot;          alloc[p, ILEAF_ST] = in.a2[p] / tot
            alloc[p, IFROOT] = in.a3[p] / tot;         alloc[p, IFROOT_ST] = in.a4[p] / tot
            alloc[p, ILIVESTEM] = in.a5[p] / tot;      alloc[p, ILIVESTEM_ST] = in.a6[p] / tot
            alloc[p, IDEADSTEM] = in.a7[p] / tot;      alloc[p, IDEADSTEM_ST] = in.a8[p] / tot
            alloc[p, ILIVECROOT] = in.a9[p] / tot;     alloc[p, ILIVECROOT_ST] = in.a10[p] / tot
            alloc[p, IDEADCROOT] = in.a11[p] / tot;    alloc[p, IDEADCROOT_ST] = in.a12[p] / tot
            if use_crop
                alloc[p, IGRAIN] = in.g1[p, 1] / tot;  alloc[p, IGRAIN_ST] = in.g2[p, 1] / tot
            end
        end
    end
end

function cn_veg_matrix_alloc_c!(cf::CNVegCarbonFluxData,
        mask_soilp::AbstractVector{Bool}, bounds_patch::UnitRange{Int};
        use_crop::Bool = false)
    fill!(cf.matrix_alloc_patch, zero(eltype(cf.matrix_alloc_patch)))
    fill!(cf.matrix_Cinput_patch, zero(eltype(cf.matrix_Cinput_patch)))
    gf(s) = getfield(cf, s)
    out = _MatAllocOut(; alloc=cf.matrix_alloc_patch, cinput=cf.matrix_Cinput_patch)
    in_ = _MatAllocIn(;
        a1=gf(_ALLOC_TARGET[1][2]), a2=gf(_ALLOC_TARGET[2][2]), a3=gf(_ALLOC_TARGET[3][2]), a4=gf(_ALLOC_TARGET[4][2]),
        a5=gf(_ALLOC_TARGET[5][2]), a6=gf(_ALLOC_TARGET[6][2]), a7=gf(_ALLOC_TARGET[7][2]), a8=gf(_ALLOC_TARGET[8][2]),
        a9=gf(_ALLOC_TARGET[9][2]), a10=gf(_ALLOC_TARGET[10][2]), a11=gf(_ALLOC_TARGET[11][2]), a12=gf(_ALLOC_TARGET[12][2]),
        g1=gf(_ALLOC_TARGET_CROP[1][2]), g2=gf(_ALLOC_TARGET_CROP[2][2]))
    backend = _kernel_backend(cf.matrix_alloc_patch)
    _mat_alloc_c_kernel!(backend)(out, in_, mask_soilp, use_crop, first(bounds_patch), last(bounds_patch); ndrange = last(bounds_patch))
    KA.synchronize(backend)
    return nothing
end

"""
    cn_veg_matrix_accumulate_fi_c!(cf, cs, mask_soilp, bounds_patch; dt, matrixcheck=false)

Register the veg-C fire transfers: 2 in-veg live->dead conversions
(m_livestemc_to_deadstemc_fire / m_livecrootc_to_deadcrootc_fire) + 18 out-of-veg
losses (each pool -> iout, flux = m_*_to_fire + m_*_to_litter_fire), via
matrix_update_fic!.
"""
# One fire register: transfer idx, donor column d. Cap = ph + gm + fi turnover for pool d.
@inline function _mat_fi_reg!(fit, fitn, phtn, gmtn, p::Int, idx::Int, d::Int, flux, donorval, dt, mc::Bool)
    rate = donorval > 0 ? flux / donorval : zero(flux)
    tot = phtn[p, d] + gmtn[p, d] + fitn[p, d]
    if mc && (tot + rate * dt >= one(rate))
        applied = max(zero(rate), (one(rate) - tot) / dt)
    else
        applied = rate
    end
    @inbounds fitn[p, d]   = fitn[p, d]   + applied * dt
    @inbounds fit[p, idx]  = fit[p, idx]  + applied
    return nothing
end

Base.@kwdef struct _MatFiOut{M}
    fit::M; fitn::M; phtn::M; gmtn::M   # phtn/gmtn read-only (for the cap)
end
Adapt.@adapt_structure _MatFiOut
Base.@kwdef struct _MatFiIn{V}
    ff1::V; ff2::V; ff3::V; ff4::V; ff5::V; ff6::V; ff7::V; ff8::V; ff9::V
    ff10::V; ff11::V; ff12::V; ff13::V; ff14::V; ff15::V; ff16::V; ff17::V; ff18::V
    fl1::V; fl2::V; fl3::V; fl4::V; fl5::V; fl6::V; fl7::V; fl8::V; fl9::V
    fl10::V; fl11::V; fl12::V; fl13::V; fl14::V; fl15::V; fl16::V; fl17::V; fl18::V
    po1::V; po2::V; po3::V; po4::V; po5::V; po6::V; po7::V; po8::V; po9::V
    po10::V; po11::V; po12::V; po13::V; po14::V; po15::V; po16::V; po17::V; po18::V
    f_ls_dead::V; f_lc_dead::V
end
Adapt.@adapt_structure _MatFiIn
Base.@kwdef struct _MatFiIdx
    i_ls::Int; i_lc::Int; matrixcheck::Bool
end

@kernel function _mat_fi_c_kernel!(out::_MatFiOut, in::_MatFiIn, @Const(mask), ix::_MatFiIdx, dt, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask[p]
        fit = out.fit; fitn = out.fitn; phtn = out.phtn; gmtn = out.gmtn; mc = ix.matrixcheck
        # live->dead conversions FIRST (their turnover contributes to the pool-7/13 cap below)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, ix.i_ls, ILIVESTEM,  in.f_ls_dead[p], in.po7[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, ix.i_lc, ILIVECROOT, in.f_lc_dead[p], in.po13[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 3,  1,  in.ff1[p]  + in.fl1[p],  in.po1[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 4,  2,  in.ff2[p]  + in.fl2[p],  in.po2[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 5,  3,  in.ff3[p]  + in.fl3[p],  in.po3[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 6,  4,  in.ff4[p]  + in.fl4[p],  in.po4[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 7,  5,  in.ff5[p]  + in.fl5[p],  in.po5[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 8,  6,  in.ff6[p]  + in.fl6[p],  in.po6[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 9,  7,  in.ff7[p]  + in.fl7[p],  in.po7[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 10, 8,  in.ff8[p]  + in.fl8[p],  in.po8[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 11, 9,  in.ff9[p]  + in.fl9[p],  in.po9[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 12, 10, in.ff10[p] + in.fl10[p], in.po10[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 13, 11, in.ff11[p] + in.fl11[p], in.po11[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 14, 12, in.ff12[p] + in.fl12[p], in.po12[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 15, 13, in.ff13[p] + in.fl13[p], in.po13[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 16, 14, in.ff14[p] + in.fl14[p], in.po14[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 17, 15, in.ff15[p] + in.fl15[p], in.po15[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 18, 16, in.ff16[p] + in.fl16[p], in.po16[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 19, 17, in.ff17[p] + in.fl17[p], in.po17[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 20, 18, in.ff18[p] + in.fl18[p], in.po18[p], dt, mc)
    end
end

function cn_veg_matrix_accumulate_fi_c!(cf::CNVegCarbonFluxData,
        cs::CNVegCarbonStateData, mask_soilp::AbstractVector{Bool},
        bounds_patch::UnitRange{Int}; dt::Real, matrixcheck::Bool = false)
    fill!(cf.matrix_fitransfer_patch, zero(eltype(cf.matrix_fitransfer_patch)))
    fill!(cf.matrix_fiturnover_patch, zero(eltype(cf.matrix_fiturnover_patch)))
    gf(sym) = getfield(cf, sym)
    out = _MatFiOut(; fit=cf.matrix_fitransfer_patch, fitn=cf.matrix_fiturnover_patch,
        phtn=cf.matrix_phturnover_patch, gmtn=cf.matrix_gmturnover_patch)
    in_ = _MatFiIn(;
        ff1=gf(_FIC_TO_FIRE[1]), ff2=gf(_FIC_TO_FIRE[2]), ff3=gf(_FIC_TO_FIRE[3]), ff4=gf(_FIC_TO_FIRE[4]),
        ff5=gf(_FIC_TO_FIRE[5]), ff6=gf(_FIC_TO_FIRE[6]), ff7=gf(_FIC_TO_FIRE[7]), ff8=gf(_FIC_TO_FIRE[8]),
        ff9=gf(_FIC_TO_FIRE[9]), ff10=gf(_FIC_TO_FIRE[10]), ff11=gf(_FIC_TO_FIRE[11]), ff12=gf(_FIC_TO_FIRE[12]),
        ff13=gf(_FIC_TO_FIRE[13]), ff14=gf(_FIC_TO_FIRE[14]), ff15=gf(_FIC_TO_FIRE[15]), ff16=gf(_FIC_TO_FIRE[16]),
        ff17=gf(_FIC_TO_FIRE[17]), ff18=gf(_FIC_TO_FIRE[18]),
        fl1=gf(_FIC_TO_LITTER[1]), fl2=gf(_FIC_TO_LITTER[2]), fl3=gf(_FIC_TO_LITTER[3]), fl4=gf(_FIC_TO_LITTER[4]),
        fl5=gf(_FIC_TO_LITTER[5]), fl6=gf(_FIC_TO_LITTER[6]), fl7=gf(_FIC_TO_LITTER[7]), fl8=gf(_FIC_TO_LITTER[8]),
        fl9=gf(_FIC_TO_LITTER[9]), fl10=gf(_FIC_TO_LITTER[10]), fl11=gf(_FIC_TO_LITTER[11]), fl12=gf(_FIC_TO_LITTER[12]),
        fl13=gf(_FIC_TO_LITTER[13]), fl14=gf(_FIC_TO_LITTER[14]), fl15=gf(_FIC_TO_LITTER[15]), fl16=gf(_FIC_TO_LITTER[16]),
        fl17=gf(_FIC_TO_LITTER[17]), fl18=gf(_FIC_TO_LITTER[18]),
        po1=cs.leafc_patch, po2=cs.leafc_storage_patch, po3=cs.leafc_xfer_patch,
        po4=cs.frootc_patch, po5=cs.frootc_storage_patch, po6=cs.frootc_xfer_patch,
        po7=cs.livestemc_patch, po8=cs.livestemc_storage_patch, po9=cs.livestemc_xfer_patch,
        po10=cs.deadstemc_patch, po11=cs.deadstemc_storage_patch, po12=cs.deadstemc_xfer_patch,
        po13=cs.livecrootc_patch, po14=cs.livecrootc_storage_patch, po15=cs.livecrootc_xfer_patch,
        po16=cs.deadcrootc_patch, po17=cs.deadcrootc_storage_patch, po18=cs.deadcrootc_xfer_patch,
        f_ls_dead=cf.m_livestemc_to_deadstemc_fire_patch, f_lc_dead=cf.m_livecrootc_to_deadcrootc_fire_patch)
    ix = _MatFiIdx(; i_ls=cf.ilivestem_to_ideadstem_fi, i_lc=cf.ilivecroot_to_ideadcroot_fi, matrixcheck=matrixcheck)
    backend = _kernel_backend(cf.matrix_fitransfer_patch)
    _mat_fi_c_kernel!(backend)(out, in_, mask_soilp, ix,
        eltype(cf.matrix_fitransfer_patch)(dt), first(bounds_patch), last(bounds_patch); ndrange = last(bounds_patch))
    KA.synchronize(backend)
    return nothing
end

# The 18 gap-mortality m_*_to_litter flux fields, in pool-index order (donor = i).
const _GMC_FLUX = Symbol[
    :m_leafc_to_litter_patch, :m_leafc_storage_to_litter_patch, :m_leafc_xfer_to_litter_patch,
    :m_frootc_to_litter_patch, :m_frootc_storage_to_litter_patch, :m_frootc_xfer_to_litter_patch,
    :m_livestemc_to_litter_patch, :m_livestemc_storage_to_litter_patch, :m_livestemc_xfer_to_litter_patch,
    :m_deadstemc_to_litter_patch, :m_deadstemc_storage_to_litter_patch, :m_deadstemc_xfer_to_litter_patch,
    :m_livecrootc_to_litter_patch, :m_livecrootc_storage_to_litter_patch, :m_livecrootc_xfer_to_litter_patch,
    :m_deadcrootc_to_litter_patch, :m_deadcrootc_storage_to_litter_patch, :m_deadcrootc_xfer_to_litter_patch]

# One gap-mortality register (transfer i = pool i → litter; doner[i]=i). Cap accounts for
# the already-accumulated phenology + gap turnover for the same pool.
@inline function _mat_gm_reg!(gmt, gmtn, phtn, p::Int, i::Int, flux, donorval, dt, mc::Bool)
    rate = donorval > 0 ? flux / donorval : zero(flux)
    if mc && (phtn[p, i] + gmtn[p, i] + rate * dt >= one(rate))
        applied = max(zero(rate), (one(rate) - phtn[p, i] - gmtn[p, i]) / dt)
    else
        applied = rate
    end
    @inbounds gmtn[p, i] = gmtn[p, i] + applied * dt
    @inbounds gmt[p, i]  = gmt[p, i]  + applied
    return nothing
end

Base.@kwdef struct _MatGmOut{M}
    gmt::M; gmtn::M; phtn::M   # phtn read-only (prior phenology turnover, for the cap)
end
Adapt.@adapt_structure _MatGmOut
# 18 gap-mortality litter fluxes (pool order) + 18 donor-pool values.
Base.@kwdef struct _MatGmIn{V}
    fl1::V; fl2::V; fl3::V; fl4::V; fl5::V; fl6::V; fl7::V; fl8::V; fl9::V
    fl10::V; fl11::V; fl12::V; fl13::V; fl14::V; fl15::V; fl16::V; fl17::V; fl18::V
    po1::V; po2::V; po3::V; po4::V; po5::V; po6::V; po7::V; po8::V; po9::V
    po10::V; po11::V; po12::V; po13::V; po14::V; po15::V; po16::V; po17::V; po18::V
end
Adapt.@adapt_structure _MatGmIn

@kernel function _mat_gm_c_kernel!(out::_MatGmOut, in::_MatGmIn, @Const(mask), mc::Bool, dt, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask[p]
        gmt = out.gmt; gmtn = out.gmtn; phtn = out.phtn
        _mat_gm_reg!(gmt, gmtn, phtn, p, 1,  in.fl1[p],  in.po1[p],  dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 2,  in.fl2[p],  in.po2[p],  dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 3,  in.fl3[p],  in.po3[p],  dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 4,  in.fl4[p],  in.po4[p],  dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 5,  in.fl5[p],  in.po5[p],  dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 6,  in.fl6[p],  in.po6[p],  dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 7,  in.fl7[p],  in.po7[p],  dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 8,  in.fl8[p],  in.po8[p],  dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 9,  in.fl9[p],  in.po9[p],  dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 10, in.fl10[p], in.po10[p], dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 11, in.fl11[p], in.po11[p], dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 12, in.fl12[p], in.po12[p], dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 13, in.fl13[p], in.po13[p], dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 14, in.fl14[p], in.po14[p], dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 15, in.fl15[p], in.po15[p], dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 16, in.fl16[p], in.po16[p], dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 17, in.fl17[p], in.po17[p], dt, mc)
        _mat_gm_reg!(gmt, gmtn, phtn, p, 18, in.fl18[p], in.po18[p], dt, mc)
    end
end

"""
    cn_veg_matrix_accumulate_gm_c!(cf, cs, mask_soilp, bounds_patch; dt, matrixcheck=false)

Register the 18 veg-C gap-mortality transfers (each veg pool -> litter, a pure
out-of-veg loss) into the matrix from the m_*_to_litter flux fields gap_mortality
already computed. rate = flux/pool; recorded via matrix_update_gmc!.
"""
function cn_veg_matrix_accumulate_gm_c!(cf::CNVegCarbonFluxData,
        cs::CNVegCarbonStateData, mask_soilp::AbstractVector{Bool},
        bounds_patch::UnitRange{Int}; dt::Real, matrixcheck::Bool = false)
    fill!(cf.matrix_gmtransfer_patch, zero(eltype(cf.matrix_gmtransfer_patch)))
    fill!(cf.matrix_gmturnover_patch, zero(eltype(cf.matrix_gmturnover_patch)))
    out = _MatGmOut(; gmt=cf.matrix_gmtransfer_patch, gmtn=cf.matrix_gmturnover_patch, phtn=cf.matrix_phturnover_patch)
    in_ = _MatGmIn(;
        fl1=getfield(cf, _GMC_FLUX[1]), fl2=getfield(cf, _GMC_FLUX[2]), fl3=getfield(cf, _GMC_FLUX[3]),
        fl4=getfield(cf, _GMC_FLUX[4]), fl5=getfield(cf, _GMC_FLUX[5]), fl6=getfield(cf, _GMC_FLUX[6]),
        fl7=getfield(cf, _GMC_FLUX[7]), fl8=getfield(cf, _GMC_FLUX[8]), fl9=getfield(cf, _GMC_FLUX[9]),
        fl10=getfield(cf, _GMC_FLUX[10]), fl11=getfield(cf, _GMC_FLUX[11]), fl12=getfield(cf, _GMC_FLUX[12]),
        fl13=getfield(cf, _GMC_FLUX[13]), fl14=getfield(cf, _GMC_FLUX[14]), fl15=getfield(cf, _GMC_FLUX[15]),
        fl16=getfield(cf, _GMC_FLUX[16]), fl17=getfield(cf, _GMC_FLUX[17]), fl18=getfield(cf, _GMC_FLUX[18]),
        po1=cs.leafc_patch, po2=cs.leafc_storage_patch, po3=cs.leafc_xfer_patch,
        po4=cs.frootc_patch, po5=cs.frootc_storage_patch, po6=cs.frootc_xfer_patch,
        po7=cs.livestemc_patch, po8=cs.livestemc_storage_patch, po9=cs.livestemc_xfer_patch,
        po10=cs.deadstemc_patch, po11=cs.deadstemc_storage_patch, po12=cs.deadstemc_xfer_patch,
        po13=cs.livecrootc_patch, po14=cs.livecrootc_storage_patch, po15=cs.livecrootc_xfer_patch,
        po16=cs.deadcrootc_patch, po17=cs.deadcrootc_storage_patch, po18=cs.deadcrootc_xfer_patch)
    backend = _kernel_backend(cf.matrix_gmtransfer_patch)
    _mat_gm_c_kernel!(backend)(out, in_, mask_soilp, matrixcheck,
        eltype(cf.matrix_gmtransfer_patch)(dt), first(bounds_patch), last(bounds_patch); ndrange = last(bounds_patch))
    KA.synchronize(backend)
    return nothing
end

# ==========================================================================
# GPU: per-patch kernelization of the veg-C phenology accumulate. The reg() calls
# are a per-patch sequential accumulation into matrix_phtransfer/phturnover — one
# KA thread per patch (the register order is preserved in-thread, so the shared
# donor turnover accumulation + the matrixcheck cap stay byte-identical). The ~17
# flux inputs + ~16 donor-pool inputs are bundled into ONE device-view struct arg
# (collapses to a single Metal buffer), mirroring the BGC-module kernelization.
# ==========================================================================

# One transfer register: rate = flux/donor (0 if donor<=0), then (matrixcheck cap)
# turnover[p,d] += applied·dt ; transfer[p,idx] += applied. `d` is the donor pool
# column (== matrix_phtransfer_doner_patch[idx] by topology construction).
@inline function _mat_ph_reg!(pht, phtn, p::Int, idx::Int, flux, d::Int, donorval, dt, matrixcheck::Bool)
    idx == 0 && return nothing
    rate = donorval > 0 ? flux / donorval : zero(flux)
    if matrixcheck && (phtn[p, d] + rate * dt >= one(rate))
        applied = max(zero(rate), (one(rate) - phtn[p, d]) / dt)
    else
        applied = rate
    end
    @inbounds phtn[p, d]   = phtn[p, d]   + applied * dt
    @inbounds pht[p, idx]  = pht[p, idx]  + applied
    return nothing
end

# Output bundle (2D accumulators) + input bundle (1D flux + donor-pool fields;
# crop grain fields are 2D → separate param). One kernel arg each.
Base.@kwdef struct _MatPhCOut{M}
    pht::M; phtn::M
end
Adapt.@adapt_structure _MatPhCOut
Base.@kwdef struct _MatPhCIn{V,M}
    f_leafc_st::V; f_leafc_xf::V; f_frootc_st::V; f_frootc_xf::V
    f_livestemc_st::V; f_livestemc_xf::V; f_deadstemc_st::V; f_deadstemc_xf::V
    f_livecrootc_st::V; f_livecrootc_xf::V; f_deadcrootc_st::V; f_deadcrootc_xf::V
    f_livestem_dead::V; f_livecroot_dead::V
    f_leaf_lit::V; f_froot_lit::V; f_livestem_lit::V
    p_leafc_st::V; p_leafc_xf::V; p_frootc_st::V; p_frootc_xf::V
    p_livestemc_st::V; p_livestemc_xf::V; p_deadstemc_st::V; p_deadstemc_xf::V
    p_livecrootc_st::V; p_livecrootc_xf::V; p_deadcrootc_st::V; p_deadcrootc_xf::V
    p_livestemc::V; p_livecrootc::V; p_leafc::V; p_frootc::V
    grain_food::M; grain_seed::M; p_repro::M
end
Adapt.@adapt_structure _MatPhCIn

# Scalar bundle: the 18 topology entry indices + flags. isbits → materializes on Metal.
Base.@kwdef struct _MatPhCIdx
    i1::Int; i2::Int; i3::Int; i4::Int; i5::Int; i6::Int; i7::Int; i8::Int; i9::Int
    i10::Int; i11::Int; i12::Int; i13::Int; i14::Int; i15::Int; i16::Int; i17::Int; i18::Int
    use_crop::Bool; matrixcheck::Bool
end

@kernel function _mat_ph_c_kernel!(out::_MatPhCOut, in::_MatPhCIn, @Const(mask), ix::_MatPhCIdx, dt, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask[p]
        pht = out.pht; phtn = out.phtn; mc = ix.matrixcheck
        _mat_ph_reg!(pht, phtn, p, ix.i1,  in.f_leafc_st[p],      ILEAF_ST,      in.p_leafc_st[p],      dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i2,  in.f_leafc_xf[p],      ILEAF_XF,      in.p_leafc_xf[p],      dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i3,  in.f_frootc_st[p],     IFROOT_ST,     in.p_frootc_st[p],     dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i4,  in.f_frootc_xf[p],     IFROOT_XF,     in.p_frootc_xf[p],     dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i5,  in.f_livestemc_st[p],  ILIVESTEM_ST,  in.p_livestemc_st[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i6,  in.f_livestemc_xf[p],  ILIVESTEM_XF,  in.p_livestemc_xf[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i7,  in.f_deadstemc_st[p],  IDEADSTEM_ST,  in.p_deadstemc_st[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i8,  in.f_deadstemc_xf[p],  IDEADSTEM_XF,  in.p_deadstemc_xf[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i9,  in.f_livecrootc_st[p], ILIVECROOT_ST, in.p_livecrootc_st[p], dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i10, in.f_livecrootc_xf[p], ILIVECROOT_XF, in.p_livecrootc_xf[p], dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i11, in.f_deadcrootc_st[p], IDEADCROOT_ST, in.p_deadcrootc_st[p], dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i12, in.f_deadcrootc_xf[p], IDEADCROOT_XF, in.p_deadcrootc_xf[p], dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i13, in.f_livestem_dead[p], ILIVESTEM,     in.p_livestemc[p],     dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i14, in.f_livecroot_dead[p], ILIVECROOT,   in.p_livecrootc[p],    dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i15, in.f_leaf_lit[p],      ILEAF,         in.p_leafc[p],         dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i16, in.f_froot_lit[p],     IFROOT,        in.p_frootc[p],        dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.i17, in.f_livestem_lit[p],  ILIVESTEM,     in.p_livestemc[p],     dt, mc)
        if ix.use_crop
            _mat_ph_reg!(pht, phtn, p, ix.i18, in.grain_food[p, 1] + in.grain_seed[p, 1], IGRAIN, in.p_repro[p, 1], dt, mc)
        end
    end
end

"""
    cn_veg_matrix_accumulate_ph_c!(cf, cs, mask_soilp, bounds_patch;
                                    dt, matrixcheck=false)

Register every veg-C phenology transfer into the matrix from the per-flux
transfer fields that phenology already computed: for each transfer, the rate is
recovered as `flux / donor_pool` (0 when the donor pool is empty, in which case
the flux is 0 too) and recorded via `matrix_update_phc!`. Zeros the
transfer/turnover accumulators first (they accumulate within the step).

Requires `cn_veg_matrix_c_topology!(cf)` to have been called (indices set).
"""
function cn_veg_matrix_accumulate_ph_c!(cf::CNVegCarbonFluxData,
        cs::CNVegCarbonStateData, mask_soilp::AbstractVector{Bool},
        bounds_patch::UnitRange{Int}; dt::Real, matrixcheck::Bool = false,
        use_crop::Bool = false)

    fill!(cf.matrix_phtransfer_patch, zero(eltype(cf.matrix_phtransfer_patch)))
    fill!(cf.matrix_phturnover_patch, zero(eltype(cf.matrix_phturnover_patch)))

    out = _MatPhCOut(; pht=cf.matrix_phtransfer_patch, phtn=cf.matrix_phturnover_patch)
    in_ = _MatPhCIn(;
        f_leafc_st=cf.leafc_storage_to_xfer_patch, f_leafc_xf=cf.leafc_xfer_to_leafc_patch,
        f_frootc_st=cf.frootc_storage_to_xfer_patch, f_frootc_xf=cf.frootc_xfer_to_frootc_patch,
        f_livestemc_st=cf.livestemc_storage_to_xfer_patch, f_livestemc_xf=cf.livestemc_xfer_to_livestemc_patch,
        f_deadstemc_st=cf.deadstemc_storage_to_xfer_patch, f_deadstemc_xf=cf.deadstemc_xfer_to_deadstemc_patch,
        f_livecrootc_st=cf.livecrootc_storage_to_xfer_patch, f_livecrootc_xf=cf.livecrootc_xfer_to_livecrootc_patch,
        f_deadcrootc_st=cf.deadcrootc_storage_to_xfer_patch, f_deadcrootc_xf=cf.deadcrootc_xfer_to_deadcrootc_patch,
        f_livestem_dead=cf.livestemc_to_deadstemc_patch, f_livecroot_dead=cf.livecrootc_to_deadcrootc_patch,
        f_leaf_lit=cf.leafc_to_litter_patch, f_froot_lit=cf.frootc_to_litter_patch, f_livestem_lit=cf.livestemc_to_litter_patch,
        p_leafc_st=cs.leafc_storage_patch, p_leafc_xf=cs.leafc_xfer_patch,
        p_frootc_st=cs.frootc_storage_patch, p_frootc_xf=cs.frootc_xfer_patch,
        p_livestemc_st=cs.livestemc_storage_patch, p_livestemc_xf=cs.livestemc_xfer_patch,
        p_deadstemc_st=cs.deadstemc_storage_patch, p_deadstemc_xf=cs.deadstemc_xfer_patch,
        p_livecrootc_st=cs.livecrootc_storage_patch, p_livecrootc_xf=cs.livecrootc_xfer_patch,
        p_deadcrootc_st=cs.deadcrootc_storage_patch, p_deadcrootc_xf=cs.deadcrootc_xfer_patch,
        p_livestemc=cs.livestemc_patch, p_livecrootc=cs.livecrootc_patch,
        p_leafc=cs.leafc_patch, p_frootc=cs.frootc_patch,
        grain_food=cf.repr_grainc_to_food_patch, grain_seed=cf.repr_grainc_to_seed_patch, p_repro=cs.reproductivec_patch)
    ix = _MatPhCIdx(;
        i1=cf.ileafst_to_ileafxf_ph, i2=cf.ileafxf_to_ileaf_ph, i3=cf.ifrootst_to_ifrootxf_ph,
        i4=cf.ifrootxf_to_ifroot_ph, i5=cf.ilivestemst_to_ilivestemxf_ph, i6=cf.ilivestemxf_to_ilivestem_ph,
        i7=cf.ideadstemst_to_ideadstemxf_ph, i8=cf.ideadstemxf_to_ideadstem_ph,
        i9=cf.ilivecrootst_to_ilivecrootxf_ph, i10=cf.ilivecrootxf_to_ilivecroot_ph,
        i11=cf.ideadcrootst_to_ideadcrootxf_ph, i12=cf.ideadcrootxf_to_ideadcroot_ph,
        i13=cf.ilivestem_to_ideadstem_ph, i14=cf.ilivecroot_to_ideadcroot_ph,
        i15=cf.ileaf_to_iout_ph, i16=cf.ifroot_to_iout_ph, i17=cf.ilivestem_to_iout_ph,
        i18=cf.igrain_to_iout_ph, use_crop=use_crop, matrixcheck=matrixcheck)
    backend = _kernel_backend(cf.matrix_phtransfer_patch)
    _mat_ph_c_kernel!(backend)(out, in_, mask_soilp, ix, eltype(cf.matrix_phtransfer_patch)(dt),
        first(bounds_patch), last(bounds_patch); ndrange = last(bounds_patch))
    KA.synchronize(backend)
    return nothing
end

# ==========================================================================
# Matrix-CN veg-N wiring (phase 2) — same separate-accumulation architecture as
# phase 1, with the retranslocation pool (iretransn = nvegnpool = 19 for natveg;
# ioutn = 20). Pools 1..18 mirror the C pools (ILEAF..IDEADCROOT_XF); pool 19 =
# retransn (ns.retransn_patch). The N solve expects the FORTRAN convention where
# nvegnpool INCLUDES retransn, so wiring passes nvegnpool = vp.iretransn (=19).
#
# Retranslocation accounting (NutrientCompetitionFlexibleCNMod.F90:880–922):
#   Ninput = npool_to_veg - retransn_to_npool          (sminn-sourced B-input)
#   nalloc[X] = npool_to_Xn / npool_to_veg             (allocation fraction)
#   retransn→X transfer rate = nalloc[X] * retransn_to_npool / retransn
# so pool X receives nalloc[X]*retransn_to_npool*dt (A-matrix) + nalloc[X]*Ninput*dt
# (B-input) = npool_to_Xn*dt, and retransn depletes by Σ_X = retransn_to_npool*dt
# (Σnalloc=1) — no double count. iretransn→iout (index 34) fires only in the
# npool_to_veg==0 edge case (rate retransn_to_npool/retransn).
# ==========================================================================

const IRETRANSN_NATVEG = NVEGPOOL_NATVEG + 1   # = 19 (Fortran iretransn, natveg)

# Pool index -> CNVegNitrogenState value at patch p. `iretransn` is the retransn pool
# index (19 natveg, 22 crop); grain (19..21) only exists for crop and is checked after
# retransn so the natveg overlap (iretransn==19==IGRAIN) resolves to retransn.
@inline function _vegn_pool_val(ns::CNVegNitrogenStateData, i::Int, p::Int,
                                iretransn::Int = IRETRANSN_NATVEG)
    @inbounds if i == iretransn          ; return ns.retransn_patch[p]
    elseif i == IGRAIN         ; return ns.reproductiven_patch[p, 1]
    elseif i == IGRAIN_ST      ; return ns.reproductiven_storage_patch[p, 1]
    elseif i == IGRAIN_XF      ; return ns.reproductiven_xfer_patch[p, 1]
    elseif i == ILEAF          ; return ns.leafn_patch[p]
    elseif i == ILEAF_ST       ; return ns.leafn_storage_patch[p]
    elseif i == ILEAF_XF       ; return ns.leafn_xfer_patch[p]
    elseif i == IFROOT         ; return ns.frootn_patch[p]
    elseif i == IFROOT_ST      ; return ns.frootn_storage_patch[p]
    elseif i == IFROOT_XF      ; return ns.frootn_xfer_patch[p]
    elseif i == ILIVESTEM      ; return ns.livestemn_patch[p]
    elseif i == ILIVESTEM_ST   ; return ns.livestemn_storage_patch[p]
    elseif i == ILIVESTEM_XF   ; return ns.livestemn_xfer_patch[p]
    elseif i == IDEADSTEM      ; return ns.deadstemn_patch[p]
    elseif i == IDEADSTEM_ST   ; return ns.deadstemn_storage_patch[p]
    elseif i == IDEADSTEM_XF   ; return ns.deadstemn_xfer_patch[p]
    elseif i == ILIVECROOT     ; return ns.livecrootn_patch[p]
    elseif i == ILIVECROOT_ST  ; return ns.livecrootn_storage_patch[p]
    elseif i == ILIVECROOT_XF  ; return ns.livecrootn_xfer_patch[p]
    elseif i == IDEADCROOT     ; return ns.deadcrootn_patch[p]
    elseif i == IDEADCROOT_ST  ; return ns.deadcrootn_storage_patch[p]
    elseif i == IDEADCROOT_XF  ; return ns.deadcrootn_xfer_patch[p]
    else                       ; return 0.0
    end
end

# 34 natveg phenology transfers, EXACT Fortran order (CNVegNitrogenFluxType.F90
# :450–585). Out-transfers (31–34) are the trailing nnphouttrans=4.
const _PHN_DONER_NATVEG = Int[
    ILEAF, ILEAF_ST, ILEAF_XF, IFROOT, IFROOT_ST, IFROOT_XF,
    ILIVESTEM, ILIVESTEM, ILIVESTEM_ST, ILIVESTEM_XF, IDEADSTEM_ST, IDEADSTEM_XF,
    ILIVECROOT, ILIVECROOT, ILIVECROOT_ST, ILIVECROOT_XF, IDEADCROOT_ST, IDEADCROOT_XF,
    IRETRANSN_NATVEG, IRETRANSN_NATVEG, IRETRANSN_NATVEG, IRETRANSN_NATVEG,
    IRETRANSN_NATVEG, IRETRANSN_NATVEG, IRETRANSN_NATVEG, IRETRANSN_NATVEG,
    IRETRANSN_NATVEG, IRETRANSN_NATVEG, IRETRANSN_NATVEG, IRETRANSN_NATVEG,
    ILEAF, IFROOT, ILIVESTEM, IRETRANSN_NATVEG]
const _PHN_RECEIVER_NATVEG = Int[
    IRETRANSN_NATVEG, ILEAF_XF, ILEAF, IRETRANSN_NATVEG, IFROOT_XF, IFROOT,
    IDEADSTEM, IRETRANSN_NATVEG, ILIVESTEM_XF, ILIVESTEM, IDEADSTEM_XF, IDEADSTEM,
    IDEADCROOT, IRETRANSN_NATVEG, ILIVECROOT_XF, ILIVECROOT, IDEADCROOT_XF, IDEADCROOT,
    ILEAF, ILEAF_ST, IFROOT, IFROOT_ST, ILIVESTEM, ILIVESTEM_ST,
    IDEADSTEM, IDEADSTEM_ST, ILIVECROOT, ILIVECROOT_ST, IDEADCROOT, IDEADCROOT_ST,
    0, 0, 0, 0]   # 31-34 receivers = ioutn (filled in topology)

"""
    cn_veg_matrix_n_topology!(nf; use_crop=false, nvegnpool=IRETRANSN_NATVEG)

Set up the veg-N transfer topology (ph 34 / gm 19 / fi 21 for natveg): assign the
named `_ph`/`_gm`/`_fi` index fields, populate the doner/receiver arrays, and size
the accumulator matrices. `nvegnpool` INCLUDES retransn (=19 natveg). Idempotent.
"""
function cn_veg_matrix_n_topology!(nf::CNVegNitrogenFluxData; use_crop::Bool = false,
                                   nvegnpool::Int = IRETRANSN_NATVEG)
    iretransn = nvegnpool               # 19 natveg, 22 crop (grain occupies 19..21)
    ioutn = nvegnpool + 1
    # in-veg transfers 1..30 (natveg structure); remap the baked-in retransn pool
    # (IRETRANSN_NATVEG) to the actual iretransn for crop (grain shifts it to 22).
    doner    = _PHN_DONER_NATVEG[1:30]
    receiver = _PHN_RECEIVER_NATVEG[1:30]
    if iretransn != IRETRANSN_NATVEG
        for k in 1:30
            doner[k]    == IRETRANSN_NATVEG && (doner[k]    = iretransn)
            receiver[k] == IRETRANSN_NATVEG && (receiver[k] = iretransn)
        end
    end
    if use_crop
        # +2 retransn→grain supply (31,32), then 5 out-transfers (33..37).
        append!(doner,    [iretransn, iretransn, ILEAF, IFROOT, ILIVESTEM, IGRAIN, iretransn])
        append!(receiver, [IGRAIN, IGRAIN_ST, ioutn, ioutn, ioutn, ioutn, ioutn])
    else
        # 4 out-transfers (31..34).
        append!(doner,    [ILEAF, IFROOT, ILIVESTEM, iretransn])
        append!(receiver, [ioutn, ioutn, ioutn, ioutn])
    end
    nf.matrix_nphtransfer_doner_patch    = doner
    nf.matrix_nphtransfer_receiver_patch = receiver

    np = length(nf.leafn_to_litter_patch)
    nf.matrix_nphtransfer_patch = zeros(np, length(doner))
    nf.matrix_nphturnover_patch = zeros(np, nvegnpool)
    nf.matrix_nalloc_patch  = zeros(np, nvegnpool)
    nf.matrix_Ninput_patch  = zeros(np)

    nf.ileaf_to_iretransn_ph            = 1
    nf.ileafst_to_ileafxf_ph            = 2
    nf.ileafxf_to_ileaf_ph              = 3
    nf.ifroot_to_iretransn_ph           = 4
    nf.ifrootst_to_ifrootxf_ph          = 5
    nf.ifrootxf_to_ifroot_ph            = 6
    nf.ilivestem_to_ideadstem_ph        = 7
    nf.ilivestem_to_iretransn_ph        = 8
    nf.ilivestemst_to_ilivestemxf_ph    = 9
    nf.ilivestemxf_to_ilivestem_ph      = 10
    nf.ideadstemst_to_ideadstemxf_ph    = 11
    nf.ideadstemxf_to_ideadstem_ph      = 12
    nf.ilivecroot_to_ideadcroot_ph      = 13
    nf.ilivecroot_to_iretransn_ph       = 14
    nf.ilivecrootst_to_ilivecrootxf_ph  = 15
    nf.ilivecrootxf_to_ilivecroot_ph    = 16
    nf.ideadcrootst_to_ideadcrootxf_ph  = 17
    nf.ideadcrootxf_to_ideadcroot_ph    = 18
    nf.iretransn_to_ileaf_ph            = 19
    nf.iretransn_to_ileafst_ph          = 20
    nf.iretransn_to_ifroot_ph           = 21
    nf.iretransn_to_ifrootst_ph         = 22
    nf.iretransn_to_ilivestem_ph        = 23
    nf.iretransn_to_ilivestemst_ph      = 24
    nf.iretransn_to_ideadstem_ph        = 25
    nf.iretransn_to_ideadstemst_ph      = 26
    nf.iretransn_to_ilivecroot_ph       = 27
    nf.iretransn_to_ilivecrootst_ph     = 28
    nf.iretransn_to_ideadcroot_ph       = 29
    nf.iretransn_to_ideadcrootst_ph     = 30
    if use_crop
        nf.iretransn_to_igrain_ph       = 31
        nf.iretransn_to_igrainst_ph     = 32
        nf.ileaf_to_iout_ph             = 33
        nf.ifroot_to_iout_ph            = 34
        nf.ilivestem_to_iout_ph         = 35
        nf.igrain_to_iout_ph            = 36
        nf.iretransn_to_iout_ph         = 37
    else
        nf.ileaf_to_iout_ph             = 31
        nf.ifroot_to_iout_ph            = 32
        nf.ilivestem_to_iout_ph         = 33
        nf.iretransn_to_iout_ph         = 34
    end

    # --- Gap mortality: 19 pure out-of-veg losses (18 natveg pools + retransn -> iout).
    # Grain (19..21 for crop) does NOT gap/burn, so the doner set is the 18 natveg pools
    # plus retransn (=iretransn), never the grain pools even when nvegnpool=22.
    ngm = NVEGPOOL_NATVEG + 1
    nf.matrix_ngmtransfer_doner_patch    = vcat(collect(1:NVEGPOOL_NATVEG), iretransn)
    nf.matrix_ngmtransfer_receiver_patch = fill(ioutn, ngm)
    nf.matrix_ngmtransfer_patch = zeros(np, ngm)
    nf.matrix_ngmturnover_patch = zeros(np, nvegnpool)
    nf.ileaf_to_iout_gm=1;  nf.ileafst_to_iout_gm=2;  nf.ileafxf_to_iout_gm=3
    nf.ifroot_to_iout_gm=4; nf.ifrootst_to_iout_gm=5; nf.ifrootxf_to_iout_gm=6
    nf.ilivestem_to_iout_gm=7;  nf.ilivestemst_to_iout_gm=8;  nf.ilivestemxf_to_iout_gm=9
    nf.ideadstem_to_iout_gm=10; nf.ideadstemst_to_iout_gm=11; nf.ideadstemxf_to_iout_gm=12
    nf.ilivecroot_to_iout_gm=13; nf.ilivecrootst_to_iout_gm=14; nf.ilivecrootxf_to_iout_gm=15
    nf.ideadcroot_to_iout_gm=16; nf.ideadcrootst_to_iout_gm=17; nf.ideadcrootxf_to_iout_gm=18
    nf.iretransn_to_iout_gm=19

    # --- Fire: 2 in-veg live->dead (indices 1,2) + 19 out (pools 1..19 -> iout,
    # indices 3..21). nnfiouttrans=19 trailing → the 2 in-veg transfers are first.
    fi_doner    = vcat([ILIVESTEM, ILIVECROOT], collect(1:NVEGPOOL_NATVEG), iretransn)
    fi_receiver = vcat([IDEADSTEM, IDEADCROOT], fill(ioutn, ngm))
    nf.matrix_nfitransfer_doner_patch    = fi_doner
    nf.matrix_nfitransfer_receiver_patch = fi_receiver
    nf.matrix_nfitransfer_patch = zeros(np, length(fi_doner))
    nf.matrix_nfiturnover_patch = zeros(np, nvegnpool)
    nf.ilivestem_to_ideadstem_fi=1; nf.ilivecroot_to_ideadcroot_fi=2
    nf.ileaf_to_iout_fi=3;  nf.ileafst_to_iout_fi=4;  nf.ileafxf_to_iout_fi=5
    nf.ifroot_to_iout_fi=6; nf.ifrootst_to_iout_fi=7; nf.ifrootxf_to_iout_fi=8
    nf.ilivestem_to_iout_fi=9;  nf.ilivestemst_to_iout_fi=10; nf.ilivestemxf_to_iout_fi=11
    nf.ideadstem_to_iout_fi=12; nf.ideadstemst_to_iout_fi=13; nf.ideadstemxf_to_iout_fi=14
    nf.ilivecroot_to_iout_fi=15; nf.ilivecrootst_to_iout_fi=16; nf.ilivecrootxf_to_iout_fi=17
    nf.ideadcroot_to_iout_fi=18; nf.ideadcrootst_to_iout_fi=19; nf.ideadcrootxf_to_iout_fi=20
    nf.iretransn_to_iout_fi=21
    return nothing
end

# The 12 allocation targets (natveg): (veg pool, npool_to_* flux field).
const _NALLOC_TARGET = Tuple{Int,Symbol}[
    (ILEAF, :npool_to_leafn_patch), (ILEAF_ST, :npool_to_leafn_storage_patch),
    (IFROOT, :npool_to_frootn_patch), (IFROOT_ST, :npool_to_frootn_storage_patch),
    (ILIVESTEM, :npool_to_livestemn_patch), (ILIVESTEM_ST, :npool_to_livestemn_storage_patch),
    (IDEADSTEM, :npool_to_deadstemn_patch), (IDEADSTEM_ST, :npool_to_deadstemn_storage_patch),
    (ILIVECROOT, :npool_to_livecrootn_patch), (ILIVECROOT_ST, :npool_to_livecrootn_storage_patch),
    (IDEADCROOT, :npool_to_deadcrootn_patch), (IDEADCROOT_ST, :npool_to_deadcrootn_storage_patch)]

"""
    cn_veg_matrix_alloc_n!(nf, mask_soilp, bounds_patch)

Assemble the veg-N B-input: matrix_nalloc[p,X] = npool_to_Xn/npool_to_veg (the
allocation fraction, shared by the A-matrix retransn→X transfers), and matrix_Ninput
= npool_to_veg − retransn_to_npool (the sminn-sourced portion; the retransn-sourced
portion enters via the A-matrix). Zero when nothing is allocated.
"""
# The 2 crop grain N allocation targets (reproductive N; 2D fields indexed [p,1]).
const _NALLOC_TARGET_CROP = Tuple{Int,Symbol}[
    (IGRAIN, :npool_to_reproductiven_patch),
    (IGRAIN_ST, :npool_to_reproductiven_storage_patch)]

# N allocation B-input: nalloc = per-pool fraction; Ninput = total − retransn_to_npool.
Base.@kwdef struct _MatNAllocIn{V,M}
    a1::V; a2::V; a3::V; a4::V; a5::V; a6::V; a7::V; a8::V; a9::V; a10::V; a11::V; a12::V
    g1::M; g2::M; r2n::V
end
Adapt.@adapt_structure _MatNAllocIn

@kernel function _mat_alloc_n_kernel!(out::_MatAllocOut, in::_MatNAllocIn, @Const(mask), use_crop::Bool, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask[p]
        alloc = out.alloc; nin = out.cinput
        tot = in.a1[p] + in.a2[p] + in.a3[p] + in.a4[p] + in.a5[p] + in.a6[p] +
              in.a7[p] + in.a8[p] + in.a9[p] + in.a10[p] + in.a11[p] + in.a12[p]
        if use_crop; tot += in.g1[p, 1] + in.g2[p, 1]; end
        if tot > 0
            alloc[p, ILEAF] = in.a1[p] / tot;          alloc[p, ILEAF_ST] = in.a2[p] / tot
            alloc[p, IFROOT] = in.a3[p] / tot;         alloc[p, IFROOT_ST] = in.a4[p] / tot
            alloc[p, ILIVESTEM] = in.a5[p] / tot;      alloc[p, ILIVESTEM_ST] = in.a6[p] / tot
            alloc[p, IDEADSTEM] = in.a7[p] / tot;      alloc[p, IDEADSTEM_ST] = in.a8[p] / tot
            alloc[p, ILIVECROOT] = in.a9[p] / tot;     alloc[p, ILIVECROOT_ST] = in.a10[p] / tot
            alloc[p, IDEADCROOT] = in.a11[p] / tot;    alloc[p, IDEADCROOT_ST] = in.a12[p] / tot
            if use_crop
                alloc[p, IGRAIN] = in.g1[p, 1] / tot;  alloc[p, IGRAIN_ST] = in.g2[p, 1] / tot
            end
        end
        nin[p] = tot - in.r2n[p]
    end
end

function cn_veg_matrix_alloc_n!(nf::CNVegNitrogenFluxData,
        mask_soilp::AbstractVector{Bool}, bounds_patch::UnitRange{Int};
        use_crop::Bool = false)
    fill!(nf.matrix_nalloc_patch, zero(eltype(nf.matrix_nalloc_patch)))
    fill!(nf.matrix_Ninput_patch, zero(eltype(nf.matrix_Ninput_patch)))
    gf(s) = getfield(nf, s)
    out = _MatAllocOut(; alloc=nf.matrix_nalloc_patch, cinput=nf.matrix_Ninput_patch)
    in_ = _MatNAllocIn(;
        a1=gf(_NALLOC_TARGET[1][2]), a2=gf(_NALLOC_TARGET[2][2]), a3=gf(_NALLOC_TARGET[3][2]), a4=gf(_NALLOC_TARGET[4][2]),
        a5=gf(_NALLOC_TARGET[5][2]), a6=gf(_NALLOC_TARGET[6][2]), a7=gf(_NALLOC_TARGET[7][2]), a8=gf(_NALLOC_TARGET[8][2]),
        a9=gf(_NALLOC_TARGET[9][2]), a10=gf(_NALLOC_TARGET[10][2]), a11=gf(_NALLOC_TARGET[11][2]), a12=gf(_NALLOC_TARGET[12][2]),
        g1=gf(_NALLOC_TARGET_CROP[1][2]), g2=gf(_NALLOC_TARGET_CROP[2][2]), r2n=nf.retransn_to_npool_patch)
    backend = _kernel_backend(nf.matrix_nalloc_patch)
    _mat_alloc_n_kernel!(backend)(out, in_, mask_soilp, use_crop, first(bounds_patch), last(bounds_patch); ndrange = last(bounds_patch))
    KA.synchronize(backend)
    return nothing
end

# The 12 retransn→pool supply transfers: (ph index field, receiver pool).
const _RETRANSN_SUPPLY = Tuple{Symbol,Int}[
    (:iretransn_to_ileaf_ph, ILEAF), (:iretransn_to_ileafst_ph, ILEAF_ST),
    (:iretransn_to_ifroot_ph, IFROOT), (:iretransn_to_ifrootst_ph, IFROOT_ST),
    (:iretransn_to_ilivestem_ph, ILIVESTEM), (:iretransn_to_ilivestemst_ph, ILIVESTEM_ST),
    (:iretransn_to_ideadstem_ph, IDEADSTEM), (:iretransn_to_ideadstemst_ph, IDEADSTEM_ST),
    (:iretransn_to_ilivecroot_ph, ILIVECROOT), (:iretransn_to_ilivecrootst_ph, ILIVECROOT_ST),
    (:iretransn_to_ideadcroot_ph, IDEADCROOT), (:iretransn_to_ideadcrootst_ph, IDEADCROOT_ST)]

"""
    cn_veg_matrix_accumulate_phn!(nf, ns, mask_soilp, bounds_patch; dt, matrixcheck=false)

Register every veg-N phenology transfer (storage↔xfer, xfer→pool, livewood turnover,
pool→retransn resorption, litterfall) from the phenology N flux fields, PLUS the 12
retransn→pool supply transfers valued from the allocation nalloc + retransn_to_npool
(requires `cn_veg_matrix_alloc_n!` called first). Requires `cn_veg_matrix_n_topology!`.
"""
Base.@kwdef struct _MatPhNOut{M}
    pht::M; phtn::M
end
Adapt.@adapt_structure _MatPhNOut
Base.@kwdef struct _MatPhNIn{V,M}
    # storage/xfer (12), livewood (2), resorption (4), litterfall (3) N transfer fluxes
    fx1::V; fx2::V; fx3::V; fx4::V; fx5::V; fx6::V; fx7::V; fx8::V; fx9::V; fx10::V; fx11::V; fx12::V
    fw1::V; fw2::V; fr1::V; fr2::V; fr3::V; fr4::V; fll1::V; fll2::V; fll3::V
    na1::V; na2::V; na3::V; na4::V; na5::V; na6::V; na7::V; na8::V; na9::V; na10::V; na11::V; na12::V  # npool_to_* (for npool_to_veg)
    r2n::V                                                                                          # retransn_to_npool
    po1::V; po2::V; po3::V; po4::V; po5::V; po6::V; po7::V; po8::V; po9::V; po10::V
    po11::V; po12::V; po13::V; po14::V; po15::V; po16::V; po17::V; po18::V; po19::V
    nalloc::M; grain_food::M; grain_seed::M; p_repro::M; na_gr1::M; na_gr2::M                        # 2D (crop)
end
Adapt.@adapt_structure _MatPhNIn
Base.@kwdef struct _MatPhNIdx
    x1::Int; x2::Int; x3::Int; x4::Int; x5::Int; x6::Int; x7::Int; x8::Int; x9::Int; x10::Int; x11::Int; x12::Int
    w1::Int; w2::Int; r1::Int; r2::Int; r3::Int; r4::Int; l1::Int; l2::Int; l3::Int; g::Int
    s1::Int; s2::Int; s3::Int; s4::Int; s5::Int; s6::Int; s7::Int; s8::Int; s9::Int; s10::Int; s11::Int; s12::Int
    e::Int; cg1::Int; cg2::Int; iretr::Int; use_crop::Bool; matrixcheck::Bool
end

@kernel function _mat_phn_kernel!(out::_MatPhNOut, in::_MatPhNIn, @Const(mask), ix::_MatPhNIdx, dt, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask[p]
        pht = out.pht; phtn = out.phtn; mc = ix.matrixcheck; iretr = ix.iretr
        # storage/xfer (12)
        _mat_ph_reg!(pht, phtn, p, ix.x1,  in.fx1[p],  ILEAF_ST,      in.po2[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.x2,  in.fx2[p],  ILEAF_XF,      in.po3[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.x3,  in.fx3[p],  IFROOT_ST,     in.po5[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.x4,  in.fx4[p],  IFROOT_XF,     in.po6[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.x5,  in.fx5[p],  ILIVESTEM_ST,  in.po8[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.x6,  in.fx6[p],  ILIVESTEM_XF,  in.po9[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.x7,  in.fx7[p],  IDEADSTEM_ST,  in.po11[p], dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.x8,  in.fx8[p],  IDEADSTEM_XF,  in.po12[p], dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.x9,  in.fx9[p],  ILIVECROOT_ST, in.po14[p], dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.x10, in.fx10[p], ILIVECROOT_XF, in.po15[p], dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.x11, in.fx11[p], IDEADCROOT_ST, in.po17[p], dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.x12, in.fx12[p], IDEADCROOT_XF, in.po18[p], dt, mc)
        # livewood turnover (2)
        _mat_ph_reg!(pht, phtn, p, ix.w1, in.fw1[p], ILIVESTEM,  in.po7[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.w2, in.fw2[p], ILIVECROOT, in.po13[p], dt, mc)
        # pool->retransn resorption (4)
        _mat_ph_reg!(pht, phtn, p, ix.r1, in.fr1[p], ILEAF,      in.po1[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.r2, in.fr2[p], IFROOT,     in.po4[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.r3, in.fr3[p], ILIVESTEM,  in.po7[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.r4, in.fr4[p], ILIVECROOT, in.po13[p], dt, mc)
        # litterfall out (3)
        _mat_ph_reg!(pht, phtn, p, ix.l1, in.fll1[p], ILEAF,     in.po1[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.l2, in.fll2[p], IFROOT,    in.po4[p],  dt, mc)
        _mat_ph_reg!(pht, phtn, p, ix.l3, in.fll3[p], ILIVESTEM, in.po7[p],  dt, mc)
        if ix.use_crop
            _mat_ph_reg!(pht, phtn, p, ix.g, in.grain_food[p, 1] + in.grain_seed[p, 1], IGRAIN, in.p_repro[p, 1], dt, mc)
        end
        # retransn->pool supply, valued from nalloc·r2n (or the iout edge when nothing allocated)
        r2n = in.r2n[p]; retr = in.po19[p]
        n2v = in.na1[p] + in.na2[p] + in.na3[p] + in.na4[p] + in.na5[p] + in.na6[p] +
              in.na7[p] + in.na8[p] + in.na9[p] + in.na10[p] + in.na11[p] + in.na12[p]
        if ix.use_crop
            n2v += in.na_gr1[p, 1] + in.na_gr2[p, 1]
        end
        if n2v > 0
            _mat_ph_reg!(pht, phtn, p, ix.s1,  in.nalloc[p, ILEAF]         * r2n, iretr, retr, dt, mc)
            _mat_ph_reg!(pht, phtn, p, ix.s2,  in.nalloc[p, ILEAF_ST]      * r2n, iretr, retr, dt, mc)
            _mat_ph_reg!(pht, phtn, p, ix.s3,  in.nalloc[p, IFROOT]        * r2n, iretr, retr, dt, mc)
            _mat_ph_reg!(pht, phtn, p, ix.s4,  in.nalloc[p, IFROOT_ST]     * r2n, iretr, retr, dt, mc)
            _mat_ph_reg!(pht, phtn, p, ix.s5,  in.nalloc[p, ILIVESTEM]     * r2n, iretr, retr, dt, mc)
            _mat_ph_reg!(pht, phtn, p, ix.s6,  in.nalloc[p, ILIVESTEM_ST]  * r2n, iretr, retr, dt, mc)
            _mat_ph_reg!(pht, phtn, p, ix.s7,  in.nalloc[p, IDEADSTEM]     * r2n, iretr, retr, dt, mc)
            _mat_ph_reg!(pht, phtn, p, ix.s8,  in.nalloc[p, IDEADSTEM_ST]  * r2n, iretr, retr, dt, mc)
            _mat_ph_reg!(pht, phtn, p, ix.s9,  in.nalloc[p, ILIVECROOT]    * r2n, iretr, retr, dt, mc)
            _mat_ph_reg!(pht, phtn, p, ix.s10, in.nalloc[p, ILIVECROOT_ST] * r2n, iretr, retr, dt, mc)
            _mat_ph_reg!(pht, phtn, p, ix.s11, in.nalloc[p, IDEADCROOT]    * r2n, iretr, retr, dt, mc)
            _mat_ph_reg!(pht, phtn, p, ix.s12, in.nalloc[p, IDEADCROOT_ST] * r2n, iretr, retr, dt, mc)
            if ix.use_crop
                _mat_ph_reg!(pht, phtn, p, ix.cg1, in.nalloc[p, IGRAIN]    * r2n, iretr, retr, dt, mc)
                _mat_ph_reg!(pht, phtn, p, ix.cg2, in.nalloc[p, IGRAIN_ST] * r2n, iretr, retr, dt, mc)
            end
        else
            _mat_ph_reg!(pht, phtn, p, ix.e, r2n, iretr, retr, dt, mc)
        end
    end
end

function cn_veg_matrix_accumulate_phn!(nf::CNVegNitrogenFluxData,
        ns::CNVegNitrogenStateData, mask_soilp::AbstractVector{Bool},
        bounds_patch::UnitRange{Int}; dt::Real, matrixcheck::Bool = false,
        use_crop::Bool = false, nvegnpool::Int = IRETRANSN_NATVEG)
    iretransn = nvegnpool
    fill!(nf.matrix_nphtransfer_patch, zero(eltype(nf.matrix_nphtransfer_patch)))
    fill!(nf.matrix_nphturnover_patch, zero(eltype(nf.matrix_nphturnover_patch)))
    gf(s) = getfield(nf, s)
    out = _MatPhNOut(; pht=nf.matrix_nphtransfer_patch, phtn=nf.matrix_nphturnover_patch)
    in_ = _MatPhNIn(;
        fx1=nf.leafn_storage_to_xfer_patch, fx2=nf.leafn_xfer_to_leafn_patch,
        fx3=nf.frootn_storage_to_xfer_patch, fx4=nf.frootn_xfer_to_frootn_patch,
        fx5=nf.livestemn_storage_to_xfer_patch, fx6=nf.livestemn_xfer_to_livestemn_patch,
        fx7=nf.deadstemn_storage_to_xfer_patch, fx8=nf.deadstemn_xfer_to_deadstemn_patch,
        fx9=nf.livecrootn_storage_to_xfer_patch, fx10=nf.livecrootn_xfer_to_livecrootn_patch,
        fx11=nf.deadcrootn_storage_to_xfer_patch, fx12=nf.deadcrootn_xfer_to_deadcrootn_patch,
        fw1=nf.livestemn_to_deadstemn_patch, fw2=nf.livecrootn_to_deadcrootn_patch,
        fr1=nf.leafn_to_retransn_patch, fr2=nf.frootn_to_retransn_patch,
        fr3=nf.livestemn_to_retransn_patch, fr4=nf.livecrootn_to_retransn_patch,
        fll1=nf.leafn_to_litter_patch, fll2=nf.frootn_to_litter_patch, fll3=nf.livestemn_to_litter_patch,
        na1=gf(_NALLOC_TARGET[1][2]), na2=gf(_NALLOC_TARGET[2][2]), na3=gf(_NALLOC_TARGET[3][2]),
        na4=gf(_NALLOC_TARGET[4][2]), na5=gf(_NALLOC_TARGET[5][2]), na6=gf(_NALLOC_TARGET[6][2]),
        na7=gf(_NALLOC_TARGET[7][2]), na8=gf(_NALLOC_TARGET[8][2]), na9=gf(_NALLOC_TARGET[9][2]),
        na10=gf(_NALLOC_TARGET[10][2]), na11=gf(_NALLOC_TARGET[11][2]), na12=gf(_NALLOC_TARGET[12][2]),
        r2n=nf.retransn_to_npool_patch,
        po1=ns.leafn_patch, po2=ns.leafn_storage_patch, po3=ns.leafn_xfer_patch,
        po4=ns.frootn_patch, po5=ns.frootn_storage_patch, po6=ns.frootn_xfer_patch,
        po7=ns.livestemn_patch, po8=ns.livestemn_storage_patch, po9=ns.livestemn_xfer_patch,
        po10=ns.deadstemn_patch, po11=ns.deadstemn_storage_patch, po12=ns.deadstemn_xfer_patch,
        po13=ns.livecrootn_patch, po14=ns.livecrootn_storage_patch, po15=ns.livecrootn_xfer_patch,
        po16=ns.deadcrootn_patch, po17=ns.deadcrootn_storage_patch, po18=ns.deadcrootn_xfer_patch, po19=ns.retransn_patch,
        nalloc=nf.matrix_nalloc_patch, grain_food=nf.repr_grainn_to_food_patch, grain_seed=nf.repr_grainn_to_seed_patch,
        p_repro=ns.reproductiven_patch, na_gr1=gf(_NALLOC_TARGET_CROP[1][2]), na_gr2=gf(_NALLOC_TARGET_CROP[2][2]))
    ix = _MatPhNIdx(;
        x1=nf.ileafst_to_ileafxf_ph, x2=nf.ileafxf_to_ileaf_ph, x3=nf.ifrootst_to_ifrootxf_ph, x4=nf.ifrootxf_to_ifroot_ph,
        x5=nf.ilivestemst_to_ilivestemxf_ph, x6=nf.ilivestemxf_to_ilivestem_ph, x7=nf.ideadstemst_to_ideadstemxf_ph,
        x8=nf.ideadstemxf_to_ideadstem_ph, x9=nf.ilivecrootst_to_ilivecrootxf_ph, x10=nf.ilivecrootxf_to_ilivecroot_ph,
        x11=nf.ideadcrootst_to_ideadcrootxf_ph, x12=nf.ideadcrootxf_to_ideadcroot_ph,
        w1=nf.ilivestem_to_ideadstem_ph, w2=nf.ilivecroot_to_ideadcroot_ph,
        r1=nf.ileaf_to_iretransn_ph, r2=nf.ifroot_to_iretransn_ph, r3=nf.ilivestem_to_iretransn_ph, r4=nf.ilivecroot_to_iretransn_ph,
        l1=nf.ileaf_to_iout_ph, l2=nf.ifroot_to_iout_ph, l3=nf.ilivestem_to_iout_ph, g=nf.igrain_to_iout_ph,
        s1=nf.iretransn_to_ileaf_ph, s2=nf.iretransn_to_ileafst_ph, s3=nf.iretransn_to_ifroot_ph, s4=nf.iretransn_to_ifrootst_ph,
        s5=nf.iretransn_to_ilivestem_ph, s6=nf.iretransn_to_ilivestemst_ph, s7=nf.iretransn_to_ideadstem_ph, s8=nf.iretransn_to_ideadstemst_ph,
        s9=nf.iretransn_to_ilivecroot_ph, s10=nf.iretransn_to_ilivecrootst_ph, s11=nf.iretransn_to_ideadcroot_ph, s12=nf.iretransn_to_ideadcrootst_ph,
        e=nf.iretransn_to_iout_ph, cg1=nf.iretransn_to_igrain_ph, cg2=nf.iretransn_to_igrainst_ph,
        iretr=iretransn, use_crop=use_crop, matrixcheck=matrixcheck)
    backend = _kernel_backend(nf.matrix_nphtransfer_patch)
    _mat_phn_kernel!(backend)(out, in_, mask_soilp, ix,
        eltype(nf.matrix_nphtransfer_patch)(dt), first(bounds_patch), last(bounds_patch); ndrange = last(bounds_patch))
    KA.synchronize(backend)
    return nothing
end

# The 19 gap-mortality N flux fields (m_*n_to_litter), pools 1..18 + retransn.
const _GMN_FLUX = Symbol[
    :m_leafn_to_litter_patch, :m_leafn_storage_to_litter_patch, :m_leafn_xfer_to_litter_patch,
    :m_frootn_to_litter_patch, :m_frootn_storage_to_litter_patch, :m_frootn_xfer_to_litter_patch,
    :m_livestemn_to_litter_patch, :m_livestemn_storage_to_litter_patch, :m_livestemn_xfer_to_litter_patch,
    :m_deadstemn_to_litter_patch, :m_deadstemn_storage_to_litter_patch, :m_deadstemn_xfer_to_litter_patch,
    :m_livecrootn_to_litter_patch, :m_livecrootn_storage_to_litter_patch, :m_livecrootn_xfer_to_litter_patch,
    :m_deadcrootn_to_litter_patch, :m_deadcrootn_storage_to_litter_patch, :m_deadcrootn_xfer_to_litter_patch,
    :m_retransn_to_litter_patch]

"""
    cn_veg_matrix_accumulate_gmn!(nf, ns, mask_soilp, bounds_patch; dt, matrixcheck=false)

Register the 19 veg-N gap-mortality transfers (pools 1..18 + retransn -> litter).
"""
# N gap-mortality/fire register with separate transfer index idx + turnover column d.
@inline function _mat_gmn_reg!(gmt, gmtn, phtn, p::Int, idx::Int, d::Int, flux, donorval, dt, mc::Bool)
    rate = donorval > 0 ? flux / donorval : zero(flux)
    if mc && (phtn[p, d] + gmtn[p, d] + rate * dt >= one(rate))
        applied = max(zero(rate), (one(rate) - phtn[p, d] - gmtn[p, d]) / dt)
    else
        applied = rate
    end
    @inbounds gmtn[p, d]  = gmtn[p, d]  + applied * dt
    @inbounds gmt[p, idx] = gmt[p, idx] + applied
    return nothing
end

Base.@kwdef struct _MatGmNOut{M}
    gmt::M; gmtn::M; phtn::M
end
Adapt.@adapt_structure _MatGmNOut
Base.@kwdef struct _MatGmNIn{V}
    fl1::V; fl2::V; fl3::V; fl4::V; fl5::V; fl6::V; fl7::V; fl8::V; fl9::V; fl10::V
    fl11::V; fl12::V; fl13::V; fl14::V; fl15::V; fl16::V; fl17::V; fl18::V; fl19::V
    po1::V; po2::V; po3::V; po4::V; po5::V; po6::V; po7::V; po8::V; po9::V; po10::V
    po11::V; po12::V; po13::V; po14::V; po15::V; po16::V; po17::V; po18::V; po19::V
end
Adapt.@adapt_structure _MatGmNIn

@kernel function _mat_gmn_kernel!(out::_MatGmNOut, in::_MatGmNIn, @Const(mask), iretr::Int, mc::Bool, dt, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask[p]
        gmt = out.gmt; gmtn = out.gmtn; phtn = out.phtn
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 1,  1,  in.fl1[p],  in.po1[p],  dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 2,  2,  in.fl2[p],  in.po2[p],  dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 3,  3,  in.fl3[p],  in.po3[p],  dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 4,  4,  in.fl4[p],  in.po4[p],  dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 5,  5,  in.fl5[p],  in.po5[p],  dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 6,  6,  in.fl6[p],  in.po6[p],  dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 7,  7,  in.fl7[p],  in.po7[p],  dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 8,  8,  in.fl8[p],  in.po8[p],  dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 9,  9,  in.fl9[p],  in.po9[p],  dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 10, 10, in.fl10[p], in.po10[p], dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 11, 11, in.fl11[p], in.po11[p], dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 12, 12, in.fl12[p], in.po12[p], dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 13, 13, in.fl13[p], in.po13[p], dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 14, 14, in.fl14[p], in.po14[p], dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 15, 15, in.fl15[p], in.po15[p], dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 16, 16, in.fl16[p], in.po16[p], dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 17, 17, in.fl17[p], in.po17[p], dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 18, 18, in.fl18[p], in.po18[p], dt, mc)
        _mat_gmn_reg!(gmt, gmtn, phtn, p, 19, iretr, in.fl19[p], in.po19[p], dt, mc)   # retransn
    end
end

function cn_veg_matrix_accumulate_gmn!(nf::CNVegNitrogenFluxData,
        ns::CNVegNitrogenStateData, mask_soilp::AbstractVector{Bool},
        bounds_patch::UnitRange{Int}; dt::Real, matrixcheck::Bool = false,
        nvegnpool::Int = IRETRANSN_NATVEG)
    iretransn = nvegnpool
    fill!(nf.matrix_ngmtransfer_patch, zero(eltype(nf.matrix_ngmtransfer_patch)))
    fill!(nf.matrix_ngmturnover_patch, zero(eltype(nf.matrix_ngmturnover_patch)))
    gf(s) = getfield(nf, s)
    out = _MatGmNOut(; gmt=nf.matrix_ngmtransfer_patch, gmtn=nf.matrix_ngmturnover_patch, phtn=nf.matrix_nphturnover_patch)
    in_ = _MatGmNIn(;
        fl1=gf(_GMN_FLUX[1]), fl2=gf(_GMN_FLUX[2]), fl3=gf(_GMN_FLUX[3]), fl4=gf(_GMN_FLUX[4]), fl5=gf(_GMN_FLUX[5]),
        fl6=gf(_GMN_FLUX[6]), fl7=gf(_GMN_FLUX[7]), fl8=gf(_GMN_FLUX[8]), fl9=gf(_GMN_FLUX[9]), fl10=gf(_GMN_FLUX[10]),
        fl11=gf(_GMN_FLUX[11]), fl12=gf(_GMN_FLUX[12]), fl13=gf(_GMN_FLUX[13]), fl14=gf(_GMN_FLUX[14]), fl15=gf(_GMN_FLUX[15]),
        fl16=gf(_GMN_FLUX[16]), fl17=gf(_GMN_FLUX[17]), fl18=gf(_GMN_FLUX[18]), fl19=gf(_GMN_FLUX[19]),
        po1=ns.leafn_patch, po2=ns.leafn_storage_patch, po3=ns.leafn_xfer_patch,
        po4=ns.frootn_patch, po5=ns.frootn_storage_patch, po6=ns.frootn_xfer_patch,
        po7=ns.livestemn_patch, po8=ns.livestemn_storage_patch, po9=ns.livestemn_xfer_patch,
        po10=ns.deadstemn_patch, po11=ns.deadstemn_storage_patch, po12=ns.deadstemn_xfer_patch,
        po13=ns.livecrootn_patch, po14=ns.livecrootn_storage_patch, po15=ns.livecrootn_xfer_patch,
        po16=ns.deadcrootn_patch, po17=ns.deadcrootn_storage_patch, po18=ns.deadcrootn_xfer_patch, po19=ns.retransn_patch)
    backend = _kernel_backend(nf.matrix_ngmtransfer_patch)
    _mat_gmn_kernel!(backend)(out, in_, mask_soilp, iretransn, matrixcheck,
        eltype(nf.matrix_ngmtransfer_patch)(dt), first(bounds_patch), last(bounds_patch); ndrange = last(bounds_patch))
    KA.synchronize(backend)
    return nothing
end

# The 19 fire out-loss N flux fields per pool (to_fire, to_litter_fire summed).
const _FIN_TO_FIRE = Symbol[
    :m_leafn_to_fire_patch, :m_leafn_storage_to_fire_patch, :m_leafn_xfer_to_fire_patch,
    :m_frootn_to_fire_patch, :m_frootn_storage_to_fire_patch, :m_frootn_xfer_to_fire_patch,
    :m_livestemn_to_fire_patch, :m_livestemn_storage_to_fire_patch, :m_livestemn_xfer_to_fire_patch,
    :m_deadstemn_to_fire_patch, :m_deadstemn_storage_to_fire_patch, :m_deadstemn_xfer_to_fire_patch,
    :m_livecrootn_to_fire_patch, :m_livecrootn_storage_to_fire_patch, :m_livecrootn_xfer_to_fire_patch,
    :m_deadcrootn_to_fire_patch, :m_deadcrootn_storage_to_fire_patch, :m_deadcrootn_xfer_to_fire_patch,
    :m_retransn_to_fire_patch]
const _FIN_TO_LITTER = Symbol[
    :m_leafn_to_litter_fire_patch, :m_leafn_storage_to_litter_fire_patch, :m_leafn_xfer_to_litter_fire_patch,
    :m_frootn_to_litter_fire_patch, :m_frootn_storage_to_litter_fire_patch, :m_frootn_xfer_to_litter_fire_patch,
    :m_livestemn_to_litter_fire_patch, :m_livestemn_storage_to_litter_fire_patch, :m_livestemn_xfer_to_litter_fire_patch,
    :m_deadstemn_to_litter_fire_patch, :m_deadstemn_storage_to_litter_fire_patch, :m_deadstemn_xfer_to_litter_fire_patch,
    :m_livecrootn_to_litter_fire_patch, :m_livecrootn_storage_to_litter_fire_patch, :m_livecrootn_xfer_to_litter_fire_patch,
    :m_deadcrootn_to_litter_fire_patch, :m_deadcrootn_storage_to_litter_fire_patch, :m_deadcrootn_xfer_to_litter_fire_patch,
    :m_retransn_to_litter_fire_patch]

"""
    cn_veg_matrix_accumulate_fin!(nf, ns, mask_soilp, bounds_patch; dt, matrixcheck=false)

Register the veg-N fire transfers: 2 in-veg live->dead conversions + 19 out-of-veg
losses (pools 1..18 + retransn, flux = m_*n_to_fire + m_*n_to_litter_fire).
"""
Base.@kwdef struct _MatFiNOut{M}
    fit::M; fitn::M; phtn::M; gmtn::M
end
Adapt.@adapt_structure _MatFiNOut
Base.@kwdef struct _MatFiNIn{V}
    ff1::V; ff2::V; ff3::V; ff4::V; ff5::V; ff6::V; ff7::V; ff8::V; ff9::V; ff10::V
    ff11::V; ff12::V; ff13::V; ff14::V; ff15::V; ff16::V; ff17::V; ff18::V; ff19::V
    fl1::V; fl2::V; fl3::V; fl4::V; fl5::V; fl6::V; fl7::V; fl8::V; fl9::V; fl10::V
    fl11::V; fl12::V; fl13::V; fl14::V; fl15::V; fl16::V; fl17::V; fl18::V; fl19::V
    po1::V; po2::V; po3::V; po4::V; po5::V; po6::V; po7::V; po8::V; po9::V; po10::V
    po11::V; po12::V; po13::V; po14::V; po15::V; po16::V; po17::V; po18::V; po19::V
    f_ls_dead::V; f_lc_dead::V
end
Adapt.@adapt_structure _MatFiNIn
Base.@kwdef struct _MatFiNIdx
    i_ls::Int; i_lc::Int; iretr::Int; matrixcheck::Bool
end

@kernel function _mat_fin_kernel!(out::_MatFiNOut, in::_MatFiNIn, @Const(mask), ix::_MatFiNIdx, dt, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask[p]
        fit = out.fit; fitn = out.fitn; phtn = out.phtn; gmtn = out.gmtn; mc = ix.matrixcheck
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, ix.i_ls, ILIVESTEM,  in.f_ls_dead[p], in.po7[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, ix.i_lc, ILIVECROOT, in.f_lc_dead[p], in.po13[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 3,  1,  in.ff1[p]  + in.fl1[p],  in.po1[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 4,  2,  in.ff2[p]  + in.fl2[p],  in.po2[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 5,  3,  in.ff3[p]  + in.fl3[p],  in.po3[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 6,  4,  in.ff4[p]  + in.fl4[p],  in.po4[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 7,  5,  in.ff5[p]  + in.fl5[p],  in.po5[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 8,  6,  in.ff6[p]  + in.fl6[p],  in.po6[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 9,  7,  in.ff7[p]  + in.fl7[p],  in.po7[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 10, 8,  in.ff8[p]  + in.fl8[p],  in.po8[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 11, 9,  in.ff9[p]  + in.fl9[p],  in.po9[p],  dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 12, 10, in.ff10[p] + in.fl10[p], in.po10[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 13, 11, in.ff11[p] + in.fl11[p], in.po11[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 14, 12, in.ff12[p] + in.fl12[p], in.po12[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 15, 13, in.ff13[p] + in.fl13[p], in.po13[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 16, 14, in.ff14[p] + in.fl14[p], in.po14[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 17, 15, in.ff15[p] + in.fl15[p], in.po15[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 18, 16, in.ff16[p] + in.fl16[p], in.po16[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 19, 17, in.ff17[p] + in.fl17[p], in.po17[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 20, 18, in.ff18[p] + in.fl18[p], in.po18[p], dt, mc)
        _mat_fi_reg!(fit, fitn, phtn, gmtn, p, 21, ix.iretr, in.ff19[p] + in.fl19[p], in.po19[p], dt, mc)  # retransn
    end
end

function cn_veg_matrix_accumulate_fin!(nf::CNVegNitrogenFluxData,
        ns::CNVegNitrogenStateData, mask_soilp::AbstractVector{Bool},
        bounds_patch::UnitRange{Int}; dt::Real, matrixcheck::Bool = false,
        nvegnpool::Int = IRETRANSN_NATVEG)
    iretransn = nvegnpool
    fill!(nf.matrix_nfitransfer_patch, zero(eltype(nf.matrix_nfitransfer_patch)))
    fill!(nf.matrix_nfiturnover_patch, zero(eltype(nf.matrix_nfiturnover_patch)))
    gf(s) = getfield(nf, s)
    out = _MatFiNOut(; fit=nf.matrix_nfitransfer_patch, fitn=nf.matrix_nfiturnover_patch,
        phtn=nf.matrix_nphturnover_patch, gmtn=nf.matrix_ngmturnover_patch)
    in_ = _MatFiNIn(;
        ff1=gf(_FIN_TO_FIRE[1]), ff2=gf(_FIN_TO_FIRE[2]), ff3=gf(_FIN_TO_FIRE[3]), ff4=gf(_FIN_TO_FIRE[4]), ff5=gf(_FIN_TO_FIRE[5]),
        ff6=gf(_FIN_TO_FIRE[6]), ff7=gf(_FIN_TO_FIRE[7]), ff8=gf(_FIN_TO_FIRE[8]), ff9=gf(_FIN_TO_FIRE[9]), ff10=gf(_FIN_TO_FIRE[10]),
        ff11=gf(_FIN_TO_FIRE[11]), ff12=gf(_FIN_TO_FIRE[12]), ff13=gf(_FIN_TO_FIRE[13]), ff14=gf(_FIN_TO_FIRE[14]), ff15=gf(_FIN_TO_FIRE[15]),
        ff16=gf(_FIN_TO_FIRE[16]), ff17=gf(_FIN_TO_FIRE[17]), ff18=gf(_FIN_TO_FIRE[18]), ff19=gf(_FIN_TO_FIRE[19]),
        fl1=gf(_FIN_TO_LITTER[1]), fl2=gf(_FIN_TO_LITTER[2]), fl3=gf(_FIN_TO_LITTER[3]), fl4=gf(_FIN_TO_LITTER[4]), fl5=gf(_FIN_TO_LITTER[5]),
        fl6=gf(_FIN_TO_LITTER[6]), fl7=gf(_FIN_TO_LITTER[7]), fl8=gf(_FIN_TO_LITTER[8]), fl9=gf(_FIN_TO_LITTER[9]), fl10=gf(_FIN_TO_LITTER[10]),
        fl11=gf(_FIN_TO_LITTER[11]), fl12=gf(_FIN_TO_LITTER[12]), fl13=gf(_FIN_TO_LITTER[13]), fl14=gf(_FIN_TO_LITTER[14]), fl15=gf(_FIN_TO_LITTER[15]),
        fl16=gf(_FIN_TO_LITTER[16]), fl17=gf(_FIN_TO_LITTER[17]), fl18=gf(_FIN_TO_LITTER[18]), fl19=gf(_FIN_TO_LITTER[19]),
        po1=ns.leafn_patch, po2=ns.leafn_storage_patch, po3=ns.leafn_xfer_patch,
        po4=ns.frootn_patch, po5=ns.frootn_storage_patch, po6=ns.frootn_xfer_patch,
        po7=ns.livestemn_patch, po8=ns.livestemn_storage_patch, po9=ns.livestemn_xfer_patch,
        po10=ns.deadstemn_patch, po11=ns.deadstemn_storage_patch, po12=ns.deadstemn_xfer_patch,
        po13=ns.livecrootn_patch, po14=ns.livecrootn_storage_patch, po15=ns.livecrootn_xfer_patch,
        po16=ns.deadcrootn_patch, po17=ns.deadcrootn_storage_patch, po18=ns.deadcrootn_xfer_patch, po19=ns.retransn_patch,
        f_ls_dead=nf.m_livestemn_to_deadstemn_fire_patch, f_lc_dead=nf.m_livecrootn_to_deadcrootn_fire_patch)
    ix = _MatFiNIdx(; i_ls=nf.ilivestem_to_ideadstem_fi, i_lc=nf.ilivecroot_to_ideadcroot_fi,
        iretr=iretransn, matrixcheck=matrixcheck)
    backend = _kernel_backend(nf.matrix_nfitransfer_patch)
    _mat_fin_kernel!(backend)(out, in_, mask_soilp, ix,
        eltype(nf.matrix_nfitransfer_patch)(dt), first(bounds_patch), last(bounds_patch); ndrange = last(bounds_patch))
    KA.synchronize(backend)
    return nothing
end
