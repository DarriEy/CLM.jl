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

function cn_veg_matrix_alloc_c!(cf::CNVegCarbonFluxData,
        mask_soilp::AbstractVector{Bool}, bounds_patch::UnitRange{Int};
        use_crop::Bool = false)
    fill!(cf.matrix_alloc_patch, 0.0)
    fill!(cf.matrix_Cinput_patch, 0.0)
    for p in bounds_patch
        mask_soilp[p] || continue
        tot = 0.0
        for (_, fld) in _ALLOC_TARGET; tot += getfield(cf, fld)[p]; end
        if use_crop
            for (_, fld) in _ALLOC_TARGET_CROP; tot += getfield(cf, fld)[p, 1]; end
        end
        cf.matrix_Cinput_patch[p] = tot
        if tot > 0.0
            for (pool, fld) in _ALLOC_TARGET
                cf.matrix_alloc_patch[p, pool] = getfield(cf, fld)[p] / tot
            end
            if use_crop
                for (pool, fld) in _ALLOC_TARGET_CROP
                    cf.matrix_alloc_patch[p, pool] = getfield(cf, fld)[p, 1] / tot
                end
            end
        end
    end
    return nothing
end

"""
    cn_veg_matrix_accumulate_fi_c!(cf, cs, mask_soilp, bounds_patch; dt, matrixcheck=false)

Register the veg-C fire transfers: 2 in-veg live->dead conversions
(m_livestemc_to_deadstemc_fire / m_livecrootc_to_deadcrootc_fire) + 18 out-of-veg
losses (each pool -> iout, flux = m_*_to_fire + m_*_to_litter_fire), via
matrix_update_fic!.
"""
function cn_veg_matrix_accumulate_fi_c!(cf::CNVegCarbonFluxData,
        cs::CNVegCarbonStateData, mask_soilp::AbstractVector{Bool},
        bounds_patch::UnitRange{Int}; dt::Real, matrixcheck::Bool = false)
    fill!(cf.matrix_fitransfer_patch, 0.0)
    fill!(cf.matrix_fiturnover_patch, 0.0)
    for p in bounds_patch
        mask_soilp[p] || continue
        reg(idx, flux, d) = begin
            donor = _vegc_pool_val(cs, d, p)
            rate = donor > 0.0 ? flux / donor : 0.0
            matrix_update_fic!(cf, p, idx, rate, dt; matrixcheck=matrixcheck, acc=true)
        end
        # in-veg live->dead conversions
        reg(cf.ilivestem_to_ideadstem_fi,   cf.m_livestemc_to_deadstemc_fire_patch[p],   ILIVESTEM)
        reg(cf.ilivecroot_to_ideadcroot_fi, cf.m_livecrootc_to_deadcrootc_fire_patch[p], ILIVECROOT)
        # out-of-veg losses (to_fire + to_litter_fire), donor = pool i, index = i+2
        for i in 1:NVEGPOOL_NATVEG
            flux = getfield(cf, _FIC_TO_FIRE[i])[p] + getfield(cf, _FIC_TO_LITTER[i])[p]
            reg(i + 2, flux, i)
        end
    end
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

"""
    cn_veg_matrix_accumulate_gm_c!(cf, cs, mask_soilp, bounds_patch; dt, matrixcheck=false)

Register the 18 veg-C gap-mortality transfers (each veg pool -> litter, a pure
out-of-veg loss) into the matrix from the m_*_to_litter flux fields gap_mortality
already computed. rate = flux/pool; recorded via matrix_update_gmc!.
"""
function cn_veg_matrix_accumulate_gm_c!(cf::CNVegCarbonFluxData,
        cs::CNVegCarbonStateData, mask_soilp::AbstractVector{Bool},
        bounds_patch::UnitRange{Int}; dt::Real, matrixcheck::Bool = false)
    fill!(cf.matrix_gmtransfer_patch, 0.0)
    fill!(cf.matrix_gmturnover_patch, 0.0)
    for p in bounds_patch
        mask_soilp[p] || continue
        for i in 1:NVEGPOOL_NATVEG
            flux = getfield(cf, _GMC_FLUX[i])[p]
            donor = _vegc_pool_val(cs, i, p)
            rate = donor > 0.0 ? flux / donor : 0.0
            matrix_update_gmc!(cf, p, i, rate, dt; matrixcheck=matrixcheck, acc=true)
        end
    end
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

function cn_veg_matrix_alloc_n!(nf::CNVegNitrogenFluxData,
        mask_soilp::AbstractVector{Bool}, bounds_patch::UnitRange{Int};
        use_crop::Bool = false)
    fill!(nf.matrix_nalloc_patch, 0.0)
    fill!(nf.matrix_Ninput_patch, 0.0)
    for p in bounds_patch
        mask_soilp[p] || continue
        tot = 0.0
        for (_, fld) in _NALLOC_TARGET; tot += getfield(nf, fld)[p]; end
        if use_crop
            for (_, fld) in _NALLOC_TARGET_CROP; tot += getfield(nf, fld)[p, 1]; end
        end
        if tot > 0.0
            for (pool, fld) in _NALLOC_TARGET
                nf.matrix_nalloc_patch[p, pool] = getfield(nf, fld)[p] / tot
            end
            if use_crop
                for (pool, fld) in _NALLOC_TARGET_CROP
                    nf.matrix_nalloc_patch[p, pool] = getfield(nf, fld)[p, 1] / tot
                end
            end
        end
        nf.matrix_Ninput_patch[p] = tot - nf.retransn_to_npool_patch[p]
    end
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
function cn_veg_matrix_accumulate_phn!(nf::CNVegNitrogenFluxData,
        ns::CNVegNitrogenStateData, mask_soilp::AbstractVector{Bool},
        bounds_patch::UnitRange{Int}; dt::Real, matrixcheck::Bool = false,
        use_crop::Bool = false, nvegnpool::Int = IRETRANSN_NATVEG)
    iretransn = nvegnpool
    fill!(nf.matrix_nphtransfer_patch, 0.0)
    fill!(nf.matrix_nphturnover_patch, 0.0)
    for p in bounds_patch
        mask_soilp[p] || continue
        reg(idx, flux, d) = begin
            donor = _vegn_pool_val(ns, d, p, iretransn)
            rate = donor > 0.0 ? flux / donor : 0.0
            matrix_update_phn!(nf, p, idx, rate, dt; matrixcheck=matrixcheck, acc=true)
        end
        # storage->xfer, xfer->pool (12)
        reg(nf.ileafst_to_ileafxf_ph,        nf.leafn_storage_to_xfer_patch[p],      ILEAF_ST)
        reg(nf.ileafxf_to_ileaf_ph,          nf.leafn_xfer_to_leafn_patch[p],        ILEAF_XF)
        reg(nf.ifrootst_to_ifrootxf_ph,      nf.frootn_storage_to_xfer_patch[p],     IFROOT_ST)
        reg(nf.ifrootxf_to_ifroot_ph,        nf.frootn_xfer_to_frootn_patch[p],      IFROOT_XF)
        reg(nf.ilivestemst_to_ilivestemxf_ph, nf.livestemn_storage_to_xfer_patch[p], ILIVESTEM_ST)
        reg(nf.ilivestemxf_to_ilivestem_ph,  nf.livestemn_xfer_to_livestemn_patch[p], ILIVESTEM_XF)
        reg(nf.ideadstemst_to_ideadstemxf_ph, nf.deadstemn_storage_to_xfer_patch[p], IDEADSTEM_ST)
        reg(nf.ideadstemxf_to_ideadstem_ph,  nf.deadstemn_xfer_to_deadstemn_patch[p], IDEADSTEM_XF)
        reg(nf.ilivecrootst_to_ilivecrootxf_ph, nf.livecrootn_storage_to_xfer_patch[p], ILIVECROOT_ST)
        reg(nf.ilivecrootxf_to_ilivecroot_ph, nf.livecrootn_xfer_to_livecrootn_patch[p], ILIVECROOT_XF)
        reg(nf.ideadcrootst_to_ideadcrootxf_ph, nf.deadcrootn_storage_to_xfer_patch[p], IDEADCROOT_ST)
        reg(nf.ideadcrootxf_to_ideadcroot_ph, nf.deadcrootn_xfer_to_deadcrootn_patch[p], IDEADCROOT_XF)
        # livewood turnover (2)
        reg(nf.ilivestem_to_ideadstem_ph,    nf.livestemn_to_deadstemn_patch[p],     ILIVESTEM)
        reg(nf.ilivecroot_to_ideadcroot_ph,  nf.livecrootn_to_deadcrootn_patch[p],   ILIVECROOT)
        # pool->retransn resorption (4)
        reg(nf.ileaf_to_iretransn_ph,        nf.leafn_to_retransn_patch[p],          ILEAF)
        reg(nf.ifroot_to_iretransn_ph,       nf.frootn_to_retransn_patch[p],         IFROOT)
        reg(nf.ilivestem_to_iretransn_ph,    nf.livestemn_to_retransn_patch[p],      ILIVESTEM)
        reg(nf.ilivecroot_to_iretransn_ph,   nf.livecrootn_to_retransn_patch[p],     ILIVECROOT)
        # litterfall out (3)
        reg(nf.ileaf_to_iout_ph,             nf.leafn_to_litter_patch[p],            ILEAF)
        reg(nf.ifroot_to_iout_ph,            nf.frootn_to_litter_patch[p],           IFROOT)
        reg(nf.ilivestem_to_iout_ph,         nf.livestemn_to_litter_patch[p],        ILIVESTEM)
        if use_crop
            # grain harvest (grain -> outside veg): food + seed.
            reg(nf.igrain_to_iout_ph,
                nf.repr_grainn_to_food_patch[p, 1] + nf.repr_grainn_to_seed_patch[p, 1], IGRAIN)
        end
        # retransn->pool supply (12, +2 grain for crop) + iretransn->iout edge (alloc-sourced)
        r2n = nf.retransn_to_npool_patch[p]
        npool_to_veg = 0.0
        for (_, fld) in _NALLOC_TARGET; npool_to_veg += getfield(nf, fld)[p]; end
        if use_crop
            for (_, fld) in _NALLOC_TARGET_CROP; npool_to_veg += getfield(nf, fld)[p, 1]; end
        end
        if npool_to_veg > 0.0
            for (idxfld, pool) in _RETRANSN_SUPPLY
                reg(getfield(nf, idxfld), nf.matrix_nalloc_patch[p, pool] * r2n, iretransn)
            end
            if use_crop
                reg(nf.iretransn_to_igrain_ph,   nf.matrix_nalloc_patch[p, IGRAIN]    * r2n, iretransn)
                reg(nf.iretransn_to_igrainst_ph, nf.matrix_nalloc_patch[p, IGRAIN_ST] * r2n, iretransn)
            end
        else
            reg(nf.iretransn_to_iout_ph, r2n, iretransn)
        end
    end
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
function cn_veg_matrix_accumulate_gmn!(nf::CNVegNitrogenFluxData,
        ns::CNVegNitrogenStateData, mask_soilp::AbstractVector{Bool},
        bounds_patch::UnitRange{Int}; dt::Real, matrixcheck::Bool = false,
        nvegnpool::Int = IRETRANSN_NATVEG)
    iretransn = nvegnpool
    fill!(nf.matrix_ngmtransfer_patch, 0.0)
    fill!(nf.matrix_ngmturnover_patch, 0.0)
    n = length(_GMN_FLUX)                       # 19 = 18 natveg + retransn
    for p in bounds_patch
        mask_soilp[p] || continue
        for i in 1:n
            flux = getfield(nf, _GMN_FLUX[i])[p]
            d = i <= NVEGPOOL_NATVEG ? i : iretransn   # transfer 19 is the retransn pool
            donor = _vegn_pool_val(ns, d, p, iretransn)
            rate = donor > 0.0 ? flux / donor : 0.0
            matrix_update_gmn!(nf, p, i, rate, dt; matrixcheck=matrixcheck, acc=true)
        end
    end
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
function cn_veg_matrix_accumulate_fin!(nf::CNVegNitrogenFluxData,
        ns::CNVegNitrogenStateData, mask_soilp::AbstractVector{Bool},
        bounds_patch::UnitRange{Int}; dt::Real, matrixcheck::Bool = false,
        nvegnpool::Int = IRETRANSN_NATVEG)
    iretransn = nvegnpool
    fill!(nf.matrix_nfitransfer_patch, 0.0)
    fill!(nf.matrix_nfiturnover_patch, 0.0)
    n = length(_FIN_TO_FIRE)                    # 19 out-losses = 18 natveg + retransn
    for p in bounds_patch
        mask_soilp[p] || continue
        reg(idx, flux, d) = begin
            donor = _vegn_pool_val(ns, d, p, iretransn)
            rate = donor > 0.0 ? flux / donor : 0.0
            matrix_update_fin!(nf, p, idx, rate, dt; matrixcheck=matrixcheck, acc=true)
        end
        reg(nf.ilivestem_to_ideadstem_fi,   nf.m_livestemn_to_deadstemn_fire_patch[p],   ILIVESTEM)
        reg(nf.ilivecroot_to_ideadcroot_fi, nf.m_livecrootn_to_deadcrootn_fire_patch[p], ILIVECROOT)
        for i in 1:n
            flux = getfield(nf, _FIN_TO_FIRE[i])[p] + getfield(nf, _FIN_TO_LITTER[i])[p]
            d = i <= NVEGPOOL_NATVEG ? i : iretransn   # out-loss 19 is the retransn pool
            reg(i + 2, flux, d)
        end
    end
    return nothing
end
