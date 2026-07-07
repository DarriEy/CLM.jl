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
        # crop adds grain->iout as a 4th out-transfer (index 18).
        push!(doner, IGRAIN); push!(receiver, NVEGPOOL_NATVEG + 1)
        cf.igrain_to_iout_ph = length(doner)
    end
    cf.matrix_phtransfer_doner_patch    = doner
    cf.matrix_phtransfer_receiver_patch = receiver

    # Size the per-step transfer/turnover accumulators (the flux init leaves them
    # empty). np from a phenology flux field the init already sized.
    np = length(cf.leafc_to_litter_patch)
    cf.matrix_phtransfer_patch = zeros(np, length(doner))
    cf.matrix_phturnover_patch = zeros(np, nvegcpool)

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
    iout = NVEGPOOL_NATVEG + 1
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
        bounds_patch::UnitRange{Int}; dt::Real, matrixcheck::Bool = false)

    fill!(cf.matrix_phtransfer_patch, 0.0)
    fill!(cf.matrix_phturnover_patch, 0.0)

    for p in bounds_patch
        mask_soilp[p] || continue
        # register(index, flux, donor_pool_index)
        reg(idx, flux, d) = begin
            donor = _vegc_pool_val(cs, d, p)
            rate = donor > 0.0 ? flux / donor : 0.0
            matrix_update_phc!(cf, p, idx, rate, dt; matrixcheck=matrixcheck, acc=true)
        end
        reg(cf.ileafst_to_ileafxf_ph,        cf.leafc_storage_to_xfer_patch[p],      ILEAF_ST)
        reg(cf.ileafxf_to_ileaf_ph,          cf.leafc_xfer_to_leafc_patch[p],        ILEAF_XF)
        reg(cf.ifrootst_to_ifrootxf_ph,      cf.frootc_storage_to_xfer_patch[p],     IFROOT_ST)
        reg(cf.ifrootxf_to_ifroot_ph,        cf.frootc_xfer_to_frootc_patch[p],      IFROOT_XF)
        reg(cf.ilivestemst_to_ilivestemxf_ph, cf.livestemc_storage_to_xfer_patch[p], ILIVESTEM_ST)
        reg(cf.ilivestemxf_to_ilivestem_ph,  cf.livestemc_xfer_to_livestemc_patch[p], ILIVESTEM_XF)
        reg(cf.ideadstemst_to_ideadstemxf_ph, cf.deadstemc_storage_to_xfer_patch[p], IDEADSTEM_ST)
        reg(cf.ideadstemxf_to_ideadstem_ph,  cf.deadstemc_xfer_to_deadstemc_patch[p], IDEADSTEM_XF)
        reg(cf.ilivecrootst_to_ilivecrootxf_ph, cf.livecrootc_storage_to_xfer_patch[p], ILIVECROOT_ST)
        reg(cf.ilivecrootxf_to_ilivecroot_ph, cf.livecrootc_xfer_to_livecrootc_patch[p], ILIVECROOT_XF)
        reg(cf.ideadcrootst_to_ideadcrootxf_ph, cf.deadcrootc_storage_to_xfer_patch[p], IDEADCROOT_ST)
        reg(cf.ideadcrootxf_to_ideadcroot_ph, cf.deadcrootc_xfer_to_deadcrootc_patch[p], IDEADCROOT_XF)
        reg(cf.ilivestem_to_ideadstem_ph,    cf.livestemc_to_deadstemc_patch[p],     ILIVESTEM)
        reg(cf.ilivecroot_to_ideadcroot_ph,  cf.livecrootc_to_deadcrootc_patch[p],   ILIVECROOT)
        reg(cf.ileaf_to_iout_ph,             cf.leafc_to_litter_patch[p],            ILEAF)
        reg(cf.ifroot_to_iout_ph,            cf.frootc_to_litter_patch[p],           IFROOT)
        reg(cf.ilivestem_to_iout_ph,         cf.livestemc_to_litter_patch[p],        ILIVESTEM)
    end
    return nothing
end
