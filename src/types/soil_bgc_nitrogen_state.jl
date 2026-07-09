# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemNitrogenStateType.F90
# Soil biogeochemistry nitrogen state data type allocation and initialization
# ==========================================================================

"""
    SoilBiogeochemNitrogenStateData

Soil biogeochemistry nitrogen state data structure. Holds decomposition N pools
at column and gridcell levels.

Ported from `soilbiogeochem_nitrogenstate_type` in `SoilBiogeochemNitrogenStateType.F90`.
"""
Base.@kwdef mutable struct SoilBiogeochemNitrogenStateData{FT<:Real,
                                               V<:AbstractVector{FT},
                                               M<:AbstractMatrix{FT},
                                               A3<:AbstractArray{FT,3}}
    # --- Vertically-resolved decomposition pools (col × nlev × npools) ---
    decomp_npools_vr_col              ::A3 = Array{Float64}(undef, 0, 0, 0)  # (gN/m3)
    decomp0_npools_vr_col             ::A3 = Array{Float64}(undef, 0, 0, 0)  # (gN/m3) baseline (initial value of this year)
    decomp_npools_vr_SASUsave_col     ::A3 = Array{Float64}(undef, 0, 0, 0)  # (gN/m3) SASU save
    decomp_soiln_vr_col               ::M = Matrix{Float64}(undef, 0, 0)    # (gN/m3) total soil N vr

    # --- Vertically-resolved mineral N and truncation ---
    sminn_vr_col                      ::M = Matrix{Float64}(undef, 0, 0)    # (gN/m3) soil mineral N vr
    ntrunc_vr_col                     ::M = Matrix{Float64}(undef, 0, 0)    # (gN/m3) N truncation vr

    # --- Nitrification/denitrification pools ---
    smin_no3_vr_col                   ::M = Matrix{Float64}(undef, 0, 0)    # (gN/m3) soil mineral NO3 vr
    smin_no3_col                      ::V = Float64[]                        # (gN/m2) soil mineral NO3
    smin_nh4_vr_col                   ::M = Matrix{Float64}(undef, 0, 0)    # (gN/m3) soil mineral NH4 vr
    smin_nh4_col                      ::V = Float64[]                        # (gN/m2) soil mineral NH4

    # --- Summary (diagnostic) state variables ---
    decomp_npools_col                 ::M = Matrix{Float64}(undef, 0, 0)    # (gN/m2) decomp N pools
    decomp_npools_1m_col              ::M = Matrix{Float64}(undef, 0, 0)    # (gN/m2) decomp N pools to 1m
    sminn_col                         ::V = Float64[]                        # (gN/m2) soil mineral N
    ntrunc_col                        ::V = Float64[]                        # (gN/m2) N truncation
    cwdn_col                          ::V = Float64[]                        # (gN/m2) coarse woody debris N
    totlitn_col                       ::V = Float64[]                        # (gN/m2) total litter N
    totmicn_col                       ::V = Float64[]                        # (gN/m2) total microbial N
    totsomn_col                       ::V = Float64[]                        # (gN/m2) total SOM N
    totlitn_1m_col                    ::V = Float64[]                        # (gN/m2) total litter N to 1m
    totsomn_1m_col                    ::V = Float64[]                        # (gN/m2) total SOM N to 1m
    dyn_nbal_adjustments_col          ::V = Float64[]                        # (gN/m2) dynamic column N adjustments
    dyn_no3bal_adjustments_col        ::V = Float64[]                        # (gN/m2) dynamic column NO3 adjustments
    dyn_nh4bal_adjustments_col        ::V = Float64[]                        # (gN/m2) dynamic column NH4 adjustments

    # --- Scalars ---
    totvegcthresh                     ::Float64 = NaN       # threshold for zeroing decomp pools

    # --- Nitrogen totals (includes soil, veg) ---
    totn_col                          ::V = Float64[]  # (gN/m2) total column N
    totecosysn_col                    ::V = Float64[]  # (gN/m2) total ecosystem N
    totn_grc                          ::V = Float64[]  # (gN/m2) total gridcell N

    # --- Matrix-CN fields ---
    matrix_cap_decomp_npools_col      ::M = Matrix{Float64}(undef, 0, 0)    # (gN/m2)
    matrix_cap_decomp_npools_vr_col   ::A3 = Array{Float64}(undef, 0, 0, 0)  # (gN/m3)
    in_nacc                           ::M = Matrix{Float64}(undef, 0, 0)    # (gN/m3/yr)
    in_nacc_2d                        ::A3 = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/yr)
    tran_nacc                         ::A3 = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/yr)
    vert_up_tran_nacc                 ::A3 = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/yr)
    vert_down_tran_nacc               ::A3 = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/yr)
    exit_nacc                         ::A3 = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/yr)
    hori_tran_nacc                    ::A3 = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/yr)
end

SoilBiogeochemNitrogenStateData{FT}(; kwargs...) where {FT<:Real} =
    SoilBiogeochemNitrogenStateData{FT, Vector{FT}, Matrix{FT}, Array{FT,3}}(; kwargs...)
Adapt.@adapt_structure SoilBiogeochemNitrogenStateData


# ---------------------------------------------------------------------------
# Helper constructors (reuse from carbon state if already defined)
# ---------------------------------------------------------------------------
if !@isdefined(nanvec)
    nanvec(n) = fill(NaN, n)
end
if !@isdefined(nanmat)
    nanmat(r, c) = fill(NaN, r, c)
end
if !@isdefined(nan3d)
    nan3d(a, b, c) = fill(NaN, a, b, c)
end

"""
    soil_bgc_nitrogen_state_init!(ns, nc, ng, nlevdecomp_full, ndecomp_pools;
                                   nlevdecomp=nlevdecomp_full,
                                   ndecomp_cascade_transitions=0,
                                   use_soil_matrixcn=false)

Allocate all fields of `SoilBiogeochemNitrogenStateData`.
Corresponds to `InitAllocate` in the Fortran source.
"""
function soil_bgc_nitrogen_state_init!(ns::SoilBiogeochemNitrogenStateData,
                                        nc::Int, ng::Int,
                                        nlevdecomp_full::Int,
                                        ndecomp_pools::Int;
                                        nlevdecomp::Int=nlevdecomp_full,
                                        ndecomp_cascade_transitions::Int=0,
                                        use_soil_matrixcn::Bool=false)

    ns.totvegcthresh = NaN

    # --- Vertically-resolved mineral N and truncation (col × nlev) ---
    ns.sminn_vr_col     = nanmat(nc, nlevdecomp_full)
    ns.ntrunc_vr_col    = nanmat(nc, nlevdecomp_full)
    ns.smin_no3_vr_col  = nanmat(nc, nlevdecomp_full)
    ns.smin_nh4_vr_col  = nanmat(nc, nlevdecomp_full)

    # --- Column-level 1D (nitrif/denitrif) ---
    ns.smin_no3_col     = nanvec(nc)
    ns.smin_nh4_col     = nanvec(nc)

    # --- Column-level 1D ---
    ns.cwdn_col         = nanvec(nc)
    ns.sminn_col        = nanvec(nc)
    ns.ntrunc_col       = nanvec(nc)
    ns.totlitn_col      = nanvec(nc)
    ns.totmicn_col      = nanvec(nc)
    ns.totsomn_col      = nanvec(nc)
    ns.totlitn_1m_col   = nanvec(nc)
    ns.totsomn_1m_col   = nanvec(nc)
    ns.dyn_nbal_adjustments_col   = nanvec(nc)
    ns.dyn_no3bal_adjustments_col = nanvec(nc)
    ns.dyn_nh4bal_adjustments_col = nanvec(nc)

    # --- Column × npools 2D ---
    ns.decomp_npools_col    = nanmat(nc, ndecomp_pools)
    ns.decomp_npools_1m_col = nanmat(nc, ndecomp_pools)

    if use_soil_matrixcn
        ns.matrix_cap_decomp_npools_col = nanmat(nc, ndecomp_pools)
    end

    # --- Vertically-resolved 3D (col × nlev × npools) ---
    ns.decomp_npools_vr_col = nan3d(nc, nlevdecomp_full, ndecomp_pools)

    if use_soil_matrixcn
        ns.matrix_cap_decomp_npools_vr_col = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        ns.decomp0_npools_vr_col           = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        ns.decomp_npools_vr_SASUsave_col   = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        ns.in_nacc      = nanmat(nc, nlevdecomp * ndecomp_pools)
        ns.tran_nacc    = nan3d(nc, nlevdecomp * ndecomp_pools, nlevdecomp * ndecomp_pools)
        ns.in_nacc_2d   = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        ns.vert_up_tran_nacc   = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        ns.vert_down_tran_nacc = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        ns.exit_nacc           = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        ns.hori_tran_nacc      = nan3d(nc, nlevdecomp_full, ndecomp_cascade_transitions)
    end

    # --- Vertically-resolved decomposing total soil N (col × nlev) ---
    ns.decomp_soiln_vr_col = nanmat(nc, nlevdecomp_full)

    # --- Totals ---
    ns.totn_col       = nanvec(nc)
    ns.totecosysn_col = nanvec(nc)
    ns.totn_grc       = nanvec(ng)

    return nothing
end

"""
    soil_bgc_nitrogen_state_set_values!(ns, mask_col, value_column;
                                         nlevdecomp_full=..., nlevdecomp=...,
                                         ndecomp_pools=..., ndecomp_cascade_transitions=0,
                                         use_soil_matrixcn=false,
                                         use_nitrif_denitrif=false)

Set nitrogen state variables for columns matching `mask_col`.
Corresponds to `SetValues` in the Fortran source.
"""
function soil_bgc_nitrogen_state_set_values!(ns::SoilBiogeochemNitrogenStateData,
                                              mask_col::BitVector,
                                              value_column::Real;
                                              nlevdecomp_full::Int=size(ns.ntrunc_vr_col, 2),
                                              nlevdecomp::Int=nlevdecomp_full,
                                              ndecomp_pools::Int=size(ns.decomp_npools_col, 2),
                                              ndecomp_cascade_transitions::Int=0,
                                              use_soil_matrixcn::Bool=false,
                                              use_nitrif_denitrif::Bool=false)
    for i in eachindex(mask_col)
        mask_col[i] || continue
        ns.sminn_col[i]       = value_column
        ns.ntrunc_col[i]      = value_column
        ns.cwdn_col[i]        = value_column
        if use_nitrif_denitrif
            ns.smin_no3_col[i] = value_column
            ns.smin_nh4_col[i] = value_column
        end
        ns.totlitn_col[i]     = value_column
        ns.totmicn_col[i]     = value_column
        ns.totsomn_col[i]     = value_column
        ns.totsomn_1m_col[i]  = value_column
        ns.totlitn_1m_col[i]  = value_column
    end

    for j in 1:nlevdecomp_full
        for i in eachindex(mask_col)
            mask_col[i] || continue
            ns.sminn_vr_col[i, j]  = value_column
            ns.ntrunc_vr_col[i, j] = value_column
            if use_nitrif_denitrif
                ns.smin_no3_vr_col[i, j] = value_column
                ns.smin_nh4_vr_col[i, j] = value_column
            end
        end
    end

    for k in 1:ndecomp_pools
        for i in eachindex(mask_col)
            mask_col[i] || continue
            ns.decomp_npools_col[i, k]    = value_column
            ns.decomp_npools_1m_col[i, k] = value_column
            if use_soil_matrixcn
                ns.matrix_cap_decomp_npools_col[i, k] = value_column
            end
        end
    end

    for j in 1:nlevdecomp_full
        for k in 1:ndecomp_pools
            for i in eachindex(mask_col)
                mask_col[i] || continue
                ns.decomp_npools_vr_col[i, j, k] = value_column
                if use_soil_matrixcn
                    ns.matrix_cap_decomp_npools_vr_col[i, j, k] = value_column
                    ns.decomp0_npools_vr_col[i, j, k] = value_column
                end
            end
        end
    end

    for i in eachindex(mask_col)
        mask_col[i] || continue
        ns.totecosysn_col[i] = value_column
        ns.totn_col[i]       = value_column
    end

    if use_soil_matrixcn
        for j in 1:nlevdecomp
            for k in 1:ndecomp_pools
                for i in eachindex(mask_col)
                    mask_col[i] || continue
                    ns.in_nacc_2d[i, j, k]          = value_column
                    ns.vert_up_tran_nacc[i, j, k]   = value_column
                    ns.vert_down_tran_nacc[i, j, k] = value_column
                    ns.exit_nacc[i, j, k]            = value_column
                end
            end
            for k in 1:ndecomp_cascade_transitions
                for i in eachindex(mask_col)
                    mask_col[i] || continue
                    ns.hori_tran_nacc[i, j, k] = value_column
                end
            end
        end
    end

    return nothing
end

"""
    soil_bgc_nitrogen_state_init_cold!(ns, bounds_col;
                                        nlevdecomp=..., nlevdecomp_full=...,
                                        ndecomp_pools=..., ndecomp_cascade_transitions=0,
                                        decomp_cpools_vr_col=..., decomp_cpools_col=...,
                                        decomp_cpools_1m_col=..., initial_cn_ratio=...,
                                        use_soil_matrixcn=false, use_nitrif_denitrif=false,
                                        mask_soil_crop=nothing)

Cold-start initialization of soil N pools from C pools via C:N ratios.
Corresponds to `InitCold` in the Fortran source.
"""
function soil_bgc_nitrogen_state_init_cold!(ns::SoilBiogeochemNitrogenStateData,
                                             bounds_col::UnitRange{Int};
                                             nlevdecomp::Int=size(ns.decomp_npools_vr_col, 2),
                                             nlevdecomp_full::Int=nlevdecomp,
                                             ndecomp_pools::Int=size(ns.decomp_npools_vr_col, 3),
                                             ndecomp_cascade_transitions::Int=0,
                                             decomp_cpools_vr_col::Array{<:Real,3}=zeros(0,0,0),
                                             decomp_cpools_col::Matrix{<:Real}=zeros(0,0),
                                             decomp_cpools_1m_col::Matrix{<:Real}=zeros(0,0),
                                             initial_cn_ratio::Vector{<:Real}=ones(ndecomp_pools),
                                             use_soil_matrixcn::Bool=false,
                                             use_nitrif_denitrif::Bool=false,
                                             mask_soil_crop::Union{BitVector,Nothing}=nothing)
    for c in bounds_col
        # matrix-spinup
        if use_soil_matrixcn
            ns.in_nacc[c, :] .= 0.0
        end

        is_soil_crop = mask_soil_crop === nothing || mask_soil_crop[c]

        if is_soil_crop
            # column nitrogen state variables
            ns.ntrunc_col[c] = 0.0
            ns.sminn_col[c]  = 0.0

            for j in 1:nlevdecomp
                for k in 1:ndecomp_pools
                    ns.decomp_npools_vr_col[c, j, k] = decomp_cpools_vr_col[c, j, k] / initial_cn_ratio[k]
                    if use_soil_matrixcn
                        ns.matrix_cap_decomp_npools_vr_col[c, j, k] = decomp_cpools_vr_col[c, j, k] / initial_cn_ratio[k]
                        ns.in_nacc_2d[c, j, k]          = 0.0
                        ns.vert_up_tran_nacc[c, j, k]   = 0.0
                        ns.vert_down_tran_nacc[c, j, k] = 0.0
                        ns.exit_nacc[c, j, k]           = 0.0
                        ns.decomp0_npools_vr_col[c, j, k]       = max(ns.decomp_npools_vr_col[c, j, k], 1.0e-30)
                        ns.decomp_npools_vr_SASUsave_col[c, j, k] = 0.0
                    end
                end
                if use_soil_matrixcn
                    for k in 1:ndecomp_cascade_transitions
                        ns.hori_tran_nacc[c, j, k] = 0.0
                    end
                end

                ns.sminn_vr_col[c, j]  = 0.0
                ns.ntrunc_vr_col[c, j] = 0.0
            end

            if nlevdecomp > 1
                for j in (nlevdecomp+1):nlevdecomp_full
                    for k in 1:ndecomp_pools
                        ns.decomp_npools_vr_col[c, j, k] = 0.0
                        if use_soil_matrixcn
                            ns.matrix_cap_decomp_npools_vr_col[c, j, k] = 0.0
                            ns.in_nacc_2d[c, j, k]          = 0.0
                            ns.vert_up_tran_nacc[c, j, k]   = 0.0
                            ns.vert_down_tran_nacc[c, j, k] = 0.0
                            ns.exit_nacc[c, j, k]           = 0.0
                            ns.decomp0_npools_vr_col[c, j, k] = ns.decomp_npools_vr_col[c, j, k]
                        end
                    end
                    if use_soil_matrixcn
                        for k in 1:ndecomp_cascade_transitions
                            ns.hori_tran_nacc[c, j, k] = 0.0
                        end
                    end
                    ns.sminn_vr_col[c, j]  = 0.0
                    ns.ntrunc_vr_col[c, j] = 0.0
                end
            end

            for k in 1:ndecomp_pools
                ns.decomp_npools_col[c, k]    = decomp_cpools_col[c, k] / initial_cn_ratio[k]
                ns.decomp_npools_1m_col[c, k] = decomp_cpools_1m_col[c, k] / initial_cn_ratio[k]
                if use_soil_matrixcn
                    ns.matrix_cap_decomp_npools_col[c, k] = decomp_cpools_col[c, k] / initial_cn_ratio[k]
                end
            end

            if use_nitrif_denitrif
                for j in 1:nlevdecomp_full
                    ns.smin_nh4_vr_col[c, j] = 0.0
                    ns.smin_no3_vr_col[c, j] = 0.0
                end
                ns.smin_nh4_col[c] = 0.0
                ns.smin_no3_col[c] = 0.0
            end

            ns.totlitn_col[c]    = 0.0
            ns.totmicn_col[c]    = 0.0
            ns.totsomn_col[c]    = 0.0
            ns.totlitn_1m_col[c] = 0.0
            ns.totsomn_1m_col[c] = 0.0
            ns.cwdn_col[c]       = 0.0

            # total nitrogen pools
            ns.totecosysn_col[c] = 0.0
            ns.totn_col[c]       = 0.0
        end
    end

    return nothing
end

# ==========================================================================
# GPU kernelization of the nitrogen summary reduction — see the carbon summary
# for the design. One thread per column, all vertical/pool integrals via LOCAL
# scalar accumulators (an in-place `arr[c,l] += …` in a nested loop miscompiles on
# the KA CPU backend under --check-bounds=yes). Byte-identical per-column sequence.
# ==========================================================================
struct _SBNSumView{V,M,A3}
    decomp_npools_vr::A3
    matrix_cap_decomp_npools_vr::A3
    sminn_vr::M
    ntrunc_vr::M
    smin_no3_vr::M
    smin_nh4_vr::M
    decomp_soiln_vr::M
    decomp_npools::M
    decomp_npools_1m::M
    matrix_cap_decomp_npools::M
    smin_no3::V
    smin_nh4::V
    sminn::V
    ntrunc::V
    cwdn::V
    totlitn::V
    totmicn::V
    totsomn::V
    totlitn_1m::V
    totsomn_1m::V
    totn::V
    totecosysn::V
end
Adapt.@adapt_structure _SBNSumView

@kernel function _sbn_summary_kernel!(@Const(mask), b,
        @Const(dzsoi), @Const(zisoi),
        @Const(is_litter), @Const(is_soil), @Const(is_microbe), @Const(is_cwd),
        @Const(totn_p2c), @Const(totvegn), @Const(is_fates),
        lo::Int, hi::Int, nlevdecomp::Int, nlevdecomp_full::Int, ndecomp_pools::Int,
        use_soil_matrixcn::Bool, use_nitrif_denitrif::Bool)
    c = @index(Global)
    @inbounds if lo <= c <= hi && mask[c]
        T = eltype(b.decomp_npools)

        # Vertically integrate NO3/NH4 pools
        if use_nitrif_denitrif
            s3 = zero(T); s4 = zero(T)
            for j in 1:nlevdecomp
                s3 += b.smin_no3_vr[c, j] * dzsoi[j]
                s4 += b.smin_nh4_vr[c, j] * dzsoi[j]
            end
            b.smin_no3[c] = s3
            b.smin_nh4[c] = s4
        end

        # Vertically integrate each decomposing N pool
        for l in 1:ndecomp_pools
            s = zero(T); sm = zero(T)
            for j in 1:nlevdecomp
                s += b.decomp_npools_vr[c, j, l] * dzsoi[j]
                if use_soil_matrixcn
                    sm += b.matrix_cap_decomp_npools_vr[c, j, l] * dzsoi[j]
                end
            end
            b.decomp_npools[c, l] = s
            if use_soil_matrixcn
                b.matrix_cap_decomp_npools[c, l] = sm
            end
        end

        # Vertically integrate to 1 meter (+ soil-N vr + 1m litter/SOM totals)
        if nlevdecomp > 1
            maxdepth = one(T)
            for l in 1:ndecomp_pools
                s = zero(T)
                for j in 1:nlevdecomp
                    if zisoi[j+1] <= maxdepth
                        s += b.decomp_npools_vr[c, j, l] * dzsoi[j]
                    elseif zisoi[j] < maxdepth
                        s += b.decomp_npools_vr[c, j, l] * (maxdepth - zisoi[j])
                    end
                end
                b.decomp_npools_1m[c, l] = s
            end

            if nlevdecomp_full > 1
                for j in 1:nlevdecomp
                    s = zero(T)
                    for l in 1:ndecomp_pools
                        if is_soil[l]
                            s += b.decomp_npools_vr[c, j, l]
                        end
                    end
                    b.decomp_soiln_vr[c, j] = s
                end
            end

            let s = zero(T)
                for l in 1:ndecomp_pools
                    if is_litter[l]
                        s += b.decomp_npools_1m[c, l]
                    end
                end
                b.totlitn_1m[c] = s
            end
            let s = zero(T)
                for l in 1:ndecomp_pools
                    if is_soil[l]
                        s += b.decomp_npools_1m[c, l]
                    end
                end
                b.totsomn_1m[c] = s
            end
        end

        # Total litter nitrogen
        let s = zero(T)
            for l in 1:ndecomp_pools
                if is_litter[l]
                    s += b.decomp_npools[c, l]
                end
            end
            b.totlitn[c] = s
        end

        # Total microbial nitrogen
        let s = zero(T)
            for l in 1:ndecomp_pools
                if is_microbe[l]
                    s += b.decomp_npools[c, l]
                end
            end
            b.totmicn[c] = s
        end

        # Total SOM nitrogen
        let s = zero(T)
            for l in 1:ndecomp_pools
                if is_soil[l]
                    s += b.decomp_npools[c, l]
                end
            end
            b.totsomn[c] = s
        end

        # Total sminn
        let s = zero(T)
            for j in 1:nlevdecomp
                s += b.sminn_vr[c, j] * dzsoi[j]
            end
            b.sminn[c] = s
        end

        # Total ntrunc
        let s = zero(T)
            for j in 1:nlevdecomp
                s += b.ntrunc_vr[c, j] * dzsoi[j]
            end
            b.ntrunc[c] = s
        end

        # CWD nitrogen, ecosystem N, total column N
        if is_fates[c]
            b.cwdn[c] = zero(T)
            tvegn   = zero(T)
            ecovegn = zero(T)
        else
            s = zero(T)
            for l in 1:ndecomp_pools
                if is_cwd[l]
                    s += b.decomp_npools[c, l]
                end
            end
            b.cwdn[c] = s
            tvegn   = totn_p2c[c]
            ecovegn = totvegn[c]
        end

        b.totecosysn[c] = b.cwdn[c] + b.totlitn[c] + b.totmicn[c] + b.totsomn[c] +
                          b.sminn[c] + ecovegn
        b.totn[c] = b.cwdn[c] + b.totlitn[c] + b.totmicn[c] + b.totsomn[c] +
                    b.sminn[c] + b.ntrunc[c] + tvegn
    end
end

"""
    soil_bgc_nitrogen_state_summary!(ns, mask_allc, bounds_col;
                                      nlevdecomp, nlevdecomp_full, ndecomp_pools,
                                      dzsoi_decomp_vals, zisoi_vals,
                                      is_litter, is_soil, is_microbe, is_cwd,
                                      totn_p2c_col, totvegn_col, is_fates_col,
                                      use_soil_matrixcn, use_fates_bgc,
                                      use_nitrif_denitrif, mask_bgc_soilc)

Perform column-level nitrogen summary calculations.
Corresponds to `Summary` in the Fortran source.
"""
function soil_bgc_nitrogen_state_summary!(ns::SoilBiogeochemNitrogenStateData,
                                           mask_allc::AbstractVector{Bool},
                                           bounds_col::UnitRange{Int};
                                           nlevdecomp::Int,
                                           nlevdecomp_full::Int=nlevdecomp,
                                           ndecomp_pools::Int,
                                           dzsoi_decomp_vals::AbstractVector{<:Real},
                                           zisoi_vals::AbstractVector{<:Real}=Float64[],
                                           is_litter::AbstractVector{Bool}=falses(ndecomp_pools),
                                           is_soil::AbstractVector{Bool}=falses(ndecomp_pools),
                                           is_microbe::AbstractVector{Bool}=falses(ndecomp_pools),
                                           is_cwd::AbstractVector{Bool}=falses(ndecomp_pools),
                                           totn_p2c_col::AbstractVector{<:Real}=zeros(length(mask_allc)),
                                           totvegn_col::AbstractVector{<:Real}=zeros(length(mask_allc)),
                                           is_fates_col::Union{AbstractVector{Bool},Nothing}=nothing,
                                           use_soil_matrixcn::Bool=false,
                                           use_fates_bgc::Bool=false,
                                           use_nitrif_denitrif::Bool=false,
                                           mask_bgc_soilc::Union{AbstractVector{Bool},Nothing}=nothing)

    isempty(bounds_col) && return nothing
    FT = eltype(ns.decomp_npools_vr_col)
    proto = ns.decomp_npools_vr_col
    bv(x) = x isa BitVector ? collect(Bool, x) : x
    isf = is_fates_col === nothing ? falses(length(mask_allc)) : is_fates_col
    zis = isempty(zisoi_vals) ? zeros(nlevdecomp + 1) : zisoi_vals

    mask_k  = _to_backend_like(proto, FT, bv(mask_allc))
    dz_k    = _to_backend_like(proto, FT, dzsoi_decomp_vals)
    zi_k    = _to_backend_like(proto, FT, zis)
    isl_k   = _to_backend_like(proto, FT, bv(is_litter))
    iss_k   = _to_backend_like(proto, FT, bv(is_soil))
    ism_k   = _to_backend_like(proto, FT, bv(is_microbe))
    isc_k   = _to_backend_like(proto, FT, bv(is_cwd))
    tp2n_k  = _to_backend_like(proto, FT, totn_p2c_col)
    tvegn_k = _to_backend_like(proto, FT, totvegn_col)
    isf_k   = _to_backend_like(proto, FT, bv(isf))

    b = _SBNSumView(ns.decomp_npools_vr_col, ns.matrix_cap_decomp_npools_vr_col,
        ns.sminn_vr_col, ns.ntrunc_vr_col, ns.smin_no3_vr_col, ns.smin_nh4_vr_col,
        ns.decomp_soiln_vr_col, ns.decomp_npools_col, ns.decomp_npools_1m_col,
        ns.matrix_cap_decomp_npools_col, ns.smin_no3_col, ns.smin_nh4_col,
        ns.sminn_col, ns.ntrunc_col, ns.cwdn_col, ns.totlitn_col, ns.totmicn_col,
        ns.totsomn_col, ns.totlitn_1m_col, ns.totsomn_1m_col, ns.totn_col, ns.totecosysn_col)

    _launch!(_sbn_summary_kernel!, mask_k, b, dz_k, zi_k, isl_k, iss_k, ism_k, isc_k,
        tp2n_k, tvegn_k, isf_k, first(bounds_col), last(bounds_col),
        nlevdecomp, nlevdecomp_full, ndecomp_pools, use_soil_matrixcn, use_nitrif_denitrif;
        ndrange=length(mask_k))
    return nothing
end

"""
    soil_bgc_nitrogen_state_set_totvegcthresh!(ns, totvegcthresh)

Set the total vegetation carbon threshold for spinup.
Corresponds to `SetTotVgCThresh` in the Fortran source.
"""
function soil_bgc_nitrogen_state_set_totvegcthresh!(ns::SoilBiogeochemNitrogenStateData,
                                                      totvegcthresh::Real)
    if totvegcthresh <= 0.0
        error("totvegcthresh is zero or negative and should be > 0")
    end
    ns.totvegcthresh = totvegcthresh
    return nothing
end

"""
    soil_bgc_nitrogen_state_init_history!(ns, bounds_col)

Stub for history field registration (no-op in Julia port).
Corresponds to `InitHistory` in the Fortran source.
"""
function soil_bgc_nitrogen_state_init_history!(ns::SoilBiogeochemNitrogenStateData,
                                                bounds_col::UnitRange{Int})
    return nothing
end

"""
    soil_bgc_nitrogen_state_restart!(ns, bounds_col)

Stub for restart read/write (no-op in Julia port).
Corresponds to `Restart` in the Fortran source.
"""
function soil_bgc_nitrogen_state_restart!(ns::SoilBiogeochemNitrogenStateData,
                                           bounds_col::UnitRange{Int})
    return nothing
end

"""
    soil_bgc_nitrogen_state_dynamic_col_adjustments!(ns, bounds_col)

Stub for dynamic column area adjustments (no-op in Julia port).
Corresponds to `DynamicColumnAdjustments` in the Fortran source.
"""
function soil_bgc_nitrogen_state_dynamic_col_adjustments!(ns::SoilBiogeochemNitrogenStateData,
                                                            bounds_col::UnitRange{Int})
    return nothing
end
