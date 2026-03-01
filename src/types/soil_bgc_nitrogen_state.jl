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
Base.@kwdef mutable struct SoilBiogeochemNitrogenStateData
    # --- Vertically-resolved decomposition pools (col × nlev × npools) ---
    decomp_npools_vr_col              ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3)
    decomp0_npools_vr_col             ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3) baseline (initial value of this year)
    decomp_npools_vr_SASUsave_col     ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3) SASU save
    decomp_soiln_vr_col               ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gN/m3) total soil N vr

    # --- Vertically-resolved mineral N and truncation ---
    sminn_vr_col                      ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gN/m3) soil mineral N vr
    ntrunc_vr_col                     ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gN/m3) N truncation vr

    # --- Nitrification/denitrification pools ---
    smin_no3_vr_col                   ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gN/m3) soil mineral NO3 vr
    smin_no3_col                      ::Vector{Float64}  = Float64[]                        # (gN/m2) soil mineral NO3
    smin_nh4_vr_col                   ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gN/m3) soil mineral NH4 vr
    smin_nh4_col                      ::Vector{Float64}  = Float64[]                        # (gN/m2) soil mineral NH4

    # --- Summary (diagnostic) state variables ---
    decomp_npools_col                 ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gN/m2) decomp N pools
    decomp_npools_1m_col              ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gN/m2) decomp N pools to 1m
    sminn_col                         ::Vector{Float64}  = Float64[]                        # (gN/m2) soil mineral N
    ntrunc_col                        ::Vector{Float64}  = Float64[]                        # (gN/m2) N truncation
    cwdn_col                          ::Vector{Float64}  = Float64[]                        # (gN/m2) coarse woody debris N
    totlitn_col                       ::Vector{Float64}  = Float64[]                        # (gN/m2) total litter N
    totmicn_col                       ::Vector{Float64}  = Float64[]                        # (gN/m2) total microbial N
    totsomn_col                       ::Vector{Float64}  = Float64[]                        # (gN/m2) total SOM N
    totlitn_1m_col                    ::Vector{Float64}  = Float64[]                        # (gN/m2) total litter N to 1m
    totsomn_1m_col                    ::Vector{Float64}  = Float64[]                        # (gN/m2) total SOM N to 1m
    dyn_nbal_adjustments_col          ::Vector{Float64}  = Float64[]                        # (gN/m2) dynamic column N adjustments
    dyn_no3bal_adjustments_col        ::Vector{Float64}  = Float64[]                        # (gN/m2) dynamic column NO3 adjustments
    dyn_nh4bal_adjustments_col        ::Vector{Float64}  = Float64[]                        # (gN/m2) dynamic column NH4 adjustments

    # --- Scalars ---
    totvegcthresh                     ::Float64 = NaN       # threshold for zeroing decomp pools

    # --- Nitrogen totals (includes soil, veg) ---
    totn_col                          ::Vector{Float64}  = Float64[]  # (gN/m2) total column N
    totecosysn_col                    ::Vector{Float64}  = Float64[]  # (gN/m2) total ecosystem N
    totn_grc                          ::Vector{Float64}  = Float64[]  # (gN/m2) total gridcell N

    # --- Matrix-CN fields ---
    matrix_cap_decomp_npools_col      ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gN/m2)
    matrix_cap_decomp_npools_vr_col   ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3)
    in_nacc                           ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gN/m3/yr)
    in_nacc_2d                        ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/yr)
    tran_nacc                         ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/yr)
    vert_up_tran_nacc                 ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/yr)
    vert_down_tran_nacc               ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/yr)
    exit_nacc                         ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/yr)
    hori_tran_nacc                    ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/yr)
end

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
                                              value_column::Float64;
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
                                             decomp_cpools_vr_col::Array{Float64,3}=zeros(0,0,0),
                                             decomp_cpools_col::Matrix{Float64}=zeros(0,0),
                                             decomp_cpools_1m_col::Matrix{Float64}=zeros(0,0),
                                             initial_cn_ratio::Vector{Float64}=ones(ndecomp_pools),
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
                                           mask_allc::BitVector,
                                           bounds_col::UnitRange{Int};
                                           nlevdecomp::Int,
                                           nlevdecomp_full::Int=nlevdecomp,
                                           ndecomp_pools::Int,
                                           dzsoi_decomp_vals::Vector{Float64},
                                           zisoi_vals::Vector{Float64}=Float64[],
                                           is_litter::BitVector=falses(ndecomp_pools),
                                           is_soil::BitVector=falses(ndecomp_pools),
                                           is_microbe::BitVector=falses(ndecomp_pools),
                                           is_cwd::BitVector=falses(ndecomp_pools),
                                           totn_p2c_col::Vector{Float64}=zeros(length(mask_allc)),
                                           totvegn_col::Vector{Float64}=zeros(length(mask_allc)),
                                           is_fates_col::Union{BitVector,Nothing}=nothing,
                                           use_soil_matrixcn::Bool=false,
                                           use_fates_bgc::Bool=false,
                                           use_nitrif_denitrif::Bool=false,
                                           mask_bgc_soilc::Union{BitVector,Nothing}=nothing)

    # Vertically integrate NO3/NH4 pools
    if use_nitrif_denitrif
        for c in bounds_col
            mask_allc[c] || continue
            ns.smin_no3_col[c] = 0.0
            ns.smin_nh4_col[c] = 0.0
        end
        for j in 1:nlevdecomp
            for c in bounds_col
                mask_allc[c] || continue
                ns.smin_no3_col[c] += ns.smin_no3_vr_col[c, j] * dzsoi_decomp_vals[j]
                ns.smin_nh4_col[c] += ns.smin_nh4_vr_col[c, j] * dzsoi_decomp_vals[j]
            end
        end
    end

    # Vertically integrate each decomposing N pool
    for l in 1:ndecomp_pools
        for c in bounds_col
            mask_allc[c] || continue
            ns.decomp_npools_col[c, l] = 0.0
            if use_soil_matrixcn
                ns.matrix_cap_decomp_npools_col[c, l] = 0.0
            end
        end
    end
    for l in 1:ndecomp_pools
        for j in 1:nlevdecomp
            for c in bounds_col
                mask_allc[c] || continue
                ns.decomp_npools_col[c, l] += ns.decomp_npools_vr_col[c, j, l] * dzsoi_decomp_vals[j]
                if use_soil_matrixcn
                    ns.matrix_cap_decomp_npools_col[c, l] += ns.matrix_cap_decomp_npools_vr_col[c, j, l] * dzsoi_decomp_vals[j]
                end
            end
        end
    end

    # Vertically integrate to 1 meter
    if nlevdecomp > 1
        maxdepth = 1.0
        for l in 1:ndecomp_pools
            for c in bounds_col
                mask_allc[c] || continue
                ns.decomp_npools_1m_col[c, l] = 0.0
            end
        end
        for l in 1:ndecomp_pools
            for j in 1:nlevdecomp
                if zisoi_vals[j+1] <= maxdepth
                    for c in bounds_col
                        mask_allc[c] || continue
                        ns.decomp_npools_1m_col[c, l] += ns.decomp_npools_vr_col[c, j, l] * dzsoi_decomp_vals[j]
                    end
                elseif zisoi_vals[j] < maxdepth
                    for c in bounds_col
                        mask_allc[c] || continue
                        ns.decomp_npools_1m_col[c, l] += ns.decomp_npools_vr_col[c, j, l] * (maxdepth - zisoi_vals[j])
                    end
                end
            end
        end

        # Vertically-resolved decomposing total soil N pool
        if nlevdecomp_full > 1
            for j in 1:nlevdecomp
                for c in bounds_col
                    mask_allc[c] || continue
                    ns.decomp_soiln_vr_col[c, j] = 0.0
                end
            end
            for l in 1:ndecomp_pools
                if is_soil[l]
                    for j in 1:nlevdecomp
                        for c in bounds_col
                            mask_allc[c] || continue
                            ns.decomp_soiln_vr_col[c, j] += ns.decomp_npools_vr_col[c, j, l]
                        end
                    end
                end
            end
        end

        # Total litter nitrogen to 1m
        for c in bounds_col
            mask_allc[c] || continue
            ns.totlitn_1m_col[c] = 0.0
        end
        for l in 1:ndecomp_pools
            if is_litter[l]
                for c in bounds_col
                    mask_allc[c] || continue
                    ns.totlitn_1m_col[c] += ns.decomp_npools_1m_col[c, l]
                end
            end
        end

        # Total SOM nitrogen to 1m
        for c in bounds_col
            mask_allc[c] || continue
            ns.totsomn_1m_col[c] = 0.0
        end
        for l in 1:ndecomp_pools
            if is_soil[l]
                for c in bounds_col
                    mask_allc[c] || continue
                    ns.totsomn_1m_col[c] += ns.decomp_npools_1m_col[c, l]
                end
            end
        end
    end

    # Total litter nitrogen
    for c in bounds_col
        mask_allc[c] || continue
        ns.totlitn_col[c] = 0.0
    end
    for l in 1:ndecomp_pools
        if is_litter[l]
            for c in bounds_col
                mask_allc[c] || continue
                ns.totlitn_col[c] += ns.decomp_npools_col[c, l]
            end
        end
    end

    # Total microbial nitrogen
    for c in bounds_col
        mask_allc[c] || continue
        ns.totmicn_col[c] = 0.0
    end
    for l in 1:ndecomp_pools
        if is_microbe[l]
            for c in bounds_col
                mask_allc[c] || continue
                ns.totmicn_col[c] += ns.decomp_npools_col[c, l]
            end
        end
    end

    # Total SOM nitrogen
    for c in bounds_col
        mask_allc[c] || continue
        ns.totsomn_col[c] = 0.0
    end
    for l in 1:ndecomp_pools
        if is_soil[l]
            for c in bounds_col
                mask_allc[c] || continue
                ns.totsomn_col[c] += ns.decomp_npools_col[c, l]
            end
        end
    end

    # Total sminn
    for c in bounds_col
        mask_allc[c] || continue
        ns.sminn_col[c] = 0.0
    end
    for j in 1:nlevdecomp
        for c in bounds_col
            mask_allc[c] || continue
            ns.sminn_col[c] += ns.sminn_vr_col[c, j] * dzsoi_decomp_vals[j]
        end
    end

    # Total ntrunc
    for c in bounds_col
        mask_allc[c] || continue
        ns.ntrunc_col[c] = 0.0
    end
    for j in 1:nlevdecomp
        for c in bounds_col
            mask_allc[c] || continue
            ns.ntrunc_col[c] += ns.ntrunc_vr_col[c, j] * dzsoi_decomp_vals[j]
        end
    end

    # CWD nitrogen, ecosystem N, total column N
    for c in bounds_col
        mask_allc[c] || continue
        ns.cwdn_col[c] = 0.0
    end

    for c in bounds_col
        mask_allc[c] || continue

        is_fates_c = is_fates_col !== nothing && is_fates_col[c]

        local ecovegn, tvegn
        if is_fates_c
            tvegn   = 0.0
            ecovegn = 0.0
        else
            for l in 1:ndecomp_pools
                if is_cwd[l]
                    ns.cwdn_col[c] += ns.decomp_npools_col[c, l]
                end
            end
            tvegn   = totn_p2c_col[c]
            ecovegn = totvegn_col[c]
        end

        # total ecosystem nitrogen (TOTECOSYSN)
        ns.totecosysn_col[c] = ns.cwdn_col[c] + ns.totlitn_col[c] +
                                ns.totmicn_col[c] + ns.totsomn_col[c] +
                                ns.sminn_col[c] + ecovegn

        # total column nitrogen (TOTCOLN)
        ns.totn_col[c] = ns.cwdn_col[c] + ns.totlitn_col[c] +
                          ns.totmicn_col[c] + ns.totsomn_col[c] +
                          ns.sminn_col[c] + ns.ntrunc_col[c] + tvegn
    end

    return nothing
end

"""
    soil_bgc_nitrogen_state_set_totvegcthresh!(ns, totvegcthresh)

Set the total vegetation carbon threshold for spinup.
Corresponds to `SetTotVgCThresh` in the Fortran source.
"""
function soil_bgc_nitrogen_state_set_totvegcthresh!(ns::SoilBiogeochemNitrogenStateData,
                                                      totvegcthresh::Float64)
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
