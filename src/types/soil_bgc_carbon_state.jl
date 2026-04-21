# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemCarbonStateType.F90
# Soil biogeochemistry carbon state data type allocation and initialization
# ==========================================================================

"""
    SoilBiogeochemCarbonStateData

Soil biogeochemistry carbon state data structure. Holds decomposition C pools
at column and gridcell levels for C12 (and optionally C13/C14 isotopes).

Ported from `soilbiogeochem_carbonstate_type` in `SoilBiogeochemCarbonStateType.F90`.
"""
Base.@kwdef mutable struct SoilBiogeochemCarbonStateData{FT<:Real}
    # --- Vertically-resolved decomposition pools (col × nlev × npools) ---
    decomp_cpools_vr_col              ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3)
    decomp0_cpools_vr_col             ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3) baseline (initial value of this year)
    decomp_cpools_vr_SASUsave_col     ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3) SASU save
    decomp_soilc_vr_col               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gC/m3) total soil C vr
    ctrunc_vr_col                     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gC/m3) C truncation vr

    # --- Summary (diagnostic) state variables ---
    ctrunc_col                        ::Vector{FT} = Float64[]  # (gC/m2) C truncation
    totmicc_col                       ::Vector{FT} = Float64[]  # (gC/m2) total microbial C
    totlitc_col                       ::Vector{FT} = Float64[]  # (gC/m2) total litter C
    totlitc_1m_col                    ::Vector{FT} = Float64[]  # (gC/m2) total litter C to 1m
    totsomc_col                       ::Vector{FT} = Float64[]  # (gC/m2) total SOM C
    totsomc_1m_col                    ::Vector{FT} = Float64[]  # (gC/m2) total SOM C to 1m
    cwdc_col                          ::Vector{FT} = Float64[]  # (gC/m2) coarse woody debris C
    decomp_cpools_1m_col              ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # (gC/m2) decomp pools to 1m
    decomp_cpools_col                 ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # (gC/m2) decomp pools
    dyn_cbal_adjustments_col          ::Vector{FT} = Float64[]  # (gC/m2) dynamic column area adjustments

    # --- Scalars ---
    restart_file_spinup_state         ::Int = typemax(Int)  # spinup state from restart file
    totvegcthresh                     ::Float64 = NaN       # threshold for zeroing decomp pools

    # --- Carbon totals (includes soil, cpool, veg) ---
    totc_col                          ::Vector{FT} = Float64[]  # (gC/m2) total column C
    totecosysc_col                    ::Vector{FT} = Float64[]  # (gC/m2) total ecosystem C
    totc_grc                          ::Vector{FT} = Float64[]  # (gC/m2) total gridcell C

    # --- Matrix-CN fields ---
    matrix_cap_decomp_cpools_col      ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gC/m2)
    matrix_cap_decomp_cpools_vr_col   ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3)
    in_acc                            ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gC/m3/yr)
    in_acc_2d                         ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3/yr)
    tran_acc                          ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3/yr)
    vert_up_tran_acc                  ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3/yr)
    vert_down_tran_acc                ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3/yr)
    exit_acc                          ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3/yr)
    hori_tran_acc                     ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3/yr)
end

# ---------------------------------------------------------------------------
# Helper constructors
# ---------------------------------------------------------------------------
nanvec(n) = fill(NaN, n)
nanmat(r, c) = fill(NaN, r, c)
nan3d(a, b, c) = fill(NaN, a, b, c)

"""
    soil_bgc_carbon_state_init!(cs, nc, ng, nlevdecomp_full, ndecomp_pools;
                                nlevdecomp=nlevdecomp_full,
                                ndecomp_cascade_transitions=0,
                                use_soil_matrixcn=false)

Allocate all fields of `SoilBiogeochemCarbonStateData`.
Corresponds to `InitAllocate` in the Fortran source.
"""
function soil_bgc_carbon_state_init!(cs::SoilBiogeochemCarbonStateData,
                                     nc::Int, ng::Int,
                                     nlevdecomp_full::Int,
                                     ndecomp_pools::Int;
                                     nlevdecomp::Int=nlevdecomp_full,
                                     ndecomp_cascade_transitions::Int=0,
                                     use_soil_matrixcn::Bool=false)

    cs.totvegcthresh = NaN

    # --- Column-level 2D (col × npools) ---
    cs.decomp_cpools_col    = nanmat(nc, ndecomp_pools)
    cs.decomp_cpools_1m_col = nanmat(nc, ndecomp_pools)

    if use_soil_matrixcn
        cs.matrix_cap_decomp_cpools_col = nanmat(nc, ndecomp_pools)
    end

    # --- Vertically-resolved 2D (col × nlev) ---
    cs.ctrunc_vr_col       = nanmat(nc, nlevdecomp_full)
    cs.decomp_soilc_vr_col = nanmat(nc, nlevdecomp_full)

    # --- Vertically-resolved 3D (col × nlev × npools) ---
    cs.decomp_cpools_vr_col = nan3d(nc, nlevdecomp_full, ndecomp_pools)

    if use_soil_matrixcn
        cs.matrix_cap_decomp_cpools_vr_col = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        cs.decomp0_cpools_vr_col           = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        cs.decomp_cpools_vr_SASUsave_col   = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        cs.in_acc      = nanmat(nc, nlevdecomp * ndecomp_pools)
        cs.tran_acc    = nan3d(nc, nlevdecomp * ndecomp_pools, nlevdecomp * ndecomp_pools)
        cs.in_acc_2d   = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        cs.vert_up_tran_acc   = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        cs.vert_down_tran_acc = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        cs.exit_acc           = nan3d(nc, nlevdecomp_full, ndecomp_pools)
        cs.hori_tran_acc      = nan3d(nc, nlevdecomp_full, ndecomp_cascade_transitions)
    end

    # --- Column-level 1D ---
    cs.ctrunc_col     = nanvec(nc)
    cs.cwdc_col       = nanvec(nc)
    cs.totmicc_col    = nanvec(nc)
    cs.totlitc_col    = nanvec(nc)
    cs.totsomc_col    = nanvec(nc)
    cs.totlitc_1m_col = nanvec(nc)
    cs.totsomc_1m_col = nanvec(nc)
    cs.dyn_cbal_adjustments_col = nanvec(nc)

    cs.totc_col       = nanvec(nc)
    cs.totecosysc_col = nanvec(nc)
    cs.totc_grc       = nanvec(ng)

    cs.restart_file_spinup_state = typemax(Int)

    return nothing
end

"""
    soil_bgc_carbon_state_set_values!(cs, mask_col, value_column;
                                      nlevdecomp_full=..., ndecomp_pools=...,
                                      use_soil_matrixcn=false, ...)

Set carbon state variables for columns matching `mask_col`.
Corresponds to `SetValues` in the Fortran source.
"""
function soil_bgc_carbon_state_set_values!(cs::SoilBiogeochemCarbonStateData,
                                           mask_col::BitVector,
                                           value_column::Real;
                                           nlevdecomp_full::Int=size(cs.ctrunc_vr_col, 2),
                                           nlevdecomp::Int=nlevdecomp_full,
                                           ndecomp_pools::Int=size(cs.decomp_cpools_col, 2),
                                           ndecomp_cascade_transitions::Int=0,
                                           use_soil_matrixcn::Bool=false,
                                           is_fates::Union{BitVector,Nothing}=nothing)
    for i in eachindex(mask_col)
        mask_col[i] || continue
        if is_fates === nothing || !is_fates[i]
            cs.cwdc_col[i] = value_column
        end
        cs.ctrunc_col[i]     = value_column
        cs.totmicc_col[i]    = value_column
        cs.totlitc_col[i]    = value_column
        cs.totlitc_1m_col[i] = value_column
        cs.totsomc_col[i]    = value_column
        cs.totsomc_1m_col[i] = value_column
        cs.totc_col[i]       = value_column
        cs.totecosysc_col[i] = value_column
    end

    for j in 1:nlevdecomp_full
        for i in eachindex(mask_col)
            mask_col[i] || continue
            cs.ctrunc_vr_col[i, j] = value_column
        end
    end

    for k in 1:ndecomp_pools
        for i in eachindex(mask_col)
            mask_col[i] || continue
            cs.decomp_cpools_col[i, k]    = value_column
            cs.decomp_cpools_1m_col[i, k] = value_column
            if use_soil_matrixcn
                cs.matrix_cap_decomp_cpools_col[i, k] = value_column
            end
        end
    end

    for j in 1:nlevdecomp_full
        for k in 1:ndecomp_pools
            for i in eachindex(mask_col)
                mask_col[i] || continue
                cs.decomp_cpools_vr_col[i, j, k] = value_column
                if use_soil_matrixcn
                    cs.matrix_cap_decomp_cpools_vr_col[i, j, k] = value_column
                    cs.decomp0_cpools_vr_col[i, j, k]           = value_column
                end
            end
        end
    end

    if use_soil_matrixcn
        for j in 1:nlevdecomp
            for k in 1:ndecomp_pools
                for i in eachindex(mask_col)
                    mask_col[i] || continue
                    cs.in_acc_2d[i, j, k]          = value_column
                    cs.vert_up_tran_acc[i, j, k]   = value_column
                    cs.vert_down_tran_acc[i, j, k] = value_column
                    cs.exit_acc[i, j, k]           = value_column
                end
            end
            for k in 1:ndecomp_cascade_transitions
                for i in eachindex(mask_col)
                    mask_col[i] || continue
                    cs.hori_tran_acc[i, j, k] = value_column
                end
            end
        end
    end

    return nothing
end

"""
    soil_bgc_carbon_state_init_cold!(cs, bounds_col, ratio;
                                     nlevdecomp=..., nlevdecomp_full=..., ndecomp_pools=...,
                                     initial_stock=nothing, initial_stock_soildepth=0.3,
                                     zsoi=nothing, use_soil_matrixcn=false, ...)

Cold-start initialization of soil C pools.
Corresponds to `InitCold` in the Fortran source.
"""
function soil_bgc_carbon_state_init_cold!(cs::SoilBiogeochemCarbonStateData,
                                          bounds_col::UnitRange{Int},
                                          ratio::Real;
                                          nlevdecomp::Int=size(cs.decomp_cpools_vr_col, 2),
                                          nlevdecomp_full::Int=nlevdecomp,
                                          ndecomp_pools::Int=size(cs.decomp_cpools_vr_col, 3),
                                          ndecomp_cascade_transitions::Int=0,
                                          initial_stock::Vector{<:Real}=zeros(ndecomp_pools),
                                          initial_stock_soildepth::Real=0.3,
                                          zsoi_vals::Vector{<:Real}=zeros(nlevdecomp),
                                          use_soil_matrixcn::Bool=false,
                                          mask_soil_crop::Union{BitVector,Nothing}=nothing,
                                          c12_inst::Union{SoilBiogeochemCarbonStateData,Nothing}=nothing)
    for c in bounds_col
        # matrix spinup initialization
        if use_soil_matrixcn
            cs.in_acc[c, :] .= 0.0
        end

        is_soil_crop = mask_soil_crop === nothing || mask_soil_crop[c]

        if is_soil_crop
            if c12_inst === nothing  # c12 initialization
                for j in 1:nlevdecomp
                    for k in 1:ndecomp_pools
                        if zsoi_vals[j] < initial_stock_soildepth
                            cs.decomp_cpools_vr_col[c, j, k] = initial_stock[k]
                            if use_soil_matrixcn
                                cs.matrix_cap_decomp_cpools_vr_col[c, j, k] = initial_stock[k]
                            end
                        else
                            cs.decomp_cpools_vr_col[c, j, k] = 0.0
                            if use_soil_matrixcn
                                cs.matrix_cap_decomp_cpools_vr_col[c, j, k] = 0.0
                            end
                        end
                    end
                    cs.ctrunc_vr_col[c, j] = 0.0
                end
                if nlevdecomp > 1
                    for j in (nlevdecomp+1):nlevdecomp_full
                        for k in 1:ndecomp_pools
                            cs.decomp_cpools_vr_col[c, j, k] = 0.0
                            if use_soil_matrixcn
                                cs.matrix_cap_decomp_cpools_vr_col[c, j, k] = 0.0
                            end
                        end
                        cs.ctrunc_vr_col[c, j] = 0.0
                    end
                end
                cs.decomp_cpools_col[c, 1:ndecomp_pools]    .= initial_stock[1:ndecomp_pools]
                cs.decomp_cpools_1m_col[c, 1:ndecomp_pools] .= initial_stock[1:ndecomp_pools]
                if use_soil_matrixcn
                    cs.matrix_cap_decomp_cpools_col[c, 1:ndecomp_pools] .= initial_stock[1:ndecomp_pools]
                end
            else  # c13/c14 initialization from c12 instance
                for j in 1:nlevdecomp
                    for k in 1:ndecomp_pools
                        cs.decomp_cpools_vr_col[c, j, k] = c12_inst.decomp_cpools_vr_col[c, j, k] * ratio
                        if use_soil_matrixcn
                            cs.matrix_cap_decomp_cpools_vr_col[c, j, k] = c12_inst.matrix_cap_decomp_cpools_vr_col[c, j, k] * ratio
                        end
                    end
                    cs.ctrunc_vr_col[c, j] = c12_inst.ctrunc_vr_col[c, j] * ratio
                end
                if nlevdecomp > 1
                    for j in (nlevdecomp+1):nlevdecomp_full
                        for k in 1:ndecomp_pools
                            cs.decomp_cpools_vr_col[c, j, k] = 0.0
                            if use_soil_matrixcn
                                cs.matrix_cap_decomp_cpools_vr_col[c, j, k] = 0.0
                            end
                        end
                        cs.ctrunc_vr_col[c, j] = 0.0
                    end
                end
                for k in 1:ndecomp_pools
                    cs.decomp_cpools_col[c, k]    = c12_inst.decomp_cpools_col[c, k] * ratio
                    cs.decomp_cpools_1m_col[c, k] = c12_inst.decomp_cpools_1m_col[c, k] * ratio
                    if use_soil_matrixcn
                        cs.matrix_cap_decomp_cpools_col[c, k] = c12_inst.matrix_cap_decomp_cpools_col[c, k] * ratio
                    end
                end
            end

            # Matrix-CN accumulated fields
            if use_soil_matrixcn
                for j in 1:nlevdecomp_full
                    for k in 1:ndecomp_pools
                        cs.in_acc_2d[c, j, k]          = 0.0
                        cs.vert_up_tran_acc[c, j, k]   = 0.0
                        cs.vert_down_tran_acc[c, j, k] = 0.0
                        cs.exit_acc[c, j, k]           = 0.0
                        cs.decomp0_cpools_vr_col[c, j, k]       = max(cs.decomp_cpools_vr_col[c, j, k], 1.0e-30)
                        cs.decomp_cpools_vr_SASUsave_col[c, j, k] = 0.0
                    end
                    for k in 1:ndecomp_cascade_transitions
                        cs.hori_tran_acc[c, j, k] = 0.0
                    end
                end
            end

            # Summary diagnostic fields
            cs.ctrunc_col[c]     = 0.0
            cs.totmicc_col[c]    = 0.0
            cs.totlitc_col[c]    = 0.0
            cs.totsomc_col[c]    = 0.0
            cs.totlitc_1m_col[c] = 0.0
            cs.totsomc_1m_col[c] = 0.0
            cs.totc_col[c]       = 0.0
            cs.totecosysc_col[c] = 0.0

            # cwdc: for isotopes use c12 ratio, for c12 set to zero
            if c12_inst !== nothing
                cs.cwdc_col[c] = c12_inst.cwdc_col[c] * ratio
            else
                cs.cwdc_col[c] = 0.0
            end
        end
    end

    return nothing
end

"""
    soil_bgc_carbon_state_summary!(cs, mask_allc, bounds_col;
                                   nlevdecomp, nlevdecomp_full, ndecomp_pools,
                                   dzsoi_decomp_vals, zisoi_vals,
                                   is_litter, is_soil, is_microbe, is_cwd,
                                   totc_p2c_col, totvegc_col, is_fates_col,
                                   use_soil_matrixcn, use_fates_bgc, ...)

Perform column-level carbon summary calculations.
Corresponds to `Summary` in the Fortran source.
"""
function soil_bgc_carbon_state_summary!(cs::SoilBiogeochemCarbonStateData,
                                        mask_allc::BitVector,
                                        bounds_col::UnitRange{Int};
                                        nlevdecomp::Int,
                                        nlevdecomp_full::Int=nlevdecomp,
                                        ndecomp_pools::Int,
                                        dzsoi_decomp_vals::Vector{<:Real},
                                        zisoi_vals::Vector{<:Real}=Float64[],
                                        is_litter::BitVector=falses(ndecomp_pools),
                                        is_soil::BitVector=falses(ndecomp_pools),
                                        is_microbe::BitVector=falses(ndecomp_pools),
                                        is_cwd::BitVector=falses(ndecomp_pools),
                                        totc_p2c_col::Vector{<:Real}=zeros(length(mask_allc)),
                                        totvegc_col::Vector{<:Real}=zeros(length(mask_allc)),
                                        is_fates_col::Union{BitVector,Nothing}=nothing,
                                        use_soil_matrixcn::Bool=false,
                                        use_fates_bgc::Bool=false,
                                        mask_bgc_soilc::Union{BitVector,Nothing}=nothing)

    # Vertically integrate each decomposing C pool
    for l in 1:ndecomp_pools
        for c in bounds_col
            mask_allc[c] || continue
            cs.decomp_cpools_col[c, l] = 0.0
            if use_soil_matrixcn
                cs.matrix_cap_decomp_cpools_col[c, l] = 0.0
            end
        end
    end
    for l in 1:ndecomp_pools
        for j in 1:nlevdecomp
            for c in bounds_col
                mask_allc[c] || continue
                cs.decomp_cpools_col[c, l] += cs.decomp_cpools_vr_col[c, j, l] * dzsoi_decomp_vals[j]
                if use_soil_matrixcn
                    cs.matrix_cap_decomp_cpools_col[c, l] += cs.matrix_cap_decomp_cpools_vr_col[c, j, l] * dzsoi_decomp_vals[j]
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
                cs.decomp_cpools_1m_col[c, l] = 0.0
            end
        end
        for l in 1:ndecomp_pools
            for j in 1:nlevdecomp
                if zisoi_vals[j+1] <= maxdepth
                    for c in bounds_col
                        mask_allc[c] || continue
                        cs.decomp_cpools_1m_col[c, l] += cs.decomp_cpools_vr_col[c, j, l] * dzsoi_decomp_vals[j]
                    end
                elseif zisoi_vals[j] < maxdepth
                    for c in bounds_col
                        mask_allc[c] || continue
                        cs.decomp_cpools_1m_col[c, l] += cs.decomp_cpools_vr_col[c, j, l] * (maxdepth - zisoi_vals[j])
                    end
                end
            end
        end
    end

    # Vertically-resolved decomposing total soil c pool
    if nlevdecomp_full > 1
        for j in 1:nlevdecomp
            for c in bounds_col
                mask_allc[c] || continue
                cs.decomp_soilc_vr_col[c, j] = 0.0
            end
        end
        for l in 1:ndecomp_pools
            if is_soil[l]
                for j in 1:nlevdecomp
                    for c in bounds_col
                        mask_allc[c] || continue
                        cs.decomp_soilc_vr_col[c, j] += cs.decomp_cpools_vr_col[c, j, l]
                    end
                end
            end
        end
    end

    # Truncation carbon
    for c in bounds_col
        mask_allc[c] || continue
        cs.ctrunc_col[c] = 0.0
    end
    for j in 1:nlevdecomp
        for c in bounds_col
            mask_allc[c] || continue
            cs.ctrunc_col[c] += cs.ctrunc_vr_col[c, j] * dzsoi_decomp_vals[j]
        end
    end

    # Total litter carbon to 1m
    if nlevdecomp > 1
        for c in bounds_col
            mask_allc[c] || continue
            cs.totlitc_1m_col[c] = 0.0
        end
        for l in 1:ndecomp_pools
            if is_litter[l]
                for c in bounds_col
                    mask_allc[c] || continue
                    cs.totlitc_1m_col[c] += cs.decomp_cpools_1m_col[c, l]
                end
            end
        end
    end

    # Total SOM carbon to 1m
    if nlevdecomp > 1
        for c in bounds_col
            mask_allc[c] || continue
            cs.totsomc_1m_col[c] = 0.0
        end
        for l in 1:ndecomp_pools
            if is_soil[l]
                for c in bounds_col
                    mask_allc[c] || continue
                    cs.totsomc_1m_col[c] += cs.decomp_cpools_1m_col[c, l]
                end
            end
        end
    end

    # Total microbial carbon
    for c in bounds_col
        mask_allc[c] || continue
        cs.totmicc_col[c] = 0.0
    end
    for l in 1:ndecomp_pools
        if is_microbe[l]
            for c in bounds_col
                mask_allc[c] || continue
                cs.totmicc_col[c] += cs.decomp_cpools_col[c, l]
            end
        end
    end

    # Total litter carbon
    for c in bounds_col
        mask_allc[c] || continue
        cs.totlitc_col[c] = 0.0
    end
    for l in 1:ndecomp_pools
        if is_litter[l]
            for c in bounds_col
                mask_allc[c] || continue
                cs.totlitc_col[c] += cs.decomp_cpools_col[c, l]
            end
        end
    end

    # Total SOM carbon
    for c in bounds_col
        mask_allc[c] || continue
        cs.totsomc_col[c] = 0.0
    end
    for l in 1:ndecomp_pools
        if is_soil[l]
            for c in bounds_col
                mask_allc[c] || continue
                cs.totsomc_col[c] += cs.decomp_cpools_col[c, l]
            end
        end
    end

    # CWD, ecosystem C, total column C
    for c in bounds_col
        mask_allc[c] || continue
        cs.cwdc_col[c] = 0.0
    end

    for c in bounds_col
        mask_allc[c] || continue

        is_fates_c = is_fates_col !== nothing && is_fates_col[c]

        local ecovegc, tvegc
        if is_fates_c
            tvegc   = 0.0
            ecovegc = 0.0
        else
            for l in 1:ndecomp_pools
                if is_cwd[l]
                    cs.cwdc_col[c] += cs.decomp_cpools_col[c, l]
                end
            end
            tvegc   = totc_p2c_col[c]
            ecovegc = totvegc_col[c]
        end

        # total ecosystem carbon (TOTECOSYSC)
        cs.totecosysc_col[c] = cs.cwdc_col[c] + cs.totmicc_col[c] +
                               cs.totlitc_col[c] + cs.totsomc_col[c] + ecovegc

        # total column carbon (TOTCOLC)
        cs.totc_col[c] = cs.cwdc_col[c] + cs.totmicc_col[c] +
                         cs.totlitc_col[c] + cs.totsomc_col[c] +
                         cs.ctrunc_col[c] + tvegc
    end

    return nothing
end

"""
    soil_bgc_carbon_state_set_totvegcthresh!(cs, totvegcthresh)

Set the total vegetation carbon threshold for spinup.
Corresponds to `SetTotVgCThresh` in the Fortran source.
"""
function soil_bgc_carbon_state_set_totvegcthresh!(cs::SoilBiogeochemCarbonStateData,
                                                   totvegcthresh::Real)
    if totvegcthresh <= 0.0
        error("totvegcthresh is zero or negative and should be > 0")
    end
    cs.totvegcthresh = totvegcthresh
    return nothing
end

"""
    soil_bgc_carbon_state_init_history!(cs, bounds_col)

Stub for history field registration (no-op in Julia port).
Corresponds to `InitHistory` in the Fortran source.
"""
function soil_bgc_carbon_state_init_history!(cs::SoilBiogeochemCarbonStateData,
                                             bounds_col::UnitRange{Int})
    return nothing
end

"""
    soil_bgc_carbon_state_restart!(cs, bounds_col)

Stub for restart read/write (no-op in Julia port).
Corresponds to `Restart` in the Fortran source.
"""
function soil_bgc_carbon_state_restart!(cs::SoilBiogeochemCarbonStateData,
                                        bounds_col::UnitRange{Int})
    return nothing
end

"""
    soil_bgc_carbon_state_dynamic_col_adjustments!(cs, bounds_col)

Stub for dynamic column area adjustments (no-op in Julia port).
Corresponds to `DynamicColumnAdjustments` in the Fortran source.
"""
function soil_bgc_carbon_state_dynamic_col_adjustments!(cs::SoilBiogeochemCarbonStateData,
                                                         bounds_col::UnitRange{Int})
    return nothing
end
