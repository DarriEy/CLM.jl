# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemCarbonFluxType.F90
# Soil biogeochemistry carbon flux data type allocation and initialization
# ==========================================================================

"""
    SoilBiogeochemCarbonFluxData

Soil biogeochemistry carbon flux data structure. Holds decomposition C fluxes,
heterotrophic respiration, vertical transport, and related diagnostic fields
at column level.

Ported from `soilbiogeochem_carbonflux_type` in `SoilBiogeochemCarbonFluxType.F90`.
"""
Base.@kwdef mutable struct SoilBiogeochemCarbonFluxData
    # --- Fire fluxes ---
    somc_fire_col                        ::Vector{Float64} = Float64[]  # (gC/m2/s) C emissions due to peat burning

    # --- Decomposition fluxes ---
    decomp_cpools_sourcesink_col         ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3/timestep) change in decomposing C pools
    c_overflow_vr                        ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3/s) C rejected by microbes
    decomp_cascade_hr_vr_col            ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3/s) het. resp. from decomposing C pools (vr)
    decomp_cascade_hr_col               ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gC/m2/s) het. resp. from decomposing C pools (integrated)
    decomp_cascade_ctransfer_vr_col     ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3/s) C transferred along decomp cascade (vr)
    decomp_cascade_ctransfer_col        ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gC/m2/s) C transferred along decomp cascade (integrated)
    cn_col                               ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gC/gN) C:N ratio by pool
    litr_lig_c_to_n_col                  ::Vector{Float64} = Float64[]  # (gC/gN) avg lignin C:N ratio
    rf_decomp_cascade_col               ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (frac) respired fraction in decomp step
    pathfrac_decomp_cascade_col         ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (frac) C fraction through given transition
    decomp_k_col                         ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (1/s) rate coefficient for decomposition

    # --- Soil matrix fields ---
    hr_vr_col                            ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gC/m3/s) total vr het. resp.
    o_scalar_col                         ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # fraction decomp limited by anoxia
    w_scalar_col                         ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # fraction decomp limited by moisture
    t_scalar_col                         ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # fraction decomp limited by temperature
    som_c_leached_col                    ::Vector{Float64} = Float64[]  # (gC/m2/s) total SOM C loss from vertical transport
    decomp_cpools_leached_col           ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gC/m2/s) C loss from vert transport per pool
    decomp_cpools_transport_tendency_col ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # (gC/m3/s) C tendency due to vert transport

    # --- Nitrif/denitrif ---
    phr_vr_col                           ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # (gC/m3/s) potential hr (not N-limited)
    fphr_col                             ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # fraction of potential het. resp.

    # --- Summary HR fluxes ---
    hr_col                               ::Vector{Float64} = Float64[]  # (gC/m2/s) total heterotrophic respiration
    michr_col                            ::Vector{Float64} = Float64[]  # (gC/m2/s) microbial het. resp.
    cwdhr_col                            ::Vector{Float64} = Float64[]  # (gC/m2/s) CWD het. resp.
    lithr_col                            ::Vector{Float64} = Float64[]  # (gC/m2/s) litter het. resp.
    somhr_col                            ::Vector{Float64} = Float64[]  # (gC/m2/s) SOM het. resp.
    soilc_change_col                     ::Vector{Float64} = Float64[]  # (gC/m2/s) FUN used soil C
    fates_litter_flux                    ::Vector{Float64} = Float64[]  # (gC/m2/s) total litter flux from FATES

    # --- Matrix-CN fields ---
    matrix_decomp_fire_k_col            ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # decomp rate due to fire
    tri_ma_vr                            ::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)    # vertical C transfer rate in sparse matrix format
    NE_AKallsoilc                        ::Int = 0
    RI_AKallsoilc                        ::Vector{Int} = Int[]
    CI_AKallsoilc                        ::Vector{Int} = Int[]
    RI_a                                 ::Vector{Int} = Int[]
    CI_a                                 ::Vector{Int} = Int[]
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
    soil_bgc_carbon_flux_init!(cf, nc, nlevdecomp_full, ndecomp_pools,
                                ndecomp_cascade_transitions;
                                nlevgrnd=nlevdecomp_full,
                                nlevdecomp=nlevdecomp_full,
                                ndecomp_cascade_outtransitions=0,
                                use_soil_matrixcn=false,
                                use_fates=false,
                                Ntri_setup=0)

Allocate all fields of `SoilBiogeochemCarbonFluxData`.
Corresponds to `InitAllocate` in the Fortran source.
"""
function soil_bgc_carbon_flux_init!(cf::SoilBiogeochemCarbonFluxData,
                                     nc::Int,
                                     nlevdecomp_full::Int,
                                     ndecomp_pools::Int,
                                     ndecomp_cascade_transitions::Int;
                                     nlevgrnd::Int=nlevdecomp_full,
                                     nlevdecomp::Int=nlevdecomp_full,
                                     ndecomp_cascade_outtransitions::Int=0,
                                     use_soil_matrixcn::Bool=false,
                                     use_fates::Bool=false,
                                     Ntri_setup::Int=0)

    # --- Environmental scalars (col × nlev) ---
    cf.t_scalar_col = nanmat(nc, nlevdecomp_full)
    cf.w_scalar_col = nanmat(nc, nlevdecomp_full)
    cf.o_scalar_col = nanmat(nc, nlevdecomp_full)

    # --- Nitrif/denitrif fields ---
    cf.phr_vr_col = nanmat(nc, nlevdecomp_full)
    cf.fphr_col   = nanmat(nc, nlevgrnd)

    # --- Column-level 1D ---
    cf.som_c_leached_col = nanvec(nc)
    cf.somc_fire_col     = nanvec(nc)

    # --- Vertically-resolved het. resp. ---
    cf.hr_vr_col = nanmat(nc, nlevdecomp_full)

    # --- 3D: col × nlev × npools ---
    cf.decomp_cpools_sourcesink_col = nan3d(nc, nlevdecomp_full, ndecomp_pools)

    # --- 3D: col × nlev × ntrans ---
    cf.c_overflow_vr                    = nan3d(nc, nlevdecomp_full, ndecomp_cascade_transitions)
    cf.decomp_cascade_hr_vr_col         = nan3d(nc, nlevdecomp_full, ndecomp_cascade_transitions)
    cf.decomp_cascade_ctransfer_vr_col  = nan3d(nc, nlevdecomp_full, ndecomp_cascade_transitions)
    cf.rf_decomp_cascade_col            = nan3d(nc, nlevdecomp_full, ndecomp_cascade_transitions)
    cf.pathfrac_decomp_cascade_col      = nan3d(nc, nlevdecomp_full, ndecomp_cascade_transitions)

    # --- 2D: col × ntrans ---
    cf.decomp_cascade_hr_col        = nanmat(nc, ndecomp_cascade_transitions)
    cf.decomp_cascade_ctransfer_col = nanmat(nc, ndecomp_cascade_transitions)

    # --- 2D: col × npools ---
    cf.cn_col                     = nanmat(nc, ndecomp_pools)
    cf.decomp_cpools_leached_col  = nanmat(nc, ndecomp_pools)

    # --- 3D: col × nlev × npools ---
    cf.decomp_k_col                          = nan3d(nc, nlevdecomp_full, ndecomp_pools)
    cf.decomp_cpools_transport_tendency_col   = nan3d(nc, nlevdecomp_full, ndecomp_pools)

    # --- Column-level 1D summary ---
    cf.hr_col            = nanvec(nc)
    cf.michr_col         = nanvec(nc)
    cf.cwdhr_col         = nanvec(nc)
    cf.lithr_col         = nanvec(nc)
    cf.somhr_col         = nanvec(nc)
    cf.soilc_change_col  = nanvec(nc)

    # --- FATES litter flux ---
    if use_fates
        cf.fates_litter_flux = nanvec(nc)
    else
        cf.fates_litter_flux = fill(NaN, 1)
    end

    # --- Matrix-CN fields ---
    if use_soil_matrixcn
        cf.matrix_decomp_fire_k_col = nanmat(nc, nlevdecomp * ndecomp_pools)
        cf.tri_ma_vr = fill(NaN, nc, Ntri_setup)

        Ntrans = (ndecomp_cascade_transitions - ndecomp_cascade_outtransitions) * nlevdecomp
        cf.NE_AKallsoilc = Ntrans + ndecomp_pools * nlevdecomp +
                            Ntri_setup + ndecomp_pools * nlevdecomp
        cf.RI_AKallsoilc = fill(-9999, cf.NE_AKallsoilc)
        cf.CI_AKallsoilc = fill(-9999, cf.NE_AKallsoilc)

        Ntrans_diag = (ndecomp_cascade_transitions - ndecomp_cascade_outtransitions) * nlevdecomp +
                       ndecomp_pools * nlevdecomp  # ndecomp_pools_vr
        cf.RI_a = fill(-9999, Ntrans_diag)
        cf.CI_a = fill(-9999, Ntrans_diag)
    else
        cf.tri_ma_vr = fill(NaN, 1, 1)
    end

    # --- Lignin C:N ratio ---
    cf.litr_lig_c_to_n_col = zeros(nc)

    return nothing
end

"""
    soil_bgc_carbon_flux_set_values!(cf, mask_col, value_column;
                                      nlevdecomp_full=..., nlevdecomp=...,
                                      ndecomp_pools=..., ndecomp_cascade_transitions=...,
                                      use_soil_matrixcn=false, Ntri_setup=0)

Set carbon flux variables for columns matching `mask_col`.
Corresponds to `SetValues` in the Fortran source.
"""
function soil_bgc_carbon_flux_set_values!(cf::SoilBiogeochemCarbonFluxData,
                                           mask_col::BitVector,
                                           value_column::Float64;
                                           nlevdecomp_full::Int=size(cf.hr_vr_col, 2),
                                           nlevdecomp::Int=nlevdecomp_full,
                                           ndecomp_pools::Int=size(cf.cn_col, 2),
                                           ndecomp_cascade_transitions::Int=size(cf.decomp_cascade_hr_col, 2),
                                           use_soil_matrixcn::Bool=false,
                                           Ntri_setup::Int=0)

    # Decomposition cascade transitions
    for l in 1:ndecomp_cascade_transitions
        for j in 1:nlevdecomp_full
            for i in eachindex(mask_col)
                mask_col[i] || continue
                cf.decomp_cascade_hr_col[i, l]             = value_column
                cf.c_overflow_vr[i, j, l]                   = value_column
                cf.decomp_cascade_hr_vr_col[i, j, l]        = value_column
                cf.decomp_cascade_ctransfer_col[i, l]       = value_column
                cf.decomp_cascade_ctransfer_vr_col[i, j, l] = value_column
                cf.pathfrac_decomp_cascade_col[i, j, l]     = value_column
                cf.rf_decomp_cascade_col[i, j, l]           = value_column
            end
        end
    end

    # Pool-indexed fields
    for k in 1:ndecomp_pools
        for i in eachindex(mask_col)
            mask_col[i] || continue
            cf.decomp_cpools_leached_col[i, k] = value_column
            cf.cn_col[i, k]                    = value_column
        end
        for j in 1:nlevdecomp_full
            for i in eachindex(mask_col)
                mask_col[i] || continue
                cf.decomp_cpools_transport_tendency_col[i, j, k] = value_column
                cf.decomp_cpools_sourcesink_col[i, j, k]         = value_column
                cf.decomp_k_col[i, j, k]                         = value_column
            end
        end
    end

    # Matrix-CN fields
    if use_soil_matrixcn
        for k in 1:ndecomp_pools
            for j in 1:nlevdecomp
                for i in eachindex(mask_col)
                    mask_col[i] || continue
                    cf.matrix_decomp_fire_k_col[i, j + nlevdecomp * (k - 1)] = value_column
                end
            end
        end
        for k in 1:Ntri_setup
            for i in eachindex(mask_col)
                mask_col[i] || continue
                cf.tri_ma_vr[i, k] = value_column
            end
        end
    end

    # Vertically-resolved HR
    for j in 1:nlevdecomp_full
        for i in eachindex(mask_col)
            mask_col[i] || continue
            cf.hr_vr_col[i, j] = value_column
        end
    end

    # Column-level 1D
    for i in eachindex(mask_col)
        mask_col[i] || continue
        cf.hr_col[i]            = value_column
        cf.somc_fire_col[i]     = value_column
        cf.som_c_leached_col[i] = value_column
        cf.somhr_col[i]         = value_column
        cf.lithr_col[i]         = value_column
        cf.cwdhr_col[i]         = value_column
        cf.michr_col[i]         = value_column
        cf.soilc_change_col[i]  = value_column
    end

    return nothing
end

"""
    soil_bgc_carbon_flux_init_cold!(cf, bounds_col; mask_special=nothing)

Cold-start initialization of soil C flux fields.
Corresponds to `InitCold` in the Fortran source.
Sets values to zero for special landunit columns.
"""
function soil_bgc_carbon_flux_init_cold!(cf::SoilBiogeochemCarbonFluxData,
                                          bounds_col::UnitRange{Int};
                                          mask_special::Union{BitVector,Nothing}=nothing)
    if mask_special !== nothing
        soil_bgc_carbon_flux_set_values!(cf, mask_special, 0.0)
    end
    return nothing
end

"""
    soil_bgc_carbon_flux_summary!(cf, mask_bgc_soilc, bounds_col;
                                   ndecomp_cascade_transitions, nlevdecomp, ndecomp_pools,
                                   dzsoi_decomp_vals, cascade_donor_pool,
                                   is_soil, is_litter, is_cwd, is_microbe)

Carbon flux summary calculations: vertically integrate HR and decomp cascade fluxes,
compute total HR by pool type (som, litter, cwd, microbial).
Corresponds to `Summary` in the Fortran source.
"""
function soil_bgc_carbon_flux_summary!(cf::SoilBiogeochemCarbonFluxData,
                                        mask_bgc_soilc::BitVector,
                                        bounds_col::UnitRange{Int};
                                        ndecomp_cascade_transitions::Int=size(cf.decomp_cascade_hr_col, 2),
                                        nlevdecomp::Int=size(cf.hr_vr_col, 2),
                                        ndecomp_pools::Int=size(cf.cn_col, 2),
                                        dzsoi_decomp_vals::Vector{Float64},
                                        cascade_donor_pool::Vector{Int}=ones(Int, ndecomp_cascade_transitions),
                                        is_soil::BitVector=falses(ndecomp_pools),
                                        is_litter::BitVector=falses(ndecomp_pools),
                                        is_cwd::BitVector=falses(ndecomp_pools),
                                        is_microbe::BitVector=falses(ndecomp_pools))

    # Initialize som_c_leached_col to zero
    for c in bounds_col
        mask_bgc_soilc[c] || continue
        cf.som_c_leached_col[c] = 0.0
    end

    # Vertically integrate HR and decomposition cascade transfer fluxes
    for k in 1:ndecomp_cascade_transitions
        for j in 1:nlevdecomp
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                cf.decomp_cascade_hr_col[c, k] =
                    cf.decomp_cascade_hr_col[c, k] +
                    cf.decomp_cascade_hr_vr_col[c, j, k] * dzsoi_decomp_vals[j]

                cf.decomp_cascade_ctransfer_col[c, k] =
                    cf.decomp_cascade_ctransfer_col[c, k] +
                    cf.decomp_cascade_ctransfer_vr_col[c, j, k] * dzsoi_decomp_vals[j]
            end
        end
    end

    # Total heterotrophic respiration, vertically resolved (HR)
    for j in 1:nlevdecomp
        for c in bounds_col
            mask_bgc_soilc[c] || continue
            cf.hr_vr_col[c, j] = 0.0
        end
    end
    for k in 1:ndecomp_cascade_transitions
        for j in 1:nlevdecomp
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                cf.hr_vr_col[c, j] =
                    cf.hr_vr_col[c, j] +
                    cf.decomp_cascade_hr_vr_col[c, j, k]
            end
        end
    end

    # Add up all vertical transport tendency terms and calculate total SOM leaching loss
    for l in 1:ndecomp_pools
        for c in bounds_col
            mask_bgc_soilc[c] || continue
            cf.decomp_cpools_leached_col[c, l] = 0.0
        end
        for j in 1:nlevdecomp
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                cf.decomp_cpools_leached_col[c, l] =
                    cf.decomp_cpools_leached_col[c, l] +
                    cf.decomp_cpools_transport_tendency_col[c, j, l] * dzsoi_decomp_vals[j]
            end
        end
        for c in bounds_col
            mask_bgc_soilc[c] || continue
            cf.som_c_leached_col[c] = cf.som_c_leached_col[c] + cf.decomp_cpools_leached_col[c, l]
        end
    end

    # Soil organic matter heterotrophic respiration (SOMHR)
    for k in 1:ndecomp_cascade_transitions
        if is_soil[cascade_donor_pool[k]]
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                cf.somhr_col[c] = cf.somhr_col[c] + cf.decomp_cascade_hr_col[c, k]
            end
        end
    end

    # Litter heterotrophic respiration (LITHR)
    for k in 1:ndecomp_cascade_transitions
        if is_litter[cascade_donor_pool[k]]
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                cf.lithr_col[c] = cf.lithr_col[c] + cf.decomp_cascade_hr_col[c, k]
            end
        end
    end

    # Coarse woody debris heterotrophic respiration (CWDHR)
    for k in 1:ndecomp_cascade_transitions
        if is_cwd[cascade_donor_pool[k]]
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                cf.cwdhr_col[c] = cf.cwdhr_col[c] + cf.decomp_cascade_hr_col[c, k]
            end
        end
    end

    # Microbial heterotrophic respiration (MICHR)
    for k in 1:ndecomp_cascade_transitions
        if is_microbe[cascade_donor_pool[k]]
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                cf.michr_col[c] = cf.michr_col[c] + cf.decomp_cascade_hr_col[c, k]
            end
        end
    end

    # Total heterotrophic respiration (HR)
    for c in bounds_col
        mask_bgc_soilc[c] || continue
        cf.hr_col[c] =
            cf.michr_col[c] +
            cf.cwdhr_col[c] +
            cf.lithr_col[c] +
            cf.somhr_col[c]
    end

    return nothing
end

"""
    soil_bgc_carbon_flux_init_history!(cf, bounds_col)

Stub for history field registration (no-op in Julia port).
Corresponds to `InitHistory` in the Fortran source.
"""
function soil_bgc_carbon_flux_init_history!(cf::SoilBiogeochemCarbonFluxData,
                                             bounds_col::UnitRange{Int})
    return nothing
end

"""
    soil_bgc_carbon_flux_restart!(cf, bounds_col)

Stub for restart read/write (no-op in Julia port).
Corresponds to `Restart` in the Fortran source.
"""
function soil_bgc_carbon_flux_restart!(cf::SoilBiogeochemCarbonFluxData,
                                        bounds_col::UnitRange{Int})
    return nothing
end
