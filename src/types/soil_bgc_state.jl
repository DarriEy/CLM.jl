# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemStateType.F90
# Soil biogeochemistry state data type allocation and initialization
# ==========================================================================

"""
    SoilBiogeochemStateData

Soil biogeochemistry state data structure. Holds vertical profiles for litter
inputs, fraction of potential immobilization/GPP, SOM transport coefficients,
N fixation/deposition profiles, and plant N demand.

Ported from `soilbiogeochem_state_type` in `SoilBiogeochemStateType.F90`.
"""
Base.@kwdef mutable struct SoilBiogeochemStateData
    # --- Patch-level vertical profiles (patch × nlevdecomp_full) ---
    leaf_prof_patch             ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # (1/m) profile of leaves
    froot_prof_patch            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # (1/m) profile of fine roots
    croot_prof_patch            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # (1/m) profile of coarse roots
    stem_prof_patch             ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # (1/m) profile of stems

    # --- Column-level vertically-resolved (col × nlevdecomp_full) ---
    fpi_vr_col                  ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # (no units) fraction of potential immobilization vr
    nfixation_prof_col          ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # (1/m) profile for N fixation additions
    ndep_prof_col               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # (1/m) profile for N deposition additions
    som_adv_coef_col            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # (m/s) SOM advective flux
    som_diffus_coef_col         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # (m2/s) SOM diffusivity due to bio/cryo-turbation

    # --- Column-level 1D ---
    fpi_col                     ::Vector{Float64} = Float64[]  # (no units) fraction of potential immobilization
    fpg_col                     ::Vector{Float64} = Float64[]  # (no units) fraction of potential gpp
    plant_ndemand_col           ::Vector{Float64} = Float64[]  # column-level plant N demand

    # --- Transition-level 1D ---
    nue_decomp_cascade_col      ::Vector{Float64} = Float64[]  # (gN/gN) N use efficiency for a given transition
end

# ---------------------------------------------------------------------------
# Helper constructors (reuse if already defined from carbon/nitrogen state)
# ---------------------------------------------------------------------------
if !@isdefined(nanvec)
    nanvec(n) = fill(NaN, n)
end
if !@isdefined(nanmat)
    nanmat(r, c) = fill(NaN, r, c)
end

"""
    soil_bgc_state_init!(st, nc, np, nlevdecomp_full, ndecomp_cascade_transitions)

Allocate all fields of `SoilBiogeochemStateData`.
Corresponds to `InitAllocate` in the Fortran source.
"""
function soil_bgc_state_init!(st::SoilBiogeochemStateData,
                               nc::Int, np::Int,
                               nlevdecomp_full::Int,
                               ndecomp_cascade_transitions::Int)

    # --- Patch-level 2D (patch × nlevdecomp_full) --- initialized to SPVAL
    st.leaf_prof_patch     = fill(SPVAL, np, nlevdecomp_full)
    st.froot_prof_patch    = fill(SPVAL, np, nlevdecomp_full)
    st.croot_prof_patch    = fill(SPVAL, np, nlevdecomp_full)
    st.stem_prof_patch     = fill(SPVAL, np, nlevdecomp_full)

    # --- Column-level 2D (col × nlevdecomp_full) ---
    st.fpi_vr_col          = nanmat(nc, nlevdecomp_full)
    st.nfixation_prof_col  = fill(SPVAL, nc, nlevdecomp_full)
    st.ndep_prof_col       = fill(SPVAL, nc, nlevdecomp_full)
    st.som_adv_coef_col    = fill(SPVAL, nc, nlevdecomp_full)
    st.som_diffus_coef_col = fill(SPVAL, nc, nlevdecomp_full)

    # --- Column-level 1D ---
    st.fpi_col             = nanvec(nc)
    st.fpg_col             = nanvec(nc)
    st.plant_ndemand_col   = nanvec(nc)

    # --- Transition-level 1D ---
    st.nue_decomp_cascade_col = nanvec(ndecomp_cascade_transitions)

    return nothing
end

"""
    soil_bgc_state_init_cold!(st, bounds_col, nlevdecomp_full;
                               mask_special=nothing, mask_soil_crop=nothing)

Cold-start initialization of soil biogeochemistry state fields.
For special landunits: sets fpi_col, fpg_col, fpi_vr_col to SPVAL.
For soil/crop columns: zeros out fpi_vr_col, som_adv_coef_col,
som_diffus_coef_col, nfixation_prof_col, ndep_prof_col.

Corresponds to `InitCold` in the Fortran source.
"""
function soil_bgc_state_init_cold!(st::SoilBiogeochemStateData,
                                    bounds_col::UnitRange{Int},
                                    nlevdecomp_full::Int;
                                    mask_special::Union{BitVector,Nothing}=nothing,
                                    mask_soil_crop::Union{BitVector,Nothing}=nothing)

    for c in bounds_col
        # Special landunits: set to SPVAL
        if mask_special !== nothing && mask_special[c]
            st.fpi_col[c] = SPVAL
            st.fpg_col[c] = SPVAL
            for j in 1:nlevdecomp_full
                st.fpi_vr_col[c, j] = SPVAL
            end
        end

        # Soil/crop landunits: zero out profiles
        is_soil_crop = mask_soil_crop === nothing || mask_soil_crop[c]
        if is_soil_crop
            for j in 1:nlevdecomp_full
                st.fpi_vr_col[c, j]          = 0.0
                st.som_adv_coef_col[c, j]    = 0.0
                st.som_diffus_coef_col[c, j] = 0.0
                st.nfixation_prof_col[c, j]  = 0.0
                st.ndep_prof_col[c, j]       = 0.0
            end
        end
    end

    return nothing
end

"""
    soil_bgc_state_init_history!(st, bounds_col)

Stub for history field registration (no-op in Julia port).
Corresponds to `InitHistory` in the Fortran source.
"""
function soil_bgc_state_init_history!(st::SoilBiogeochemStateData,
                                      bounds_col::UnitRange{Int})
    return nothing
end

"""
    soil_bgc_state_restart!(st, bounds_col)

Stub for restart read/write (no-op in Julia port).
Corresponds to `Restart` in the Fortran source.
"""
function soil_bgc_state_restart!(st::SoilBiogeochemStateData,
                                  bounds_col::UnitRange{Int})
    return nothing
end

"""
    get_spinup_latitude_term(latitude)

Calculate a logistic function to scale spinup factors so that spinup is more
accelerated in high latitude regions.

Ported from `get_spinup_latitude_term` in `SoilBiogeochemStateType.F90`.
"""
function get_spinup_latitude_term(latitude::Float64)
    return 1.0 + 50.0 / (1.0 + exp(-0.15 * (abs(latitude) - 60.0)))
end
