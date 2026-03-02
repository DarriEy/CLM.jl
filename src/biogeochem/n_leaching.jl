# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemNLeachingMod.F90
# Mineral nitrogen leaching: soluble mineral N loss via subsurface drainage
# and (with nitrif_denitrif) NO3 loss via surface runoff.
#
# Public functions:
#   n_leaching_read_params!  — Read/set N-leaching parameters
#   n_leaching!              — Calculate N leaching and runoff losses
# ==========================================================================

# ---------------------------------------------------------------------------
# NLeachingParams — module parameters
# Ported from params_type in SoilBiogeochemNLeachingMod.F90
# ---------------------------------------------------------------------------

"""
    NLeachingParams

N-leaching parameters. Holds soluble fractions of mineral N and NO3.

Ported from `params_type` in `SoilBiogeochemNLeachingMod.F90`.
"""
Base.@kwdef mutable struct NLeachingParams
    sf     ::Float64 = 0.1    # soluble fraction of mineral N (unitless)
    sf_no3 ::Float64 = 1.0    # soluble fraction of NO3 (unitless)
end

# ---------------------------------------------------------------------------
# n_leaching_read_params! — read/set parameters
# Ported from readParams in SoilBiogeochemNLeachingMod.F90
# ---------------------------------------------------------------------------

"""
    n_leaching_read_params!(params; sf, sf_no3)

Set N-leaching parameters from keyword arguments.
Corresponds to `readParams` in the Fortran source (reads from params file).
"""
function n_leaching_read_params!(params::NLeachingParams;
                                  sf::Float64,
                                  sf_no3::Float64)
    params.sf     = sf
    params.sf_no3 = sf_no3
    return nothing
end

# ---------------------------------------------------------------------------
# n_leaching! — main calculation
# Ported from SoilBiogeochemNLeaching in SoilBiogeochemNLeachingMod.F90
# ---------------------------------------------------------------------------

"""
    n_leaching!(nf, ns, params; kwargs...)

Calculate nitrogen leaching rates as a function of soluble mineral N and
total soil water outflow. Updates `sminn_leached_vr` (when
`use_nitrif_denitrif=false`) or `smin_no3_leached_vr` and
`smin_no3_runoff_vr` (when `use_nitrif_denitrif=true`).

Corresponds to `SoilBiogeochemNLeaching` in the Fortran source.
"""
function n_leaching!(
    nf::SoilBiogeochemNitrogenFluxData,
    ns::SoilBiogeochemNitrogenStateData,
    params::NLeachingParams;
    mask_bgc_soilc::BitVector,
    bounds::UnitRange{Int},
    nlevdecomp::Int,
    nlevsoi::Int,
    dt::Float64,
    # Water state arrays
    h2osoi_liq::Matrix{Float64},
    # Water flux arrays
    qflx_drain::Vector{Float64},
    qflx_surf::Vector{Float64},
    # Column geometry
    col_dz::Matrix{Float64},
    # Vertical coordinate
    zisoi::Vector{Float64},
    # Control flag
    use_nitrif_denitrif::Bool)

    # Constant: depth over which runoff mixes with soil water for N loss to runoff
    depth_runoff_Nloss = 0.05  # (m)

    nc = length(bounds)

    # --- Allocate local work arrays ---
    tot_water     = zeros(Float64, last(bounds))
    surface_water = zeros(Float64, last(bounds))
    drain_tot     = zeros(Float64, last(bounds))

    # --- Select soluble fraction based on mode ---
    if !use_nitrif_denitrif
        sf = params.sf
    else
        sf_no3 = params.sf_no3
    end

    # --- Calculate the total soil water ---
    for j in 1:nlevsoi
        for c in bounds
            mask_bgc_soilc[c] || continue
            tot_water[c] = tot_water[c] + h2osoi_liq[c, j]
        end
    end

    # --- For runoff calculation; calculate total water to a given depth ---
    # Note: zisoi is 1-indexed with zisoi[1]=0 (surface), zisoi[j+1]=bottom of layer j
    # Fortran zisoi(j) → Julia zisoi[j+1], Fortran zisoi(j-1) → Julia zisoi[j]
    for j in 1:nlevsoi
        if zisoi[j+1] <= depth_runoff_Nloss
            for c in bounds
                mask_bgc_soilc[c] || continue
                surface_water[c] = surface_water[c] + h2osoi_liq[c, j]
            end
        elseif zisoi[j] < depth_runoff_Nloss
            for c in bounds
                mask_bgc_soilc[c] || continue
                surface_water[c] = surface_water[c] + h2osoi_liq[c, j] * ((depth_runoff_Nloss - zisoi[j]) / col_dz[c, j])
            end
        end
    end

    # --- Set drain_tot ---
    for c in bounds
        mask_bgc_soilc[c] || continue
        drain_tot[c] = qflx_drain[c]
    end

    if !use_nitrif_denitrif

        # ----------------------------------------
        # --------- NITRIF_DENITRIF OFF ----------
        # ----------------------------------------

        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue

                # Calculate the dissolved mineral N concentration (gN/kg water)
                # assumes that 10% of mineral nitrogen is soluble
                disn_conc = 0.0
                if h2osoi_liq[c, j] > 0.0
                    disn_conc = (sf * ns.sminn_vr_col[c, j] * col_dz[c, j]) / h2osoi_liq[c, j]
                end

                # Calculate the N leaching flux as a function of the dissolved
                # concentration and the sub-surface drainage flux
                nf.sminn_leached_vr_col[c, j] = disn_conc * drain_tot[c] * h2osoi_liq[c, j] / (tot_water[c] * col_dz[c, j])

                # Limit the flux based on current sminn state
                # only let at most the assumed soluble fraction
                # of sminn be leached on any given timestep
                nf.sminn_leached_vr_col[c, j] = min(nf.sminn_leached_vr_col[c, j], (sf * ns.sminn_vr_col[c, j]) / dt)

                # Limit the flux to a positive value
                nf.sminn_leached_vr_col[c, j] = max(nf.sminn_leached_vr_col[c, j], 0.0)
            end
        end

    else

        # ----------------------------------------
        # --------- NITRIF_DENITRIF ON -----------
        # ----------------------------------------

        for j in 1:nlevdecomp
            for c in bounds
                mask_bgc_soilc[c] || continue

                # Calculate the dissolved mineral N concentration (gN/kg water)
                disn_conc = 0.0
                if h2osoi_liq[c, j] > 0.0
                    disn_conc = (sf_no3 * ns.smin_no3_vr_col[c, j] * col_dz[c, j]) / h2osoi_liq[c, j]
                end

                # Calculate the N leaching flux as a function of the dissolved
                # concentration and the sub-surface drainage flux
                nf.smin_no3_leached_vr_col[c, j] = disn_conc * drain_tot[c] * h2osoi_liq[c, j] / (tot_water[c] * col_dz[c, j])

                # Ensure that leaching rate isn't larger than soil N pool
                nf.smin_no3_leached_vr_col[c, j] = min(nf.smin_no3_leached_vr_col[c, j], ns.smin_no3_vr_col[c, j] / dt)

                # Limit the leaching flux to a positive value
                nf.smin_no3_leached_vr_col[c, j] = max(nf.smin_no3_leached_vr_col[c, j], 0.0)

                # Calculate the N loss from surface runoff, assuming a shallow
                # mixing of surface waters into soil and removal based on runoff
                # Note: zisoi[j+1] = Fortran zisoi(j), zisoi[j] = Fortran zisoi(j-1)
                if zisoi[j+1] <= depth_runoff_Nloss
                    nf.smin_no3_runoff_vr_col[c, j] = disn_conc * qflx_surf[c] *
                        h2osoi_liq[c, j] / (surface_water[c] * col_dz[c, j])
                elseif zisoi[j] < depth_runoff_Nloss
                    nf.smin_no3_runoff_vr_col[c, j] = disn_conc * qflx_surf[c] *
                        h2osoi_liq[c, j] * ((depth_runoff_Nloss - zisoi[j]) /
                        col_dz[c, j]) / (surface_water[c] * (depth_runoff_Nloss - zisoi[j]))
                else
                    nf.smin_no3_runoff_vr_col[c, j] = 0.0
                end

                # Ensure that runoff rate isn't larger than soil N pool
                nf.smin_no3_runoff_vr_col[c, j] = min(nf.smin_no3_runoff_vr_col[c, j],
                    ns.smin_no3_vr_col[c, j] / dt - nf.smin_no3_leached_vr_col[c, j])

                # Limit the flux to a positive value
                nf.smin_no3_runoff_vr_col[c, j] = max(nf.smin_no3_runoff_vr_col[c, j], 0.0)

                # Limit the flux based on current smin_no3 state
                # only let at most the assumed soluble fraction
                # of smin_no3 be leached on any given timestep
                nf.smin_no3_leached_vr_col[c, j] = min(nf.smin_no3_leached_vr_col[c, j], (sf_no3 * ns.smin_no3_vr_col[c, j]) / dt)

                # Limit the flux to a positive value
                nf.smin_no3_leached_vr_col[c, j] = max(nf.smin_no3_leached_vr_col[c, j], 0.0)
            end
        end
    end

    return nothing
end
