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
                                  sf::Real,
                                  sf_no3::Real)
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
# --- (A) Total soil water (per-column REDUCTION over levels) ---
# Each thread accumulates ascending j into a LOCAL (preserving host add order),
# then writes once. cmin/cmax gate the active bounds range.
@kernel function _nleach_tot_water_kernel!(tot_water, @Const(mask), @Const(h2osoi_liq),
                                           nlevsoi::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        acc = tot_water[c]
        for j in 1:nlevsoi
            acc = acc + h2osoi_liq[c, j]
        end
        tot_water[c] = acc
    end
end

nleach_tot_water!(tot_water, mask, h2osoi_liq, nlevsoi::Int, cmin::Int, cmax::Int) =
    _launch!(_nleach_tot_water_kernel!, tot_water, mask, h2osoi_liq, nlevsoi, cmin, cmax)

# --- (B) Surface water to runoff-mixing depth (per-column REDUCTION over levels) ---
# Each thread accumulates ascending j into a LOCAL with the same depth-conditional
# branches as the host (zisoi/depth comparisons are uniform per-j, computed in-thread).
@kernel function _nleach_surface_water_kernel!(surface_water, @Const(mask),
                                               @Const(h2osoi_liq), @Const(col_dz),
                                               @Const(zisoi), depth_runoff_Nloss,
                                               nlevsoi::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        acc = surface_water[c]
        for j in 1:nlevsoi
            if zisoi[j+1] <= depth_runoff_Nloss
                acc = acc + h2osoi_liq[c, j]
            elseif zisoi[j] < depth_runoff_Nloss
                acc = acc + h2osoi_liq[c, j] * ((depth_runoff_Nloss - zisoi[j]) / col_dz[c, j])
            end
        end
        surface_water[c] = acc
    end
end

nleach_surface_water!(surface_water, mask, h2osoi_liq, col_dz, zisoi,
                      depth_runoff_Nloss, nlevsoi::Int, cmin::Int, cmax::Int) =
    _launch!(_nleach_surface_water_kernel!, surface_water, mask, h2osoi_liq, col_dz, zisoi,
             depth_runoff_Nloss, nlevsoi, cmin, cmax)

# --- (C) drain_tot gather (per-column, masked) ---
@kernel function _nleach_drain_tot_kernel!(drain_tot, @Const(mask), @Const(qflx_drain),
                                           cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        drain_tot[c] = qflx_drain[c]
    end
end

nleach_drain_tot!(drain_tot, mask, qflx_drain, cmin::Int, cmax::Int) =
    _launch!(_nleach_drain_tot_kernel!, drain_tot, mask, qflx_drain, cmin, cmax)

# --- (D) NITRIF_DENITRIF OFF leaching (per-(c,j), independent) ---
@kernel function _nleach_off_kernel!(sminn_leached_vr_col, @Const(mask),
                                     @Const(h2osoi_liq), @Const(sminn_vr_col),
                                     @Const(col_dz), @Const(tot_water), @Const(drain_tot),
                                     sf, dt, cmin::Int, cmax::Int)
    c, j = @index(Global, NTuple)
    @inbounds if cmin <= c <= cmax && mask[c]
        # Calculate the dissolved mineral N concentration (gN/kg water)
        # assumes that 10% of mineral nitrogen is soluble
        disn_conc = 0.0
        if h2osoi_liq[c, j] > 0.0
            disn_conc = (sf * sminn_vr_col[c, j] * col_dz[c, j]) / h2osoi_liq[c, j]
        end

        # Calculate the N leaching flux as a function of the dissolved
        # concentration and the sub-surface drainage flux
        sminn_leached_vr_col[c, j] = disn_conc * drain_tot[c] * h2osoi_liq[c, j] / (tot_water[c] * col_dz[c, j])

        # Limit the flux based on current sminn state
        # only let at most the assumed soluble fraction
        # of sminn be leached on any given timestep
        sminn_leached_vr_col[c, j] = min(sminn_leached_vr_col[c, j], (sf * sminn_vr_col[c, j]) / dt)

        # Limit the flux to a positive value
        sminn_leached_vr_col[c, j] = max(sminn_leached_vr_col[c, j], 0.0)
    end
end

function nleach_off!(sminn_leached_vr_col, mask, h2osoi_liq, sminn_vr_col, col_dz,
                     tot_water, drain_tot, sf, dt, nlevdecomp::Int, cmin::Int, cmax::Int)
    _launch!(_nleach_off_kernel!, sminn_leached_vr_col, mask, h2osoi_liq, sminn_vr_col,
             col_dz, tot_water, drain_tot, sf, dt, cmin, cmax;
             ndrange = (size(sminn_leached_vr_col, 1), nlevdecomp))
end

# --- (E) NITRIF_DENITRIF ON leaching + runoff (per-(c,j), independent) ---
@kernel function _nleach_on_kernel!(smin_no3_leached_vr_col, smin_no3_runoff_vr_col,
                                    @Const(mask), @Const(h2osoi_liq),
                                    @Const(smin_no3_vr_col), @Const(col_dz),
                                    @Const(tot_water), @Const(surface_water),
                                    @Const(drain_tot), @Const(qflx_surf), @Const(zisoi),
                                    sf_no3, dt, depth_runoff_Nloss, cmin::Int, cmax::Int)
    c, j = @index(Global, NTuple)
    @inbounds if cmin <= c <= cmax && mask[c]
        # Calculate the dissolved mineral N concentration (gN/kg water)
        disn_conc = 0.0
        if h2osoi_liq[c, j] > 0.0
            disn_conc = (sf_no3 * smin_no3_vr_col[c, j] * col_dz[c, j]) / h2osoi_liq[c, j]
        end

        # Calculate the N leaching flux as a function of the dissolved
        # concentration and the sub-surface drainage flux
        smin_no3_leached_vr_col[c, j] = disn_conc * drain_tot[c] * h2osoi_liq[c, j] / (tot_water[c] * col_dz[c, j])

        # Ensure that leaching rate isn't larger than soil N pool
        smin_no3_leached_vr_col[c, j] = min(smin_no3_leached_vr_col[c, j], smin_no3_vr_col[c, j] / dt)

        # Limit the leaching flux to a positive value
        smin_no3_leached_vr_col[c, j] = max(smin_no3_leached_vr_col[c, j], 0.0)

        # Calculate the N loss from surface runoff, assuming a shallow
        # mixing of surface waters into soil and removal based on runoff
        # Note: zisoi[j+1] = Fortran zisoi(j), zisoi[j] = Fortran zisoi(j-1)
        if zisoi[j+1] <= depth_runoff_Nloss
            smin_no3_runoff_vr_col[c, j] = disn_conc * qflx_surf[c] *
                h2osoi_liq[c, j] / (surface_water[c] * col_dz[c, j])
        elseif zisoi[j] < depth_runoff_Nloss
            smin_no3_runoff_vr_col[c, j] = disn_conc * qflx_surf[c] *
                h2osoi_liq[c, j] * ((depth_runoff_Nloss - zisoi[j]) /
                col_dz[c, j]) / (surface_water[c] * (depth_runoff_Nloss - zisoi[j]))
        else
            smin_no3_runoff_vr_col[c, j] = 0.0
        end

        # Ensure that runoff rate isn't larger than soil N pool
        smin_no3_runoff_vr_col[c, j] = min(smin_no3_runoff_vr_col[c, j],
            smin_no3_vr_col[c, j] / dt - smin_no3_leached_vr_col[c, j])

        # Limit the flux to a positive value
        smin_no3_runoff_vr_col[c, j] = max(smin_no3_runoff_vr_col[c, j], 0.0)

        # Limit the flux based on current smin_no3 state
        # only let at most the assumed soluble fraction
        # of smin_no3 be leached on any given timestep
        smin_no3_leached_vr_col[c, j] = min(smin_no3_leached_vr_col[c, j], (sf_no3 * smin_no3_vr_col[c, j]) / dt)

        # Limit the flux to a positive value
        smin_no3_leached_vr_col[c, j] = max(smin_no3_leached_vr_col[c, j], 0.0)
    end
end

function nleach_on!(smin_no3_leached_vr_col, smin_no3_runoff_vr_col, mask, h2osoi_liq,
                    smin_no3_vr_col, col_dz, tot_water, surface_water, drain_tot,
                    qflx_surf, zisoi, sf_no3, dt, depth_runoff_Nloss,
                    nlevdecomp::Int, cmin::Int, cmax::Int)
    _launch!(_nleach_on_kernel!, smin_no3_leached_vr_col, smin_no3_runoff_vr_col, mask,
             h2osoi_liq, smin_no3_vr_col, col_dz, tot_water, surface_water, drain_tot,
             qflx_surf, zisoi, sf_no3, dt, depth_runoff_Nloss, cmin, cmax;
             ndrange = (size(smin_no3_leached_vr_col, 1), nlevdecomp))
end

function n_leaching!(
    nf::SoilBiogeochemNitrogenFluxData,
    ns::SoilBiogeochemNitrogenStateData,
    params::NLeachingParams;
    mask_bgc_soilc::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevdecomp::Int,
    nlevsoi::Int,
    dt::Real,
    # Water state arrays
    h2osoi_liq::AbstractMatrix{<:Real},
    # Water flux arrays
    qflx_drain::AbstractVector{<:Real},
    qflx_surf::AbstractVector{<:Real},
    # Column geometry
    col_dz::AbstractMatrix{<:Real},
    # Vertical coordinate
    zisoi::AbstractVector{<:Real},
    # Control flag
    use_nitrif_denitrif::Bool)

    isempty(bounds) && return nothing

    # Constant: depth over which runoff mixes with soil water for N loss to runoff
    depth_runoff_Nloss = 0.05  # (m)

    cmin = first(bounds)
    cmax = last(bounds)

    # --- Allocate local work arrays ---
    FT = eltype(h2osoi_liq)
    tot_water     = zeros(FT, last(bounds))
    surface_water = zeros(FT, last(bounds))
    drain_tot     = zeros(FT, last(bounds))

    # --- Select soluble fraction based on mode (config flag host-resolved) ---
    if !use_nitrif_denitrif
        sf = params.sf
    else
        sf_no3 = params.sf_no3
    end

    # --- (A) Calculate the total soil water (per-column reduction over levels) ---
    nleach_tot_water!(tot_water, mask_bgc_soilc, h2osoi_liq, nlevsoi, cmin, cmax)

    # --- (B) For runoff calculation; total water to a given depth (reduction) ---
    nleach_surface_water!(surface_water, mask_bgc_soilc, h2osoi_liq, col_dz, zisoi,
                          depth_runoff_Nloss, nlevsoi, cmin, cmax)

    # --- (C) Set drain_tot ---
    nleach_drain_tot!(drain_tot, mask_bgc_soilc, qflx_drain, cmin, cmax)

    if !use_nitrif_denitrif

        # ----------------------------------------
        # --------- NITRIF_DENITRIF OFF ----------
        # ----------------------------------------
        nleach_off!(nf.sminn_leached_vr_col, mask_bgc_soilc, h2osoi_liq,
                    ns.sminn_vr_col, col_dz, tot_water, drain_tot, sf, dt,
                    nlevdecomp, cmin, cmax)

    else

        # ----------------------------------------
        # --------- NITRIF_DENITRIF ON -----------
        # ----------------------------------------
        nleach_on!(nf.smin_no3_leached_vr_col, nf.smin_no3_runoff_vr_col, mask_bgc_soilc,
                   h2osoi_liq, ns.smin_no3_vr_col, col_dz, tot_water, surface_water,
                   drain_tot, qflx_surf, zisoi, sf_no3, dt, depth_runoff_Nloss,
                   nlevdecomp, cmin, cmax)
    end

    return nothing
end
