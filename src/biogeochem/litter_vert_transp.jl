# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemLittVertTranspMod.F90 (524 lines)
# Calculate vertical mixing of all decomposing C and N pools.
# Advection-diffusion code based on algorithm in Patankar (1980).
# Initial Fortran code by C. Koven and W. Riley.
#
# Public functions:
#   litter_vert_transp_readparams! -- Read parameters from params struct
#   litter_vert_transp!            -- Main vertical transport routine
#
# Public module-level variables:
#   som_adv_flux           -- SOM advective flux (m/s)
#   max_depth_cryoturb     -- maximum depth of cryoturbation (m)
# ==========================================================================

# ---------------------------------------------------------------------------
# LitterVertTranspParams -- parameters for vertical transport
# Ported from params_type in SoilBiogeochemLittVertTranspMod.F90
# ---------------------------------------------------------------------------

Base.@kwdef mutable struct LitterVertTranspParams
    som_diffus::Float64 = 0.0                   # Soil organic matter diffusion (m^2/s)
    cryoturb_diffusion_k::Float64 = 0.0         # Cryoturbation diffusive constant (m^2/s)
    max_altdepth_cryoturbation::Float64 = 0.0    # Maximum active layer thickness for cryoturbation (m)
end

# Module-level public parameters (matches Fortran module variables)
const LITTER_VERT_TRANSP_PARAMS = LitterVertTranspParams()

# Module-level public constants (set at module load; can be overridden by user)
const som_adv_flux_ref = Ref(0.0)           # SOM advective flux (m/s)
const max_depth_cryoturb_ref = Ref(3.0)     # Maximum depth of cryoturbation (m)

# ---------------------------------------------------------------------------
# litter_vert_transp_readparams! -- read/set parameters
# Ported from readParams in SoilBiogeochemLittVertTranspMod.F90
# ---------------------------------------------------------------------------

"""
    litter_vert_transp_readparams!(params; som_diffus, cryoturb_diffusion_k,
                                   max_altdepth_cryoturbation)

Set parameters for litter vertical transport.
Ported from `readParams` in `SoilBiogeochemLittVertTranspMod.F90`.
"""
function litter_vert_transp_readparams!(params::LitterVertTranspParams;
                                        som_diffus::Real,
                                        cryoturb_diffusion_k::Real,
                                        max_altdepth_cryoturbation::Real)
    params.som_diffus = som_diffus
    params.cryoturb_diffusion_k = cryoturb_diffusion_k
    params.max_altdepth_cryoturbation = max_altdepth_cryoturbation
    return nothing
end

# ---------------------------------------------------------------------------
# patankar_A -- the "A" function from Patankar Table 5.2 pg 95
# ---------------------------------------------------------------------------

"""
    patankar_A(pe)

Patankar's "A" function (Table 5.2, pg 95).
Returns max(0, (1 - 0.1*|pe|)^5).
"""
@inline function patankar_A(pe::Real)
    return max(0.0, (1.0 - 0.1 * abs(pe))^5)
end

# ---------------------------------------------------------------------------
# litter_vert_transp! -- main vertical transport routine
# Ported from SoilBiogeochemLittVertTransp in SoilBiogeochemLittVertTranspMod.F90
# ---------------------------------------------------------------------------

"""
    litter_vert_transp!(cs, cf, ns, nf, st, col, grc, cascade_con, params;
                         mask_bgc_soilc, bounds, dtime, nlevdecomp, ndecomp_pools,
                         zsoi_vals, dzsoi_decomp_vals, zisoi_vals,
                         spinup_state, use_soil_matrixcn, som_adv_flux_val,
                         max_depth_cryoturb_val)

Calculate vertical mixing (advection + diffusion) of soil and litter pools
for both C and N. Also reconcile sources and sinks calculated in CStateUpdate1
and NStateUpdate1.

Ported from `SoilBiogeochemLittVertTransp` in `SoilBiogeochemLittVertTranspMod.F90`.

# Arguments
- `cs::SoilBiogeochemCarbonStateData`    : carbon state (in/out)
- `cf::SoilBiogeochemCarbonFluxData`     : carbon flux (in/out)
- `ns::SoilBiogeochemNitrogenStateData`  : nitrogen state (in/out)
- `nf::SoilBiogeochemNitrogenFluxData`   : nitrogen flux (in/out)
- `st::SoilBiogeochemStateData`          : soil biogeochem state (in/out)
- `col::ColumnData`                      : column data (in)
- `grc::GridcellData`                    : gridcell data (in)
- `cascade_con::DecompCascadeConData`    : decomposition cascade configuration (in)
- `active_layer::ActiveLayerData`        : active layer data (in)
- `params::LitterVertTranspParams`       : parameters (in)
- `mask_bgc_soilc::BitVector`            : soil column mask (in)
- `bounds::UnitRange{Int}`               : column index range (in)
- `dtime::Float64`                       : timestep size (s) (in)
- `nlevdecomp::Int`                      : number of decomposition levels (in)
- `ndecomp_pools::Int`                   : number of decomposition pools (in)
- `zsoi_vals::Vector{<:Real}`           : soil node depths, 1-indexed (m) (in)
- `dzsoi_decomp_vals::Vector{<:Real}`   : decomposition layer thicknesses (m) (in)
- `zisoi_vals::Vector{<:Real}`          : soil interface depths (m), 1-indexed with zisoi_vals[1]=0 (surface), zisoi_vals[j+1]=bottom of layer j (in)
- `spinup_state::Int`                    : spinup state flag (in)
- `use_soil_matrixcn::Bool`              : whether to use soil matrix solution (in)
- `som_adv_flux_val::Float64`            : SOM advective flux value (m/s) (in)
- `max_depth_cryoturb_val::Float64`      : max depth of cryoturbation (m) (in)
"""
function litter_vert_transp!(
        cs::SoilBiogeochemCarbonStateData,
        cf::SoilBiogeochemCarbonFluxData,
        ns::SoilBiogeochemNitrogenStateData,
        nf::SoilBiogeochemNitrogenFluxData,
        st::SoilBiogeochemStateData,
        col::ColumnData,
        grc::GridcellData,
        cascade_con::DecompCascadeConData,
        active_layer::ActiveLayerData,
        params::LitterVertTranspParams;
        mask_bgc_soilc::BitVector,
        bounds::UnitRange{Int},
        dtime::Real,
        nlevdecomp::Int,
        ndecomp_pools::Int,
        zsoi_vals::Vector{<:Real},
        dzsoi_decomp_vals::Vector{<:Real},
        zisoi_vals::Vector{<:Real},
        spinup_state::Int = 0,
        use_soil_matrixcn::Bool = false,
        som_adv_flux_val::Real = 0.0,
        max_depth_cryoturb_val::Real = 3.0)

    # Unpack cascade configuration
    is_cwd        = cascade_con.is_cwd
    spinup_factor = cascade_con.spinup_factor

    # Unpack active layer data
    altmax          = active_layer.altmax_col
    altmax_lastyear = active_layer.altmax_lastyear_col

    # Unpack SoilBiogeochemState output arrays
    som_adv_coef    = st.som_adv_coef_col
    som_diffus_coef = st.som_diffus_coef_col

    # Local parameters from params struct
    som_diffus_val                 = params.som_diffus
    cryoturb_diffusion_k           = params.cryoturb_diffusion_k
    max_altdepth_cryoturbation     = params.max_altdepth_cryoturbation

    nc = length(bounds)
    epsilon = 1.0e-30

    # Number of tracer types: always C and N; isotopes not implemented in this port
    ntype = 2

    # ------ Compute diffusivity / advection coefficients ------
    for c in bounds
        mask_bgc_soilc[c] || continue

        alt_max_val = max(altmax[c], altmax_lastyear[c])

        if alt_max_val <= max_altdepth_cryoturbation && alt_max_val > 0.0
            # Cryoturbation-dominated mixing
            for j in 1:(nlevdecomp + 1)
                if j <= col.nbedrock[c] + 1
                    if zisoi_vals[j + 1] < alt_max_val
                        som_diffus_coef[c, j] = cryoturb_diffusion_k
                        som_adv_coef[c, j] = 0.0
                    else
                        som_diffus_coef[c, j] = max(
                            cryoturb_diffusion_k *
                            (1.0 - (zisoi_vals[j + 1] - alt_max_val) /
                            (min(max_depth_cryoturb_val, zisoi_vals[col.nbedrock[c] + 2]) - alt_max_val)),
                            0.0)
                        som_adv_coef[c, j] = 0.0
                    end
                else
                    som_adv_coef[c, j] = 0.0
                    som_diffus_coef[c, j] = 0.0
                end
            end
        elseif alt_max_val > 0.0
            # Bioturbation: constant advection and diffusion
            for j in 1:(nlevdecomp + 1)
                if j <= col.nbedrock[c] + 1
                    som_adv_coef[c, j] = som_adv_flux_val
                    som_diffus_coef[c, j] = som_diffus_val
                else
                    som_adv_coef[c, j] = 0.0
                    som_diffus_coef[c, j] = 0.0
                end
            end
        else
            # Completely frozen soils: no mixing
            for j in 1:(nlevdecomp + 1)
                som_adv_coef[c, j] = 0.0
                som_diffus_coef[c, j] = 0.0
            end
        end
    end

    # Infer floating-point type for AD compatibility
    FT = eltype(zsoi_vals)

    # Extend zsoi with a virtual bottom node for boundary calculations
    # (Fortran typically has zsoi dimensioned to nlevgrnd > nlevdecomp)
    zsoi_ext = zeros(FT, nlevdecomp + 1)
    zsoi_ext[1:nlevdecomp] .= zsoi_vals[1:nlevdecomp]
    zsoi_ext[nlevdecomp + 1] = zisoi_vals[nlevdecomp + 1]

    # Set the distance between nodes
    dz_node = zeros(FT, nlevdecomp + 1)
    dz_node[1] = zsoi_ext[1]
    for j in 2:(nlevdecomp + 1)
        dz_node[j] = zsoi_ext[j] - zsoi_ext[j - 1]
    end

    # ------ Allocate local work arrays ------
    nj = nlevdecomp + 2  # levels 0:(nlevdecomp+1), stored as 1:(nlevdecomp+2)
    ncols = length(mask_bgc_soilc)

    diffus     = zeros(FT, ncols, nlevdecomp + 1)
    adv_flux   = zeros(FT, ncols, nlevdecomp + 1)
    a_tri      = zeros(FT, ncols, nj)  # indices 0..nlevdecomp+1 mapped to 1..nj
    b_tri      = zeros(FT, ncols, nj)
    c_tri      = zeros(FT, ncols, nj)
    r_tri      = zeros(FT, ncols, nj)
    d_p1_zp1   = zeros(FT, ncols, nlevdecomp + 1)
    d_m1_zm1   = zeros(FT, ncols, nlevdecomp + 1)
    f_p1       = zeros(FT, ncols, nlevdecomp + 1)
    f_m1       = zeros(FT, ncols, nlevdecomp + 1)
    pe_p1      = zeros(FT, ncols, nlevdecomp + 1)
    pe_m1      = zeros(FT, ncols, nlevdecomp + 1)
    conc_trcr  = zeros(FT, ncols, nj)  # indices 0..nlevdecomp+1 mapped to 1..nj

    # Offset helper: Fortran j=0..nlevdecomp+1 -> Julia index j+1=1..nlevdecomp+2
    # So conc_trcr(c,0) -> conc_trcr[c,1], conc_trcr(c,j) -> conc_trcr[c,j+1]
    # a_tri(c,0) -> a_tri[c,1], etc.

    # ------ Loop over tracer types (C, N) ------
    for i_type in 1:ntype

        # Select the appropriate state/flux arrays
        if i_type == 1  # Carbon
            conc_ptr          = cs.decomp_cpools_vr_col
            source            = cf.decomp_cpools_sourcesink_col
            trcr_tendency_ptr = cf.decomp_cpools_transport_tendency_col
        else  # i_type == 2, Nitrogen
            conc_ptr          = ns.decomp_npools_vr_col
            source            = nf.decomp_npools_sourcesink_col
            trcr_tendency_ptr = nf.decomp_npools_transport_tendency_col
        end

        for s in 1:ndecomp_pools
            if !is_cwd[s]
                # --- Compute diffusivity and Peclet numbers (only for first non-CWD pool
                #     when using soil matrix, otherwise for every pool) ---
                if !use_soil_matrixcn || s == 1

                    # Build diffus and adv_flux with spinup correction
                    for j in 1:(nlevdecomp + 1)
                        for c in bounds
                            mask_bgc_soilc[c] || continue

                            spinup_term = 1.0
                            if spinup_state >= 1
                                spinup_term = spinup_factor[s]
                            end
                            if abs(spinup_term - 1.0) > 1.0e-6
                                spinup_term = spinup_term * get_spinup_latitude_term(grc.latdeg[col.gridcell[c]])
                            end

                            if abs(som_adv_coef[c, j]) * spinup_term < epsilon
                                adv_flux[c, j] = epsilon
                            else
                                adv_flux[c, j] = som_adv_coef[c, j] * spinup_term
                            end

                            if abs(som_diffus_coef[c, j]) * spinup_term < epsilon
                                diffus[c, j] = epsilon
                            else
                                diffus[c, j] = som_diffus_coef[c, j] * spinup_term
                            end
                        end
                    end

                    # Set boundary values for conc_trcr
                    for c in bounds
                        mask_bgc_soilc[c] || continue
                        conc_trcr[c, 1] = 0.0  # j=0 -> index 1
                        for jj in (col.nbedrock[c] + 1):(nlevdecomp + 1)
                            conc_trcr[c, jj + 1] = 0.0  # j -> index j+1
                        end
                    end

                    # Set conc_trcr and compute D/dz, F, Pe throughout column
                    for j in 1:(nlevdecomp + 1)
                        for c in bounds
                            mask_bgc_soilc[c] || continue

                            if j <= nlevdecomp
                            conc_trcr[c, j + 1] = conc_ptr[c, j, s]  # j -> index j+1
                        end
                        # j > nlevdecomp: conc_trcr already set to 0 by boundary conditions above

                            if j == 1
                                d_m1_zm1[c, j] = 0.0
                                w_p1 = (zsoi_ext[j + 1] - zisoi_vals[j + 1]) / dz_node[j + 1]
                                if diffus[c, j + 1] > 0.0 && diffus[c, j] > 0.0
                                    d_p1 = 1.0 / ((1.0 - w_p1) / diffus[c, j] + w_p1 / diffus[c, j + 1])
                                else
                                    d_p1 = 0.0
                                end
                                d_p1_zp1[c, j] = d_p1 / dz_node[j + 1]
                                f_m1[c, j] = adv_flux[c, j]
                                f_p1[c, j] = adv_flux[c, j + 1]
                                pe_m1[c, j] = 0.0
                                pe_p1[c, j] = f_p1[c, j] / d_p1_zp1[c, j]

                            elseif j >= col.nbedrock[c] + 1
                                # Bottom boundary
                                w_m1 = (zisoi_vals[j] - zsoi_ext[j - 1]) / dz_node[j]
                                if diffus[c, j] > 0.0 && diffus[c, j - 1] > 0.0
                                    d_m1 = 1.0 / ((1.0 - w_m1) / diffus[c, j] + w_m1 / diffus[c, j - 1])
                                else
                                    d_m1 = 0.0
                                end
                                d_m1_zm1[c, j] = d_m1 / dz_node[j]
                                d_p1_zp1[c, j] = d_m1_zm1[c, j]
                                f_m1[c, j] = adv_flux[c, j]
                                f_p1[c, j] = 0.0
                                pe_m1[c, j] = f_m1[c, j] / d_m1_zm1[c, j]
                                pe_p1[c, j] = f_p1[c, j] / d_p1_zp1[c, j]

                            else
                                # Interior
                                w_m1 = (zisoi_vals[j] - zsoi_ext[j - 1]) / dz_node[j]
                                if diffus[c, j - 1] > 0.0 && diffus[c, j] > 0.0
                                    d_m1 = 1.0 / ((1.0 - w_m1) / diffus[c, j] + w_m1 / diffus[c, j - 1])
                                else
                                    d_m1 = 0.0
                                end
                                w_p1 = (zsoi_ext[j + 1] - zisoi_vals[j + 1]) / dz_node[j + 1]
                                if diffus[c, j + 1] > 0.0 && diffus[c, j] > 0.0
                                    d_p1 = 1.0 / ((1.0 - w_p1) / diffus[c, j] + w_p1 / diffus[c, j + 1])
                                else
                                    d_p1 = (1.0 - w_m1) * diffus[c, j] + w_p1 * diffus[c, j + 1]
                                end
                                d_m1_zm1[c, j] = d_m1 / dz_node[j]
                                d_p1_zp1[c, j] = d_p1 / dz_node[j + 1]
                                f_m1[c, j] = adv_flux[c, j]
                                f_p1[c, j] = adv_flux[c, j + 1]
                                pe_m1[c, j] = f_m1[c, j] / d_m1_zm1[c, j]
                                pe_p1[c, j] = f_p1[c, j] / d_p1_zp1[c, j]
                            end
                        end
                    end
                end  # !use_soil_matrixcn || s == 1

                # --- Calculate tridiagonal coefficients ---
                # j=0..nlevdecomp+1 mapped to Julia index j+1=1..nlevdecomp+2
                for j in 0:(nlevdecomp + 1)
                    jj = j + 1  # Julia index
                    for c in bounds
                        mask_bgc_soilc[c] || continue

                        local a_p_0::Float64
                        if j > 0 && j < nlevdecomp + 1
                            a_p_0 = dzsoi_decomp_vals[j] / dtime
                        else
                            a_p_0 = 0.0  # not used but needs a value
                        end

                        if j == 0
                            # Top layer (atmosphere boundary)
                            a_tri[c, jj] = 0.0
                            b_tri[c, jj] = 1.0
                            c_tri[c, jj] = -1.0
                            r_tri[c, jj] = 0.0
                        elseif j == 1
                            a_tri[c, jj] = -(d_m1_zm1[c, j] * patankar_A(pe_m1[c, j]) + max(f_m1[c, j], 0.0))
                            c_tri[c, jj] = -(d_p1_zp1[c, j] * patankar_A(pe_p1[c, j]) + max(-f_p1[c, j], 0.0))
                            b_tri[c, jj] = -a_tri[c, jj] - c_tri[c, jj] + a_p_0
                            r_tri[c, jj] = source[c, j, s] * dzsoi_decomp_vals[j] / dtime +
                                           (a_p_0 - adv_flux[c, j]) * conc_trcr[c, jj]
                        elseif j < nlevdecomp + 1
                            a_tri[c, jj] = -(d_m1_zm1[c, j] * patankar_A(pe_m1[c, j]) + max(f_m1[c, j], 0.0))
                            c_tri[c, jj] = -(d_p1_zp1[c, j] * patankar_A(pe_p1[c, j]) + max(-f_p1[c, j], 0.0))
                            b_tri[c, jj] = -a_tri[c, jj] - c_tri[c, jj] + a_p_0
                            r_tri[c, jj] = source[c, j, s] * dzsoi_decomp_vals[j] / dtime +
                                           a_p_0 * conc_trcr[c, jj]
                        else
                            # j == nlevdecomp+1; zero concentration gradient at bottom
                            a_tri[c, jj] = -1.0
                            b_tri[c, jj] = 1.0
                            c_tri[c, jj] = 0.0
                            r_tri[c, jj] = 0.0
                        end
                    end
                end

                # Subtract initial concentration and source for tendency calculation
                for c in bounds
                    mask_bgc_soilc[c] || continue
                    for j in 1:nlevdecomp
                        if !use_soil_matrixcn
                            trcr_tendency_ptr[c, j, s] = 0.0 - (conc_trcr[c, j + 1] + source[c, j, s])
                        else
                            trcr_tendency_ptr[c, j, s] = 0.0
                        end
                    end
                end

                if !use_soil_matrixcn
                    # Solve tridiagonal system for each column
                    # The Fortran Tridiagonal uses levels 0:nlevdecomp+1
                    # In our Julia mapping, these are indices 1:nlevdecomp+2 in the work arrays
                    nlev_tri = nlevdecomp + 2  # total levels (0 through nlevdecomp+1)

                    for c in bounds
                        mask_bgc_soilc[c] || continue

                        # Extract column vectors for the tridiagonal solve
                        a_col = a_tri[c, 1:nlev_tri]
                        b_col = b_tri[c, 1:nlev_tri]
                        c_col = c_tri[c, 1:nlev_tri]
                        r_col = r_tri[c, 1:nlev_tri]
                        u_col = zeros(FT, nlev_tri)

                        # Solve using the single-column tridiagonal solver
                        # jtop = 1 (index 1 = Fortran level 0)
                        tridiagonal_solve!(u_col, a_col, b_col, c_col, r_col, 1, nlev_tri)

                        # Copy solution back to conc_trcr
                        for jj in 1:nlev_tri
                            conc_trcr[c, jj] = u_col[jj]
                        end
                    end

                    # Add post-transport concentration to calculate tendency
                    for c in bounds
                        mask_bgc_soilc[c] || continue
                        for j in 1:nlevdecomp
                            trcr_tendency_ptr[c, j, s] = trcr_tendency_ptr[c, j, s] + conc_trcr[c, j + 1]
                            trcr_tendency_ptr[c, j, s] = trcr_tendency_ptr[c, j, s] / dtime
                        end
                    end
                end  # !use_soil_matrixcn

            else
                # CWD pools: just add source, no vertical transport
                for j in 1:nlevdecomp
                    for c in bounds
                        mask_bgc_soilc[c] || continue
                        if !use_soil_matrixcn
                            conc_trcr[c, j + 1] = conc_ptr[c, j, s] + source[c, j, s]
                        end
                    end
                end
            end  # !is_cwd

            if !use_soil_matrixcn
                # Write concentrations back to state, correct for bedrock leakage
                for j in 1:nlevdecomp
                    for c in bounds
                        mask_bgc_soilc[c] || continue
                        conc_ptr[c, j, s] = conc_trcr[c, j + 1]
                        # Correct for small amounts of tracer that leak into bedrock
                        if j > col.nbedrock[c]
                            conc_ptr[c, col.nbedrock[c], s] = conc_ptr[c, col.nbedrock[c], s] +
                                conc_trcr[c, j + 1] * (dzsoi_decomp_vals[j] / dzsoi_decomp_vals[col.nbedrock[c]])
                            conc_ptr[c, j, s] = 0.0
                        end
                    end
                end
            end  # !use_soil_matrixcn
        end  # s (pool loop)
    end  # i_type

    return nothing
end
