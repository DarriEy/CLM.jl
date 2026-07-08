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
@inline function patankar_A(pe::T) where {T<:Real}
    # eltype-generic so it lowers to valid Metal IR (no Float64 literals); on Float64
    # this is byte-identical (T(0.1)===0.1, etc.).
    return max(zero(T), (one(T) - T(0.1) * abs(pe))^5)
end

# ===========================================================================
# GPU kernelization. The only layer-coupled step is the per-column tridiagonal
# solve (Thomas), which a single column thread runs sequentially on its own row
# — so the WHOLE routine is two per-column kernels: a one-shot coefficient kernel,
# and a transport kernel launched per (tracer-type, pool) that does the diffusivity/
# Peclet build, tridiagonal assembly, in-thread Thomas solve, post-transport
# tendency, and write-back (incl. the bedrock-leak accumulation into [c,nbedrock]).
# Every write targets the thread's own column row → race-free, byte-identical.
# Thomas is inlined to match tridiagonal_solve! EXACTLY (no singularity guard, unlike
# tridiagonal_multi!) so the result is byte-identical to the host loop.
# ===========================================================================

# Per-column scratch bundle (all [ncols × *] matrices, disjoint rows per thread).
Base.@kwdef struct _LvtScr{M}
    diffus::M; adv_flux::M; conc_trcr::M
    d_p1_zp1::M; d_m1_zm1::M; f_p1::M; f_m1::M; pe_p1::M; pe_m1::M
    a_tri::M; b_tri::M; c_tri::M; r_tri::M; cp::M; dp::M
end
Adapt.@adapt_structure _LvtScr

# Diffusivity / advection coefficients (per column; internal j-loop, own-index).
@kernel function _lvt_coef_kernel!(som_diffus_coef, som_adv_coef, @Const(mask),
        @Const(altmax), @Const(altmax_lastyear), @Const(nbedrock), @Const(zisoi),
        nlevdecomp::Int, max_altdepth_cryoturbation, cryoturb_diffusion_k,
        som_diffus_val, som_adv_flux_val, max_depth_cryoturb_val)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(som_diffus_coef)
        alt_max_val = max(altmax[c], altmax_lastyear[c])
        if alt_max_val <= max_altdepth_cryoturbation && alt_max_val > zero(T)
            for j in 1:(nlevdecomp + 1)
                if j <= nbedrock[c] + 1
                    if zisoi[j + 1] < alt_max_val
                        som_diffus_coef[c, j] = cryoturb_diffusion_k
                        som_adv_coef[c, j] = zero(T)
                    else
                        som_diffus_coef[c, j] = max(
                            cryoturb_diffusion_k *
                            (one(T) - (zisoi[j + 1] - alt_max_val) /
                            (min(max_depth_cryoturb_val, zisoi[nbedrock[c] + 2]) - alt_max_val)),
                            zero(T))
                        som_adv_coef[c, j] = zero(T)
                    end
                else
                    som_adv_coef[c, j] = zero(T)
                    som_diffus_coef[c, j] = zero(T)
                end
            end
        elseif alt_max_val > zero(T)
            for j in 1:(nlevdecomp + 1)
                if j <= nbedrock[c] + 1
                    som_adv_coef[c, j] = som_adv_flux_val
                    som_diffus_coef[c, j] = som_diffus_val
                else
                    som_adv_coef[c, j] = zero(T)
                    som_diffus_coef[c, j] = zero(T)
                end
            end
        else
            for j in 1:(nlevdecomp + 1)
                som_adv_coef[c, j] = zero(T)
                som_diffus_coef[c, j] = zero(T)
            end
        end
    end
end

# Per-column transport for one (tracer-type, pool s). conc_ptr/source/trcr_tendency
# are the C- or N- arrays selected on the host; `s` indexes the pool, `is_cwd_s` the
# per-pool CWD flag, `spinup_factor_s` the per-pool spinup factor.
@kernel function _lvt_column_kernel!(conc_ptr, source, trcr_tendency,
        @Const(som_diffus_coef), @Const(som_adv_coef), scr, tri_ma_vr,
        @Const(zsoi_ext), @Const(dz_node), @Const(zisoi), @Const(dzsoi_decomp),
        @Const(nbedrock), @Const(gridcell), @Const(latdeg), @Const(mask),
        nlevdecomp::Int, s::Int, is_cwd_s::Bool, use_soil_matrixcn::Bool,
        build_tri_ma::Bool, ndecomp_pools::Int,
        spinup_state::Int, spinup_factor_s, dtime, epsilon)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(conc_ptr)
        nbr = nbedrock[c]
      if use_soil_matrixcn
        # Matrix-CN path: the transport is solved by the matrix infrastructure; here
        # only the tendency is zeroed and the state is left untouched (mirrors the host).
        for j in 1:nlevdecomp
            trcr_tendency[c, j, s] = zero(T)
        end
        # Build the vertical-transport sparse matrix tri_ma_vr ONCE (first non-CWD pool
        # of the C pass): it is pool-independent (in the non-spinup live run), so the
        # column's tridiagonal coefficients are repackaged into every non-CWD pool block.
        # Port of SoilBiogeochemLittVertTransp:401-429. Only a_tri/b_tri/c_tri are needed
        # (not conc_trcr / r_tri / the Thomas solve).
        if build_tri_ma && !is_cwd_s
            stride = nlevdecomp * 3 - 2
            for k in 1:(stride * (ndecomp_pools - 1))
                tri_ma_vr[c, k] = zero(T)
            end
            spinup_term = one(T)
            if spinup_state >= 1
                spinup_term = spinup_factor_s
            end
            if abs(spinup_term - one(T)) > T(1.0e-6)
                spinup_term = spinup_term * _spinup_lat_term(latdeg[gridcell[c]], T)
            end
            for j in 1:(nlevdecomp + 1)
                scr.adv_flux[c, j] = abs(som_adv_coef[c, j]) * spinup_term < epsilon ?
                    epsilon : som_adv_coef[c, j] * spinup_term
                scr.diffus[c, j] = abs(som_diffus_coef[c, j]) * spinup_term < epsilon ?
                    epsilon : som_diffus_coef[c, j] * spinup_term
            end
            # D/F/Pe throughout the column (identical to the sequential branch, minus
            # the conc_trcr writes which the coefficients do not use).
            for j in 1:(nlevdecomp + 1)
                if j == 1
                    scr.d_m1_zm1[c, j] = zero(T)
                    w_p1 = (zsoi_ext[j + 1] - zisoi[j + 1]) / dz_node[j + 1]
                    if scr.diffus[c, j + 1] > zero(T) && scr.diffus[c, j] > zero(T)
                        d_p1 = one(T) / ((one(T) - w_p1) / scr.diffus[c, j] + w_p1 / scr.diffus[c, j + 1])
                    else
                        d_p1 = zero(T)
                    end
                    scr.d_p1_zp1[c, j] = d_p1 / dz_node[j + 1]
                    scr.f_m1[c, j] = scr.adv_flux[c, j]
                    scr.f_p1[c, j] = scr.adv_flux[c, j + 1]
                    scr.pe_m1[c, j] = zero(T)
                    scr.pe_p1[c, j] = scr.f_p1[c, j] / scr.d_p1_zp1[c, j]
                elseif j >= nbr + 1
                    w_m1 = (zisoi[j] - zsoi_ext[j - 1]) / dz_node[j]
                    if scr.diffus[c, j] > zero(T) && scr.diffus[c, j - 1] > zero(T)
                        d_m1 = one(T) / ((one(T) - w_m1) / scr.diffus[c, j] + w_m1 / scr.diffus[c, j - 1])
                    else
                        d_m1 = zero(T)
                    end
                    scr.d_m1_zm1[c, j] = d_m1 / dz_node[j]
                    scr.d_p1_zp1[c, j] = scr.d_m1_zm1[c, j]
                    scr.f_m1[c, j] = scr.adv_flux[c, j]
                    scr.f_p1[c, j] = zero(T)
                    scr.pe_m1[c, j] = scr.f_m1[c, j] / scr.d_m1_zm1[c, j]
                    scr.pe_p1[c, j] = scr.f_p1[c, j] / scr.d_p1_zp1[c, j]
                else
                    w_m1 = (zisoi[j] - zsoi_ext[j - 1]) / dz_node[j]
                    if scr.diffus[c, j - 1] > zero(T) && scr.diffus[c, j] > zero(T)
                        d_m1 = one(T) / ((one(T) - w_m1) / scr.diffus[c, j] + w_m1 / scr.diffus[c, j - 1])
                    else
                        d_m1 = zero(T)
                    end
                    w_p1 = (zsoi_ext[j + 1] - zisoi[j + 1]) / dz_node[j + 1]
                    if scr.diffus[c, j + 1] > zero(T) && scr.diffus[c, j] > zero(T)
                        d_p1 = one(T) / ((one(T) - w_p1) / scr.diffus[c, j] + w_p1 / scr.diffus[c, j + 1])
                    else
                        d_p1 = (one(T) - w_m1) * scr.diffus[c, j] + w_p1 * scr.diffus[c, j + 1]
                    end
                    scr.d_m1_zm1[c, j] = d_m1 / dz_node[j]
                    scr.d_p1_zp1[c, j] = d_p1 / dz_node[j + 1]
                    scr.f_m1[c, j] = scr.adv_flux[c, j]
                    scr.f_p1[c, j] = scr.adv_flux[c, j + 1]
                    scr.pe_m1[c, j] = scr.f_m1[c, j] / scr.d_m1_zm1[c, j]
                    scr.pe_p1[c, j] = scr.f_p1[c, j] / scr.d_p1_zp1[c, j]
                end
            end
            # Repackage a_tri/b_tri/c_tri into tri_ma_vr for each non-CWD pool block.
            for j in 1:nlevdecomp
                a_p_0 = dzsoi_decomp[j] / dtime
                at = -(scr.d_m1_zm1[c, j] * patankar_A(scr.pe_m1[c, j]) + max(scr.f_m1[c, j], zero(T)))
                ct = -(scr.d_p1_zp1[c, j] * patankar_A(scr.pe_p1[c, j]) + max(-scr.f_p1[c, j], zero(T)))
                bt = -at - ct + a_p_0
                if j == 1
                    for i in 1:(ndecomp_pools - 1)
                        base = (i - 1) * stride
                        tri_ma_vr[c, base + 1] = (bt - a_p_0) / dzsoi_decomp[j] * (-dtime)
                        tri_ma_vr[c, base + 3] = ct / dzsoi_decomp[j] * (-dtime)
                    end
                elseif j <= nbr
                    for i in 1:(ndecomp_pools - 1)
                        base = (i - 1) * stride
                        tri_ma_vr[c, base + j * 3 - 4] = at / dzsoi_decomp[j] * (-dtime)
                        if j != nlevdecomp
                            tri_ma_vr[c, base + j * 3] = ct / dzsoi_decomp[j] * (-dtime)
                        end
                        tri_ma_vr[c, base + j * 3 - 2] = (bt - a_p_0) / dzsoi_decomp[j] * (-dtime)
                    end
                elseif j == nbr + 1 && j != nlevdecomp && j > 1
                    for i in 1:(ndecomp_pools - 1)
                        base = (i - 1) * stride
                        tri_ma_vr[c, base + (j - 1) * 3 - 2] =
                            tri_ma_vr[c, base + (j - 1) * 3 - 2] + at / dzsoi_decomp[j - 1] * (-dtime)
                    end
                end
            end
        end
      else
        if is_cwd_s
            # CWD pools: no transport, just add source.
            for j in 1:nlevdecomp
                scr.conc_trcr[c, j + 1] = conc_ptr[c, j, s] + source[c, j, s]
            end
        else
            spinup_term = one(T)
            if spinup_state >= 1
                spinup_term = spinup_factor_s
            end
            if abs(spinup_term - one(T)) > T(1.0e-6)
                spinup_term = spinup_term * _spinup_lat_term(latdeg[gridcell[c]], T)
            end

            for j in 1:(nlevdecomp + 1)
                if abs(som_adv_coef[c, j]) * spinup_term < epsilon
                    scr.adv_flux[c, j] = epsilon
                else
                    scr.adv_flux[c, j] = som_adv_coef[c, j] * spinup_term
                end
                if abs(som_diffus_coef[c, j]) * spinup_term < epsilon
                    scr.diffus[c, j] = epsilon
                else
                    scr.diffus[c, j] = som_diffus_coef[c, j] * spinup_term
                end
            end

            # Boundary values for conc_trcr (j=0 -> index 1; below bedrock -> 0).
            scr.conc_trcr[c, 1] = zero(T)
            for jj in (nbr + 1):(nlevdecomp + 1)
                scr.conc_trcr[c, jj + 1] = zero(T)
            end

            # conc_trcr + D/dz, F, Pe throughout the column.
            for j in 1:(nlevdecomp + 1)
                if j <= nlevdecomp
                    scr.conc_trcr[c, j + 1] = conc_ptr[c, j, s]
                end
                if j == 1
                    scr.d_m1_zm1[c, j] = zero(T)
                    w_p1 = (zsoi_ext[j + 1] - zisoi[j + 1]) / dz_node[j + 1]
                    if scr.diffus[c, j + 1] > zero(T) && scr.diffus[c, j] > zero(T)
                        d_p1 = one(T) / ((one(T) - w_p1) / scr.diffus[c, j] + w_p1 / scr.diffus[c, j + 1])
                    else
                        d_p1 = zero(T)
                    end
                    scr.d_p1_zp1[c, j] = d_p1 / dz_node[j + 1]
                    scr.f_m1[c, j] = scr.adv_flux[c, j]
                    scr.f_p1[c, j] = scr.adv_flux[c, j + 1]
                    scr.pe_m1[c, j] = zero(T)
                    scr.pe_p1[c, j] = scr.f_p1[c, j] / scr.d_p1_zp1[c, j]
                elseif j >= nbr + 1
                    w_m1 = (zisoi[j] - zsoi_ext[j - 1]) / dz_node[j]
                    if scr.diffus[c, j] > zero(T) && scr.diffus[c, j - 1] > zero(T)
                        d_m1 = one(T) / ((one(T) - w_m1) / scr.diffus[c, j] + w_m1 / scr.diffus[c, j - 1])
                    else
                        d_m1 = zero(T)
                    end
                    scr.d_m1_zm1[c, j] = d_m1 / dz_node[j]
                    scr.d_p1_zp1[c, j] = scr.d_m1_zm1[c, j]
                    scr.f_m1[c, j] = scr.adv_flux[c, j]
                    scr.f_p1[c, j] = zero(T)
                    scr.pe_m1[c, j] = scr.f_m1[c, j] / scr.d_m1_zm1[c, j]
                    scr.pe_p1[c, j] = scr.f_p1[c, j] / scr.d_p1_zp1[c, j]
                else
                    w_m1 = (zisoi[j] - zsoi_ext[j - 1]) / dz_node[j]
                    if scr.diffus[c, j - 1] > zero(T) && scr.diffus[c, j] > zero(T)
                        d_m1 = one(T) / ((one(T) - w_m1) / scr.diffus[c, j] + w_m1 / scr.diffus[c, j - 1])
                    else
                        d_m1 = zero(T)
                    end
                    w_p1 = (zsoi_ext[j + 1] - zisoi[j + 1]) / dz_node[j + 1]
                    if scr.diffus[c, j + 1] > zero(T) && scr.diffus[c, j] > zero(T)
                        d_p1 = one(T) / ((one(T) - w_p1) / scr.diffus[c, j] + w_p1 / scr.diffus[c, j + 1])
                    else
                        d_p1 = (one(T) - w_m1) * scr.diffus[c, j] + w_p1 * scr.diffus[c, j + 1]
                    end
                    scr.d_m1_zm1[c, j] = d_m1 / dz_node[j]
                    scr.d_p1_zp1[c, j] = d_p1 / dz_node[j + 1]
                    scr.f_m1[c, j] = scr.adv_flux[c, j]
                    scr.f_p1[c, j] = scr.adv_flux[c, j + 1]
                    scr.pe_m1[c, j] = scr.f_m1[c, j] / scr.d_m1_zm1[c, j]
                    scr.pe_p1[c, j] = scr.f_p1[c, j] / scr.d_p1_zp1[c, j]
                end
            end

            # Tridiagonal coefficients (Fortran j=0..nlevdecomp+1 -> jj=j+1).
            for j in 0:(nlevdecomp + 1)
                jj = j + 1
                a_p_0 = (j > 0 && j < nlevdecomp + 1) ? dzsoi_decomp[j] / dtime : zero(T)
                if j == 0
                    scr.a_tri[c, jj] = zero(T)
                    scr.b_tri[c, jj] = one(T)
                    scr.c_tri[c, jj] = -one(T)
                    scr.r_tri[c, jj] = zero(T)
                elseif j == 1
                    scr.a_tri[c, jj] = -(scr.d_m1_zm1[c, j] * patankar_A(scr.pe_m1[c, j]) + max(scr.f_m1[c, j], zero(T)))
                    scr.c_tri[c, jj] = -(scr.d_p1_zp1[c, j] * patankar_A(scr.pe_p1[c, j]) + max(-scr.f_p1[c, j], zero(T)))
                    scr.b_tri[c, jj] = -scr.a_tri[c, jj] - scr.c_tri[c, jj] + a_p_0
                    scr.r_tri[c, jj] = source[c, j, s] * dzsoi_decomp[j] / dtime +
                                       (a_p_0 - scr.adv_flux[c, j]) * scr.conc_trcr[c, jj]
                elseif j < nlevdecomp + 1
                    scr.a_tri[c, jj] = -(scr.d_m1_zm1[c, j] * patankar_A(scr.pe_m1[c, j]) + max(scr.f_m1[c, j], zero(T)))
                    scr.c_tri[c, jj] = -(scr.d_p1_zp1[c, j] * patankar_A(scr.pe_p1[c, j]) + max(-scr.f_p1[c, j], zero(T)))
                    scr.b_tri[c, jj] = -scr.a_tri[c, jj] - scr.c_tri[c, jj] + a_p_0
                    scr.r_tri[c, jj] = source[c, j, s] * dzsoi_decomp[j] / dtime +
                                       a_p_0 * scr.conc_trcr[c, jj]
                else
                    scr.a_tri[c, jj] = -one(T)
                    scr.b_tri[c, jj] = one(T)
                    scr.c_tri[c, jj] = zero(T)
                    scr.r_tri[c, jj] = zero(T)
                end
            end

            # Tendency: subtract initial concentration + source (pre-solve).
            for j in 1:nlevdecomp
                trcr_tendency[c, j, s] = zero(T) - (scr.conc_trcr[c, j + 1] + source[c, j, s])
            end

            # In-thread Thomas solve over jj=1..nlevdecomp+2 (jtop=1), matching
            # tridiagonal_solve! exactly (no singularity guard) for byte-identity.
            nlev_tri = nlevdecomp + 2
            scr.cp[c, 1] = scr.c_tri[c, 1] / scr.b_tri[c, 1]
            scr.dp[c, 1] = scr.r_tri[c, 1] / scr.b_tri[c, 1]
            for jj in 2:nlev_tri
                denom = scr.b_tri[c, jj] - scr.a_tri[c, jj] * scr.cp[c, jj - 1]
                scr.cp[c, jj] = scr.c_tri[c, jj] / denom
                scr.dp[c, jj] = (scr.r_tri[c, jj] - scr.a_tri[c, jj] * scr.dp[c, jj - 1]) / denom
            end
            scr.conc_trcr[c, nlev_tri] = scr.dp[c, nlev_tri]
            for jj in (nlev_tri - 1):-1:1
                scr.conc_trcr[c, jj] = scr.dp[c, jj] - scr.cp[c, jj] * scr.conc_trcr[c, jj + 1]
            end

            # Post-transport tendency.
            for j in 1:nlevdecomp
                trcr_tendency[c, j, s] = (trcr_tendency[c, j, s] + scr.conc_trcr[c, j + 1]) / dtime
            end
        end  # is_cwd_s

        # Write concentrations back; correct for tracer leaking into bedrock.
        for j in 1:nlevdecomp
            conc_ptr[c, j, s] = scr.conc_trcr[c, j + 1]
            if j > nbr
                conc_ptr[c, nbr, s] = conc_ptr[c, nbr, s] +
                    scr.conc_trcr[c, j + 1] * (dzsoi_decomp[j] / dzsoi_decomp[nbr])
                conc_ptr[c, j, s] = zero(T)
            end
        end
      end  # use_soil_matrixcn
    end
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
- `mask_bgc_soilc::AbstractVector{Bool}`            : soil column mask (in)
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
        mask_bgc_soilc::AbstractVector{Bool},
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
        max_depth_cryoturb_val::Real = 3.0,
        # Carbon-isotope transport (optional). When use_c13/use_c14 is true the
        # corresponding isotope SoilBiogeochem carbon state/flux instances must be
        # supplied; their decomp_cpools are then transported alongside C and N
        # (extra i_type passes, matching the Fortran ntype = 2 + use_c13 + use_c14).
        use_c13::Bool = false,
        use_c14::Bool = false,
        c13_cs::Union{SoilBiogeochemCarbonStateData, Nothing} = nothing,
        c13_cf::Union{SoilBiogeochemCarbonFluxData, Nothing} = nothing,
        c14_cs::Union{SoilBiogeochemCarbonStateData, Nothing} = nothing,
        c14_cf::Union{SoilBiogeochemCarbonFluxData, Nothing} = nothing)

    # Unpack cascade configuration. is_cwd (BitVector) + spinup_factor (host vector)
    # are indexed per-pool ON THE HOST and passed as scalars to the column kernel, so
    # they need no device move.
    is_cwd        = cascade_con.is_cwd
    spinup_factor = cascade_con.spinup_factor

    # Unpack SoilBiogeochemState output arrays + active-layer / column data (device-resident).
    altmax          = active_layer.altmax_col
    altmax_lastyear = active_layer.altmax_lastyear_col
    som_adv_coef    = st.som_adv_coef_col
    som_diffus_coef = st.som_diffus_coef_col

    nc = length(bounds)
    # i_type passes: always C (1) and N (2); + C13 and/or C14 when enabled, matching
    # Fortran ntype = 2 + use_c13 + use_c14. Build the per-pass (state, sourcesink,
    # transport_tendency) array triples in Fortran's case order (C13 before C14).
    if use_c13 && (c13_cs === nothing || c13_cf === nothing)
        error("litter_vert_transp!: use_c13=true requires c13_cs and c13_cf")
    end
    if use_c14 && (c14_cs === nothing || c14_cf === nothing)
        error("litter_vert_transp!: use_c14=true requires c14_cs and c14_cf")
    end
    type_arrays = Vector{NTuple{3, Any}}()
    push!(type_arrays, (cs.decomp_cpools_vr_col,
                        cf.decomp_cpools_sourcesink_col,
                        cf.decomp_cpools_transport_tendency_col))
    push!(type_arrays, (ns.decomp_npools_vr_col,
                        nf.decomp_npools_sourcesink_col,
                        nf.decomp_npools_transport_tendency_col))
    if use_c13
        push!(type_arrays, (c13_cs.decomp_cpools_vr_col,
                            c13_cf.decomp_cpools_sourcesink_col,
                            c13_cf.decomp_cpools_transport_tendency_col))
    end
    if use_c14
        push!(type_arrays, (c14_cs.decomp_cpools_vr_col,
                            c14_cf.decomp_cpools_sourcesink_col,
                            c14_cf.decomp_cpools_transport_tendency_col))
    end
    ntype = length(type_arrays)

    FT  = eltype(cs.decomp_cpools_vr_col)
    ref = cs.decomp_cpools_vr_col          # backend prototype for device scratch
    epsilon = FT(1.0e-30)

    # Working-precision scalar params (no Float64 reaches a Metal kernel).
    som_diffus_val_ft        = FT(params.som_diffus)
    cryoturb_diffusion_k_ft  = FT(params.cryoturb_diffusion_k)
    max_altdepth_cryo_ft     = FT(params.max_altdepth_cryoturbation)
    som_adv_flux_ft          = FT(som_adv_flux_val)
    max_depth_cryoturb_ft    = FT(max_depth_cryoturb_val)
    dtime_ft                 = FT(dtime)

    # Move the host node-geometry / interface vectors onto the state backend (a host
    # Vector is a non-bitstype kernel arg on Metal). zsoi_ext/dz_node are built on the
    # host from the (host) zsoi/zisoi args, then moved.
    _dev(v) = (d = similar(ref, FT, length(v)); copyto!(d, collect(FT, v)); d)
    zsoi_ext_h = zeros(FT, nlevdecomp + 1)
    zsoi_ext_h[1:nlevdecomp] .= FT.(zsoi_vals[1:nlevdecomp])
    zsoi_ext_h[nlevdecomp + 1] = FT(zisoi_vals[nlevdecomp + 1])
    dz_node_h = zeros(FT, nlevdecomp + 1)
    dz_node_h[1] = zsoi_ext_h[1]
    for j in 2:(nlevdecomp + 1)
        dz_node_h[j] = zsoi_ext_h[j] - zsoi_ext_h[j - 1]
    end
    zsoi_ext     = _dev(zsoi_ext_h)
    dz_node      = _dev(dz_node_h)
    zisoi        = _dev(zisoi_vals)
    dzsoi_decomp = _dev(dzsoi_decomp_vals)

    # mask onto the backend as a plain Bool vector (BitVector is non-bitstype on device).
    mask = similar(ref, Bool, length(mask_bgc_soilc))
    copyto!(mask, collect(Bool, mask_bgc_soilc))

    # ------ Compute diffusivity / advection coefficients (one thread per column) ------
    _launch!(_lvt_coef_kernel!, som_diffus_coef, som_adv_coef, mask,
             altmax, altmax_lastyear, col.nbedrock, zisoi, nlevdecomp,
             max_altdepth_cryo_ft, cryoturb_diffusion_k_ft, som_diffus_val_ft,
             som_adv_flux_ft, max_depth_cryoturb_ft; ndrange = nc)

    # ------ Device-resident per-column scratch (Fortran j=0..nlevdecomp+1 -> Julia
    # index j+1; a_tri/conc_trcr/cp/dp span nj = nlevdecomp+2; the D/F/Pe arrays span
    # nlevdecomp+1). Allocated like the state arrays so they live on its backend. ------
    nj   = nlevdecomp + 2
    n1   = nlevdecomp + 1
    ncols = length(mask_bgc_soilc)
    _m(w) = fill!(similar(ref, FT, ncols, w), zero(FT))
    scr = _LvtScr(; diffus = _m(n1), adv_flux = _m(n1), conc_trcr = _m(nj),
                    d_p1_zp1 = _m(n1), d_m1_zm1 = _m(n1), f_p1 = _m(n1), f_m1 = _m(n1),
                    pe_p1 = _m(n1), pe_m1 = _m(n1),
                    a_tri = _m(nj), b_tri = _m(nj), c_tri = _m(nj), r_tri = _m(nj),
                    cp = _m(nj), dp = _m(nj))

    # tri_ma_vr must be sized (ncols, Ntri_setup) for the matrix build to write it;
    # size it here if the flux instance was not initialised with use_soil_matrixcn.
    if use_soil_matrixcn
        Ntri = (nlevdecomp * 3 - 2) * (ndecomp_pools - 1)
        if size(cf.tri_ma_vr, 1) != ncols || size(cf.tri_ma_vr, 2) != Ntri
            cf.tri_ma_vr = fill!(similar(ref, FT, ncols, Ntri), zero(FT))
        end
    end

    # ------ Loop over tracer types (C, N) ------
    for i_type in 1:ntype

        # Select the appropriate state/flux arrays for this pass
        # (1=C, 2=N, then C13 and/or C14 when enabled).
        conc_ptr, source, trcr_tendency_ptr = type_arrays[i_type]

        for s in 1:ndecomp_pools
            # is_cwd[s] / spinup_factor[s] are read on the host and passed as scalars.
            spinup_factor_s = FT(spinup_factor[s])
            # Build the pool-independent vertical-transport matrix tri_ma_vr once, on the
            # first (C, non-CWD) pool pass (Fortran: s==1 .and. i_type==1).
            build_tri_ma = use_soil_matrixcn && i_type == 1 && s == 1
            _launch!(_lvt_column_kernel!, conc_ptr, source, trcr_tendency_ptr,
                     som_diffus_coef, som_adv_coef, scr, cf.tri_ma_vr,
                     zsoi_ext, dz_node, zisoi, dzsoi_decomp,
                     col.nbedrock, col.gridcell, grc.latdeg, mask,
                     nlevdecomp, s, is_cwd[s], use_soil_matrixcn,
                     build_tri_ma, ndecomp_pools, spinup_state,
                     spinup_factor_s, dtime_ft, epsilon; ndrange = nc)
        end  # s (pool loop)
    end  # i_type

    return nothing
end
