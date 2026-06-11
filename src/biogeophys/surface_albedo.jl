# ==========================================================================
# Ported from: src/biogeophys/SurfaceAlbedoMod.F90
# Performs surface albedo calculations
#
# Public functions:
#   surface_albedo_init_time_const!  — initialize time-constant albedo data
#   soil_albedo!                     — ground surface albedo
#   two_stream!                      — two-stream canopy radiative transfer
#   surface_albedo!                  — top-level surface albedo driver
# ==========================================================================

# ---- Module-level constants and data ----

# albedo land ice by waveband (1=vis, 2=nir)
const ALBICE = [0.80, 0.55]

# albedo frozen lakes by waveband (1=vis, 2=nir)
const ALBLAK = [0.60, 0.40]

# Coefficient for calculating ice "fraction" for lake surface albedo
# From D. Mironov (2010) Boreal Env. Research
const CALB = 95.6

# prevents overflow for division by zero
const MPE_ALBEDO = 1.0e-06

"""
    SurfaceAlbedoConstants

Module-level data that is initialized once and remains constant during
the simulation. Corresponds to Fortran module variables `albsat`, `albdry`,
`isoicol`, `alblakwi`, and `snowveg_affects_radiation`.
"""
Base.@kwdef mutable struct SurfaceAlbedoConstants{FT<:Real}
    albsat::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # wet soil albedo by color class and waveband
    albdry::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # dry soil albedo by color class and waveband
    isoicol::Vector{Int} = Int[]                              # column soil color class
    alblakwi::Vector{Float64} = [0.10, 0.10]                 # albedo of melting lakes (namelist-settable)
    lake_melt_icealb::Vector{Float64} = [0.10, 0.10]         # namelist default for alblakwi
    snowveg_affects_radiation::Bool = true                    # whether canopy snow affects radiation
end

const surfalb_con = SurfaceAlbedoConstants()

# --------------------------------------------------------------------------
# KernelAbstractions kernels for simple, iteration-independent per-element
# loops in this file. Each replaces an inline loop with identical semantics.
# Assumes single-clump indexing (begc==1/begp==1) so output array length
# equals the iteration range, consistent with src/infrastructure/kernels.jl.
# --------------------------------------------------------------------------

# coszen_col_arr[c] = coszen_grc[gridcell[c]]  (both downscale branches are
# identical no-ops, so the result is unconditional)
@kernel function _surfalb_coszen_col_kernel!(coszen_col_arr, @Const(gridcell), @Const(coszen_grc))
    c = @index(Global)
    @inbounds coszen_col_arr[c] = coszen_grc[gridcell[c]]
end

surfalb_coszen_col!(coszen_col_arr, gridcell, coszen_grc) =
    _launch!(_surfalb_coszen_col_kernel!, coszen_col_arr, gridcell, coszen_grc)

# coszen_patch_arr[p] = coszen_col_arr[column[p]]  (masked per patch)
@kernel function _surfalb_coszen_patch_kernel!(coszen_patch_arr, @Const(mask),
                                               @Const(column), @Const(coszen_col_arr))
    p = @index(Global)
    @inbounds if mask[p]
        coszen_patch_arr[p] = coszen_col_arr[column[p]]
    end
end

surfalb_coszen_patch!(coszen_patch_arr, mask, column, coszen_col_arr) =
    _launch!(_surfalb_coszen_patch_kernel!, coszen_patch_arr, mask, column, coszen_col_arr)

# Leaf/stem weighting fractions wl, ws_arr (masked per patch)
@kernel function _surfalb_lai_weight_kernel!(wl, ws_arr, @Const(mask),
                                             @Const(elai), @Const(esai), mpe)
    p = @index(Global)
    @inbounds if mask[p]
        denom = smooth_max(elai[p] + esai[p], mpe)
        wl[p] = elai[p] / denom
        ws_arr[p] = esai[p] / denom
    end
end

surfalb_lai_weight!(wl, ws_arr, mask, elai, esai, mpe) =
    _launch!(_surfalb_lai_weight_kernel!, wl, ws_arr, mask, elai, esai, mpe; ndrange = length(wl))

# rho/tau weighted reflectance/transmittance (2D: patch x waveband, masked per patch)
@kernel function _surfalb_rho_tau_kernel!(rho, tau, @Const(mask), @Const(itype_p),
                                          @Const(wl), @Const(ws_arr),
                                          @Const(pftcon_rhol), @Const(pftcon_rhos),
                                          @Const(pftcon_taul), @Const(pftcon_taus), mpe)
    p, ib = @index(Global, NTuple)
    @inbounds if mask[p]
        itype = itype_p[p] + 1  # Fortran 0-based PFT -> Julia 1-based
        rho[p, ib] = smooth_max(pftcon_rhol[itype, ib] * wl[p] + pftcon_rhos[itype, ib] * ws_arr[p], mpe)
        tau[p, ib] = smooth_max(pftcon_taul[itype, ib] * wl[p] + pftcon_taus[itype, ib] * ws_arr[p], mpe)
    end
end

function surfalb_rho_tau!(rho, tau, mask, itype_p, wl, ws_arr,
                          pftcon_rhol, pftcon_rhos, pftcon_taul, pftcon_taus, mpe, numrad::Int)
    _launch!(_surfalb_rho_tau_kernel!, rho, tau, mask, itype_p, wl, ws_arr,
             pftcon_rhol, pftcon_rhos, pftcon_taul, pftcon_taus, mpe;
             ndrange = (size(rho, 1), numrad))
end

# --------------------------------------------------------------------------
# Orchestrator (surface_albedo!) per-element kernels.
# Each replaces an inline host loop with identical semantics. Field-array
# bundles (immutable, isbits-friendly + @adapt_structure) group the many
# output arrays a kernel touches so the call stays under the arg limit.
# --------------------------------------------------------------------------

# ---- Output-initialization (per (ib, col)) ----
Base.@kwdef struct _SaInitColOut{M}
    albsod::M; albsoi::M; albgrd::M; albgri::M
    albgrd_pur::M; albgri_pur::M; albgrd_bc::M; albgri_bc::M
    albgrd_oc::M; albgri_oc::M; albgrd_dst::M; albgri_dst::M
    albgrd_hst::M; albgri_hst::M; albgrd_pur_hst::M; albgri_pur_hst::M
    albgrd_bc_hst::M; albgri_bc_hst::M; albgrd_oc_hst::M; albgri_oc_hst::M
    albgrd_dst_hst::M; albgri_dst_hst::M; albsnd_hst2::M; albsni_hst2::M
    flx_absdv::M; flx_absdn::M; flx_absiv::M; flx_absin::M
end
Adapt.@adapt_structure _SaInitColOut

@kernel function _sa_init_col_kernel!(o::_SaInitColOut, @Const(mask), nflx::Int)
    c, ib = @index(Global, NTuple)
    @inbounds if mask[c]
        o.albsod[c, ib] = 0.0
        o.albsoi[c, ib] = 0.0
        o.albgrd[c, ib] = 0.0
        o.albgri[c, ib] = 0.0
        o.albgrd_pur[c, ib] = 0.0
        o.albgri_pur[c, ib] = 0.0
        o.albgrd_bc[c, ib] = 0.0
        o.albgri_bc[c, ib] = 0.0
        o.albgrd_oc[c, ib] = 0.0
        o.albgri_oc[c, ib] = 0.0
        o.albgrd_dst[c, ib] = 0.0
        o.albgri_dst[c, ib] = 0.0
        o.albgrd_hst[c, ib] = SPVAL
        o.albgri_hst[c, ib] = SPVAL
        o.albgrd_pur_hst[c, ib] = SPVAL
        o.albgri_pur_hst[c, ib] = SPVAL
        o.albgrd_bc_hst[c, ib] = SPVAL
        o.albgri_bc_hst[c, ib] = SPVAL
        o.albgrd_oc_hst[c, ib] = SPVAL
        o.albgri_oc_hst[c, ib] = SPVAL
        o.albgrd_dst_hst[c, ib] = SPVAL
        o.albgri_dst_hst[c, ib] = SPVAL
        o.albsnd_hst2[c, ib] = SPVAL
        o.albsni_hst2[c, ib] = SPVAL
        # Zero the flux-absorption profiles (nflx layers) once (on ib==1).
        if ib == 1
            for j in 1:nflx
                o.flx_absdv[c, j] = 0.0
                o.flx_absdn[c, j] = 0.0
                o.flx_absiv[c, j] = 0.0
                o.flx_absin[c, j] = 0.0
            end
        end
    end
end

# ---- Output-initialization (per (ib, patch)) ----
Base.@kwdef struct _SaInitPatchOut{M}
    albd::M; albi::M; albd_hst::M; albi_hst::M
    albdSF::M; albiSF::M
    fabd::M; fabd_sun::M; fabd_sha::M
    fabi::M; fabi_sun::M; fabi_sha::M
    ftdd::M; ftid::M; ftii::M
end
Adapt.@adapt_structure _SaInitPatchOut

@kernel function _sa_init_patch_kernel!(o::_SaInitPatchOut, @Const(mask), use_SSRE::Bool)
    p, ib = @index(Global, NTuple)
    @inbounds if mask[p]
        o.albd[p, ib] = 1.0
        o.albi[p, ib] = 1.0
        o.albd_hst[p, ib] = SPVAL
        o.albi_hst[p, ib] = SPVAL
        if use_SSRE
            o.albdSF[p, ib] = 1.0
            o.albiSF[p, ib] = 1.0
        end
        o.fabd[p, ib] = 0.0
        o.fabd_sun[p, ib] = 0.0
        o.fabd_sha[p, ib] = 0.0
        o.fabi[p, ib] = 0.0
        o.fabi_sun[p, ib] = 0.0
        o.fabi_sha[p, ib] = 0.0
        o.ftdd[p, ib] = 0.0
        o.ftid[p, ib] = 0.0
        o.ftii[p, ib] = 0.0
    end
end

# ---- Fallback age-based snow albedo (per col) ----
@kernel function _sa_snow_fallback_kernel!(albsnd, albsni, @Const(mask),
        @Const(coszen_col), @Const(frac_sno), @Const(snow_persist), npersist::Int)
    c = @index(Global)
    T = eltype(albsnd)
    @inbounds if mask[c]
        if coszen_col[c] > 0.0 && frac_sno[c] > 0.0
            age_days = c <= npersist ? snow_persist[c] / T(86400.0) : zero(T)
            fage = one(T) - exp(-age_days / T(5.0))
            albsnd[c, 1] = clamp(T(0.95) * (one(T) - T(0.15) * fage), T(0.1), T(0.99))
            albsni[c, 1] = clamp(T(0.95) * (one(T) - T(0.15) * fage), T(0.1), T(0.99))
            albsnd[c, 2] = clamp(T(0.65) * (one(T) - T(0.50) * fage), T(0.1), T(0.99))
            albsni[c, 2] = clamp(T(0.65) * (one(T) - T(0.50) * fage), T(0.1), T(0.99))
        end
    end
end

# ---- Ground albedos: snow-fraction weighting (per (ib, col)) ----
Base.@kwdef struct _SaGrndOut{M}
    albgrd::M; albgri::M
    albgrd_bc::M; albgri_bc::M; albgrd_dst::M; albgri_dst::M
    albgrd_pur::M; albgri_pur::M; albgrd_oc::M; albgri_oc::M
end
Adapt.@adapt_structure _SaGrndOut
Base.@kwdef struct _SaGrndIn{M,V}
    albsod::M; albsoi::M; albsnd::M; albsni::M; frac_sno::V
end
Adapt.@adapt_structure _SaGrndIn

@kernel function _sa_grnd_kernel!(o::_SaGrndOut, in::_SaGrndIn, @Const(mask),
        @Const(coszen_col), use_snicar_frc::Bool, do_sno_oc::Bool)
    c, ib = @index(Global, NTuple)
    T = eltype(o.albgrd)
    @inbounds if mask[c]
        if coszen_col[c] > 0.0
            fs = in.frac_sno[c]
            o.albgrd[c, ib] = in.albsod[c, ib] * (one(T) - fs) + in.albsnd[c, ib] * fs
            o.albgri[c, ib] = in.albsoi[c, ib] * (one(T) - fs) + in.albsni[c, ib] * fs
            if use_snicar_frc
                o.albgrd_bc[c, ib] = in.albsod[c, ib] * (one(T) - fs)
                o.albgri_bc[c, ib] = in.albsoi[c, ib] * (one(T) - fs)
                o.albgrd_dst[c, ib] = in.albsod[c, ib] * (one(T) - fs)
                o.albgri_dst[c, ib] = in.albsoi[c, ib] * (one(T) - fs)
                o.albgrd_pur[c, ib] = in.albsod[c, ib] * (one(T) - fs)
                o.albgri_pur[c, ib] = in.albsoi[c, ib] * (one(T) - fs)
                if do_sno_oc
                    o.albgrd_oc[c, ib] = in.albsod[c, ib] * (one(T) - fs)
                    o.albgri_oc[c, ib] = in.albsoi[c, ib] * (one(T) - fs)
                end
            end
        end
    end
end

# ---- SNICAR flux-absorption mapping to surfalb fields (per col) ----
Base.@kwdef struct _SaFluxMap{M3,M,V,VI,VB}
    flx_absdv::M; flx_absdn::M; flx_absiv::M; flx_absin::M
    flx_absd_snw::M3; flx_absi_snw::M3
    albsnd::M; albsni::M; albsod::M; albsoi::M
    frac_sno::V; landunit::VI; lakpoi::VB
end
Adapt.@adapt_structure _SaFluxMap

@kernel function _sa_fluxmap_kernel!(o::_SaFluxMap, @Const(mask), nlev1::Int,
        use_subgrid_fluxes::Bool, ivis::Int, inir::Int)
    c = @index(Global)
    T = eltype(o.flx_absdv)
    @inbounds if mask[c]
        l = o.landunit[c]
        is_lake = o.lakpoi[l]
        for i in 1:nlev1
            if use_subgrid_fluxes && !is_lake
                o.flx_absdv[c, i] = o.flx_absd_snw[c, i, ivis]
                o.flx_absdn[c, i] = o.flx_absd_snw[c, i, inir]
                o.flx_absiv[c, i] = o.flx_absi_snw[c, i, ivis]
                o.flx_absin[c, i] = o.flx_absi_snw[c, i, inir]
            else
                fsnow = o.frac_sno[c]
                # flx_absdv (direct VIS)
                f_snw = o.flx_absd_snw[c, i, ivis]; a_snw = o.albsnd[c, ivis]; a_soil = o.albsod[c, ivis]
                if a_snw < 1.0
                    o.flx_absdv[c, i] = f_snw * fsnow + (one(T) - fsnow) * (one(T) - a_soil) * (f_snw / max(one(T) - a_snw, T(1.0e-6)))
                else
                    o.flx_absdv[c, i] = 0.0
                end
                # flx_absdn (direct NIR)
                f_snw = o.flx_absd_snw[c, i, inir]; a_snw = o.albsnd[c, inir]; a_soil = o.albsod[c, inir]
                if a_snw < 1.0
                    o.flx_absdn[c, i] = f_snw * fsnow + (one(T) - fsnow) * (one(T) - a_soil) * (f_snw / max(one(T) - a_snw, T(1.0e-6)))
                else
                    o.flx_absdn[c, i] = 0.0
                end
                # flx_absiv (diffuse VIS)
                f_snw = o.flx_absi_snw[c, i, ivis]; a_snw = o.albsni[c, ivis]; a_soil = o.albsoi[c, ivis]
                if a_snw < 1.0
                    o.flx_absiv[c, i] = f_snw * fsnow + (one(T) - fsnow) * (one(T) - a_soil) * (f_snw / max(one(T) - a_snw, T(1.0e-6)))
                else
                    o.flx_absiv[c, i] = 0.0
                end
                # flx_absin (diffuse NIR)
                f_snw = o.flx_absi_snw[c, i, inir]; a_snw = o.albsni[c, inir]; a_soil = o.albsoi[c, inir]
                if a_snw < 1.0
                    o.flx_absin[c, i] = f_snw * fsnow + (one(T) - fsnow) * (one(T) - a_soil) * (f_snw / max(one(T) - a_snw, T(1.0e-6)))
                else
                    o.flx_absin[c, i] = 0.0
                end
            end
        end
    end
end

# ---- Snow albedo history diagnostics (per (ib, col)) ----
@kernel function _sa_snowhist_kernel!(albsnd_hst, albsni_hst, @Const(mask),
        @Const(coszen_col), @Const(h2osno_total), @Const(albsnd), @Const(albsni))
    c, ib = @index(Global, NTuple)
    @inbounds if mask[c]
        if coszen_col[c] > 0.0 && h2osno_total[c] > 0.0
            albsnd_hst[c, ib] = albsnd[c, ib]
            albsni_hst[c, ib] = albsni[c, ib]
        else
            albsnd_hst[c, ib] = 0.0
            albsni_hst[c, ib] = 0.0
        end
    end
end

# ---- Canopy layering (per patch; internal canopy-layer iv loops) ----
Base.@kwdef struct _SaCanLayer{VI,V,M}
    nrad::VI; ncan::VI; tlai_z::M; tsai_z::M
    elai::V; esai::V; tlai::V; tsai::V
end
Adapt.@adapt_structure _SaCanLayer

@kernel function _sa_canlayer_kernel!(o::_SaCanLayer, @Const(mask),
        nlevcan::Int, dincmax, mpe)
    p = @index(Global)
    T = eltype(o.tlai_z)
    @inbounds if mask[p]
        dincmax_T = oftype(o.elai[p], dincmax)
        mpe_T = oftype(o.elai[p], mpe)
        if nlevcan == 1
            o.nrad[p] = 1
            o.ncan[p] = 1
            o.tlai_z[p, 1] = o.elai[p]
            o.tsai_z[p, 1] = o.esai[p]
        elseif nlevcan > 1
            if o.elai[p] + o.esai[p] == zero(T)
                o.nrad[p] = 0
            else
                dincmax_sum = zero(T)
                for iv in 1:nlevcan
                    dincmax_sum += dincmax_T
                    if ((o.elai[p] + o.esai[p]) - dincmax_sum) > T(1.0e-06)
                        o.nrad[p] = iv
                        dinc = dincmax_T
                        o.tlai_z[p, iv] = dinc * o.elai[p] / smooth_max(o.elai[p] + o.esai[p], mpe_T)
                        o.tsai_z[p, iv] = dinc * o.esai[p] / smooth_max(o.elai[p] + o.esai[p], mpe_T)
                    else
                        o.nrad[p] = iv
                        dinc = dincmax_T - (dincmax_sum - (o.elai[p] + o.esai[p]))
                        o.tlai_z[p, iv] = dinc * o.elai[p] / smooth_max(o.elai[p] + o.esai[p], mpe_T)
                        o.tsai_z[p, iv] = dinc * o.esai[p] / smooth_max(o.elai[p] + o.esai[p], mpe_T)
                        break
                    end
                end

                # Minimum of 4 canopy layers
                if o.nrad[p] < 4
                    o.nrad[p] = 4
                    for iv in 1:o.nrad[p]
                        o.tlai_z[p, iv] = o.elai[p] / o.nrad[p]
                        o.tsai_z[p, iv] = o.esai[p] / o.nrad[p]
                    end
                end
            end

            # Repeat for buried canopy layers
            blai = o.tlai[p] - o.elai[p]
            bsai = o.tsai[p] - o.esai[p]
            if blai + bsai == zero(T)
                o.ncan[p] = o.nrad[p]
            else
                dincmax_sum = zero(T)
                for iv in (o.nrad[p] + 1):nlevcan
                    dincmax_sum += dincmax_T
                    if ((blai + bsai) - dincmax_sum) > T(1.0e-06)
                        o.ncan[p] = iv
                        dinc = dincmax_T
                        o.tlai_z[p, iv] = dinc * blai / smooth_max(blai + bsai, mpe_T)
                        o.tsai_z[p, iv] = dinc * bsai / smooth_max(blai + bsai, mpe_T)
                    else
                        o.ncan[p] = iv
                        dinc = dincmax_T - (dincmax_sum - (blai + bsai))
                        o.tlai_z[p, iv] = dinc * blai / smooth_max(blai + bsai, mpe_T)
                        o.tsai_z[p, iv] = dinc * bsai / smooth_max(blai + bsai, mpe_T)
                        break
                    end
                end
            end
        end
    end
end

# ---- Zero fluxes for active canopy layers (per patch; iv loop) ----
@kernel function _sa_zerolayer_kernel!(fabd_sun_z, fabd_sha_z, fabi_sun_z, fabi_sha_z,
        fsun_z, @Const(mask), @Const(nrad))
    p = @index(Global)
    @inbounds if mask[p]
        for iv in 1:nrad[p]
            fabd_sun_z[p, iv] = 0.0
            fabd_sha_z[p, iv] = 0.0
            fabi_sun_z[p, iv] = 0.0
            fabi_sha_z[p, iv] = 0.0
            fsun_z[p, iv] = 0.0
        end
    end
end

# ---- Default vcmax scaling (coszen <= 0 case) (per patch) ----
@kernel function _sa_vcmaxdef_kernel!(vcmaxcintsun, vcmaxcintsha, @Const(mask),
        @Const(elai), nlevcan::Int, extkn)
    p = @index(Global)
    T = eltype(vcmaxcintsun)
    @inbounds if mask[p]
        extkn_T = oftype(elai[p], extkn)
        if nlevcan == 1
            vcmaxcintsun[p] = 0.0
            vcmaxcintsha[p] = (one(T) - exp(-extkn_T * elai[p])) / extkn_T
            if elai[p] > 0.0
                vcmaxcintsha[p] /= elai[p]
            else
                vcmaxcintsha[p] = 0.0
            end
        elseif nlevcan > 1
            vcmaxcintsun[p] = 0.0
            vcmaxcintsha[p] = 0.0
        end
    end
end

# ---- Non-vegetated patches where coszen > 0 (per (ib, patch)) ----
Base.@kwdef struct _SaNovegOut{M}
    fabd::M; fabd_sun::M; fabd_sha::M
    fabi::M; fabi_sun::M; fabi_sha::M
    ftdd::M; ftid::M; ftii::M
    albd::M; albi::M; albdSF::M; albiSF::M
end
Adapt.@adapt_structure _SaNovegOut

@kernel function _sa_noveg_kernel!(o::_SaNovegOut, @Const(mask), @Const(column),
        @Const(albgrd), @Const(albgri), @Const(albsod), @Const(albsoi), use_SSRE::Bool)
    p, ib = @index(Global, NTuple)
    @inbounds if mask[p]
        c = column[p]
        o.fabd[p, ib] = 0.0
        o.fabd_sun[p, ib] = 0.0
        o.fabd_sha[p, ib] = 0.0
        o.fabi[p, ib] = 0.0
        o.fabi_sun[p, ib] = 0.0
        o.fabi_sha[p, ib] = 0.0
        o.ftdd[p, ib] = 1.0
        o.ftid[p, ib] = 0.0
        o.ftii[p, ib] = 1.0
        o.albd[p, ib] = albgrd[c, ib]
        o.albi[p, ib] = albgri[c, ib]
        if use_SSRE
            o.albdSF[p, ib] = albsod[c, ib]
            o.albiSF[p, ib] = albsoi[c, ib]
        end
    end
end

# ---- History output: ground albedos (per (ib, col)) ----
Base.@kwdef struct _SaHistColOut{M}
    albgrd_hst::M; albgri_hst::M
    albgrd_pur_hst::M; albgri_pur_hst::M
    albgrd_bc_hst::M; albgri_bc_hst::M
    albgrd_oc_hst::M; albgri_oc_hst::M
    albgrd_dst_hst::M; albgri_dst_hst::M
    albsnd_hst2::M; albsni_hst2::M
end
Adapt.@adapt_structure _SaHistColOut
Base.@kwdef struct _SaHistColIn{M}
    albgrd::M; albgri::M
    albgrd_pur::M; albgri_pur::M; albgrd_bc::M; albgri_bc::M
    albgrd_oc::M; albgri_oc::M; albgrd_dst::M; albgri_dst::M
    albsnd_hst::M; albsni_hst::M
end
Adapt.@adapt_structure _SaHistColIn

@kernel function _sa_histcol_kernel!(o::_SaHistColOut, in::_SaHistColIn, @Const(mask),
        @Const(coszen_col), @Const(h2osno_total))
    c, ib = @index(Global, NTuple)
    @inbounds if mask[c]
        if coszen_col[c] > 0.0
            o.albgrd_hst[c, ib] = in.albgrd[c, ib]
            o.albgri_hst[c, ib] = in.albgri[c, ib]
            o.albgrd_pur_hst[c, ib] = in.albgrd_pur[c, ib]
            o.albgri_pur_hst[c, ib] = in.albgri_pur[c, ib]
            o.albgrd_bc_hst[c, ib] = in.albgrd_bc[c, ib]
            o.albgri_bc_hst[c, ib] = in.albgri_bc[c, ib]
            o.albgrd_oc_hst[c, ib] = in.albgrd_oc[c, ib]
            o.albgri_oc_hst[c, ib] = in.albgri_oc[c, ib]
            o.albgrd_dst_hst[c, ib] = in.albgrd_dst[c, ib]
            o.albgri_dst_hst[c, ib] = in.albgri_dst[c, ib]
            if h2osno_total[c] > 0.0
                o.albsnd_hst2[c, ib] = in.albsnd_hst[c, ib]
                o.albsni_hst2[c, ib] = in.albsni_hst[c, ib]
            end
        end
    end
end

# ---- History output: patch albedos (per (ib, patch)) ----
@kernel function _sa_histpatch_kernel!(albd_hst, albi_hst, @Const(mask),
        @Const(coszen_patch), @Const(albd), @Const(albi))
    p, ib = @index(Global, NTuple)
    @inbounds if mask[p]
        if coszen_patch[p] > 0.0
            albd_hst[p, ib] = albd[p, ib]
            albi_hst[p, ib] = albi[p, ib]
        end
    end
end

# --------------------------------------------------------------------------
# surface_albedo_init_time_const!
# --------------------------------------------------------------------------

"""
    surface_albedo_init_time_const!(con, mxsoil_color, soic2d,
                                    col_gridcell, bounds_col, bounds_grc)

Initialize time-constant albedo data: soil color indices, saturated/dry soil
albedos, and melting lake albedos.

Ported from `SurfaceAlbedoInitTimeConst` in `SurfaceAlbedoMod.F90`.

# Arguments
- `con::SurfaceAlbedoConstants` : module constants struct to populate
- `mxsoil_color::Int`          : maximum number of soil color classes (8 or 20)
- `soic2d::Vector{Int}`        : soil color class per gridcell
- `col_gridcell::Vector{Int}`  : gridcell index for each column
- `bounds_col::UnitRange{Int}` : column bounds
- `bounds_grc::UnitRange{Int}` : gridcell bounds (unused but kept for API symmetry)
"""
function surface_albedo_init_time_const!(con::SurfaceAlbedoConstants,
                                          mxsoil_color::Int,
                                          soic2d::Vector{Int},
                                          col_gridcell::Vector{Int},
                                          bounds_col::UnitRange{Int},
                                          bounds_grc::UnitRange{Int})
    # Map gridcell soil color to column
    con.isoicol = zeros(Int, length(bounds_col))
    for c in bounds_col
        g = col_gridcell[c]
        con.isoicol[c] = soic2d[g]
    end

    # Allocate and fill saturated/dry soil albedos
    con.albsat = zeros(mxsoil_color, NUMRAD)
    con.albdry = zeros(mxsoil_color, NUMRAD)

    if mxsoil_color == 8
        con.albsat[1:8, 1] = [0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05]
        con.albsat[1:8, 2] = [0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10]
        con.albdry[1:8, 1] = [0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10]
        con.albdry[1:8, 2] = [0.48, 0.44, 0.40, 0.36, 0.32, 0.28, 0.24, 0.20]
    elseif mxsoil_color == 20
        con.albsat[1:20, 1] = [0.25, 0.23, 0.21, 0.20, 0.19, 0.18, 0.17, 0.16,
                               0.15, 0.14, 0.13, 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04]
        con.albsat[1:20, 2] = [0.50, 0.46, 0.42, 0.40, 0.38, 0.36, 0.34, 0.32,
                               0.30, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.08]
        con.albdry[1:20, 1] = [0.36, 0.34, 0.32, 0.31, 0.30, 0.29, 0.28, 0.27,
                               0.26, 0.25, 0.24, 0.23, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.08]
        con.albdry[1:20, 2] = [0.61, 0.57, 0.53, 0.51, 0.49, 0.48, 0.45, 0.43,
                               0.41, 0.39, 0.37, 0.35, 0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.16]
    else
        error("maximum color class = $mxsoil_color is not supported")
    end

    # Set alblakwi from namelist value
    con.alblakwi .= con.lake_melt_icealb

    return nothing
end

# --------------------------------------------------------------------------
# soil_albedo! — device kernel + bundles
# --------------------------------------------------------------------------

# Outputs / state inputs (per-column / per-(col,band) arrays).
Base.@kwdef struct _SaSoilOut{M}
    albsod::M; albsoi::M
end
Adapt.@adapt_structure _SaSoilOut

Base.@kwdef struct _SaSoilIn{V,M,VI}
    coszen_col::V; landunit::VI; itype::VI; isoicol::VI
    h2osoi_vol::M; t_grnd::V; snl::VI; lake_icefrac::M
end
Adapt.@adapt_structure _SaSoilIn

# Color/band lookup tables copied onto the working backend at working eltype.
Base.@kwdef struct _SaSoilCon{V,M}
    albsat::M; albdry::M; alblakwi::V; albice::V; ablak::V
end
Adapt.@adapt_structure _SaSoilCon

@kernel function _sa_soil_kernel!(o::_SaSoilOut, in::_SaSoilIn, cn::_SaSoilCon,
        @Const(mask), nband::Int, lakepuddling::Bool, calb, tfrz)
    c = @index(Global)
    T = eltype(o.albsod)
    @inbounds if mask[c]
        if in.coszen_col[c] > zero(T)
            l = in.landunit[c]
            it = in.itype[l]
            for ib in 1:nband
                if it == ISTSOIL || it == ISTCROP
                    # Soil
                    inc = smooth_max(T(0.11) - T(0.40) * in.h2osoi_vol[c, 1], zero(T))
                    soilcol = in.isoicol[c]
                    o.albsod[c, ib] = min(cn.albsat[soilcol, ib] + inc,
                                          cn.albdry[soilcol, ib])
                    o.albsoi[c, ib] = o.albsod[c, ib]

                elseif it == ISTICE
                    # Land ice
                    o.albsod[c, ib] = cn.albice[ib]
                    o.albsoi[c, ib] = o.albsod[c, ib]

                elseif in.t_grnd[c] > tfrz ||
                       (lakepuddling && it == ISTDLAK &&
                        in.t_grnd[c] == tfrz &&
                        in.lake_icefrac[c, 1] < one(T) &&
                        in.lake_icefrac[c, 2] > zero(T))
                    # Unfrozen lake/wetland
                    albsod_unfrozen = T(0.05) / (smooth_max(T(0.001), in.coszen_col[c]) + T(0.15))
                    if it == ISTDLAK
                        albsoi_unfrozen = T(0.10)
                    else
                        albsoi_unfrozen = albsod_unfrozen
                    end

                    # Blend with frozen albedo via smooth_heaviside for AD
                    if it == ISTDLAK && !lakepuddling && in.snl[c] == 0
                        sicefr = one(T) - exp(-calb * (tfrz - in.t_grnd[c]) / tfrz)
                        albsod_frozen = sicefr * cn.ablak[ib] +
                            (one(T) - sicefr) * smooth_max(cn.alblakwi[ib],
                                                  T(0.05) / (smooth_max(T(0.001), in.coszen_col[c]) + T(0.15)))
                        albsoi_frozen = sicefr * cn.ablak[ib] +
                            (one(T) - sicefr) * smooth_max(cn.alblakwi[ib], T(0.10))
                    else
                        albsod_frozen = cn.ablak[ib]
                        albsoi_frozen = albsod_frozen
                    end

                    w_unfrozen = smooth_heaviside(in.t_grnd[c] - tfrz)
                    o.albsod[c, ib] = w_unfrozen * albsod_unfrozen + (one(T) - w_unfrozen) * albsod_frozen
                    o.albsoi[c, ib] = w_unfrozen * albsoi_unfrozen + (one(T) - w_unfrozen) * albsoi_frozen

                else
                    # Frozen lake/wetland
                    if it == ISTDLAK && !lakepuddling && in.snl[c] == 0
                        sicefr = one(T) - exp(-calb * (tfrz - in.t_grnd[c]) / tfrz)
                        o.albsod[c, ib] = sicefr * cn.ablak[ib] +
                            (one(T) - sicefr) * smooth_max(cn.alblakwi[ib],
                                                  T(0.05) / (smooth_max(T(0.001), in.coszen_col[c]) + T(0.15)))
                        o.albsoi[c, ib] = sicefr * cn.ablak[ib] +
                            (one(T) - sicefr) * smooth_max(cn.alblakwi[ib], T(0.10))
                    else
                        o.albsod[c, ib] = cn.ablak[ib]
                        o.albsoi[c, ib] = o.albsod[c, ib]
                    end
                end
            end
        end
    end
end

# --------------------------------------------------------------------------
# soil_albedo!
# --------------------------------------------------------------------------

"""
    soil_albedo!(surfalb, con, col, lun,
                 coszen_col, temperature, waterstatebulk, lakestate,
                 mask_nourbanc, bounds_col;
                 lakepuddling=false)

Determine ground surface albedo, accounting for snow/lake/glacier/wetland.

Ported from subroutine `SoilAlbedo` in `SurfaceAlbedoMod.F90`.
"""
function soil_albedo!(surfalb::SurfaceAlbedoData,
                      con::SurfaceAlbedoConstants,
                      col::ColumnData,
                      lun::LandunitData,
                      coszen_col::AbstractVector{<:Real},
                      temperature::TemperatureData,
                      waterstatebulk::WaterStateBulkData,
                      lakestate::LakeStateData,
                      mask_nourbanc::AbstractVector{Bool},
                      bounds_col::UnitRange{Int};
                      lakepuddling::Bool = false)

    nband = NUMRAD

    # Copy host color/band lookup tables + module-global albedo vectors onto the
    # working backend at the working eltype (no-op identity on CPU Arrays).
    T = eltype(surfalb.albsod_col)
    _dev(a) = surfalb.albsod_col isa Array ?
        a : copyto!(similar(surfalb.albsod_col, T, size(a)), T.(a))
    # Int color index: copy onto the backend keeping Int eltype (model off snl).
    isoicol_dev = col.snl isa Array ?
        con.isoicol : copyto!(similar(col.snl, Int, length(con.isoicol)), con.isoicol)

    o = _SaSoilOut(; albsod = surfalb.albsod_col, albsoi = surfalb.albsoi_col)
    in = _SaSoilIn(; coszen_col = coszen_col, landunit = col.landunit,
                   itype = lun.itype, isoicol = isoicol_dev,
                   h2osoi_vol = waterstatebulk.ws.h2osoi_vol_col,
                   t_grnd = temperature.t_grnd_col, snl = col.snl,
                   lake_icefrac = lakestate.lake_icefrac_col)
    cn = _SaSoilCon(; albsat = _dev(con.albsat), albdry = _dev(con.albdry),
                    alblakwi = _dev(con.alblakwi), albice = _dev(ALBICE),
                    ablak = _dev(ALBLAK))

    be = _kernel_backend(surfalb.albsod_col)
    nc = length(bounds_col)
    if nc != 0
        _sa_soil_kernel!(be)(o, in, cn, mask_nourbanc, nband, lakepuddling,
                             T(CALB), T(TFRZ); ndrange = nc)
        KA.synchronize(be)
    end

    return nothing
end

# --------------------------------------------------------------------------
# two_stream! — device kernels + bundles
#
# Fully per-patch independent: each patch computes its own two-stream solution
# from its own inputs and writes its own albedo/flux outputs (no cross-patch
# coupling). Two kernels, both one-thread-per-patch:
#   _ts_params_kernel!   — waveband-independent params (writes per-patch scratch)
#   _ts_waveband_kernel! — internal `for ib in 1:numrad` loop + ib==1 per-layer
# --------------------------------------------------------------------------

# Per-patch waveband-independent scratch (written by params, read by waveband).
Base.@kwdef struct _TsScratch{V,M}
    chil::V; gdir::V; twostext::V; avmu::V; temp0::V; temp2::V; omega::M
end
Adapt.@adapt_structure _TsScratch

@kernel function _ts_params_kernel!(s::_TsScratch, @Const(mask),
        @Const(itype), @Const(coszen_patch), @Const(pftcon_xl))
    p = @index(Global)
    T = eltype(s.chil)
    @inbounds if mask[p]
        cosz = smooth_max(T(0.001), coszen_patch[p])

        s.chil[p] = smooth_clamp(pftcon_xl[itype[p] + 1], T(-0.4), T(0.6))
        if abs(s.chil[p]) <= T(0.01)
            s.chil[p] = T(0.01)
        end
        phi1 = T(0.5) - T(0.633) * s.chil[p] - T(0.330) * s.chil[p] * s.chil[p]
        phi2 = T(0.877) * (T(1.0) - T(2.0) * phi1)
        s.gdir[p] = phi1 + phi2 * cosz
        s.twostext[p] = s.gdir[p] / cosz
        s.avmu[p] = (T(1.0) - phi1 / phi2 * log((phi1 + phi2) / phi1)) / phi2
        s.temp0[p] = smooth_max(s.gdir[p] + phi2 * cosz, T(1.0e-6))
        temp1_val = phi1 * cosz
        s.temp2[p] = (T(1.0) - temp1_val / s.temp0[p] * log((temp1_val + s.temp0[p]) / temp1_val))
    end
end

# Outputs (per-patch / per-(patch,band) flux + albedo arrays).
Base.@kwdef struct _TsOut{V,M}
    albd::M; ftid::M; ftdd::M; fabd::M; albdSF::M
    fabd_sun::M; fabd_sha::M
    albi::M; ftii::M; fabi::M; albiSF::M
    fabi_sun::M; fabi_sha::M
    fsun_z::M; fabd_sun_z::M; fabd_sha_z::M; fabi_sun_z::M; fabi_sha_z::M
    vcmaxcintsun::V; vcmaxcintsha::V
end
Adapt.@adapt_structure _TsOut

# State inputs read by the waveband kernel.
Base.@kwdef struct _TsIn{V,M,VI,Vo}
    column::VI; nrad::VI
    rho_in::M; tau_in::M
    fcansno::V; fwet::V; t_veg::V
    albgrd::M; albgri::M; albsod::M; albsoi::M
    elai::V; esai::V; tlai_z::M; tsai_z::M
    omegas::Vo   # SNICAR snow single-scatter constant — stays Float64 even under AD
                 # (Dual state fields), so it needs its OWN type param to unify.
end
Adapt.@adapt_structure _TsIn

@kernel function _ts_waveband_kernel!(o::_TsOut, in::_TsIn, s::_TsScratch,
        @Const(mask), @Const(coszen_patch), numrad::Int, nlevcan::Int,
        lSFonly::Bool, snowveg_affects_radiation::Bool, betads, betais, tfrz)
    p = @index(Global)
    T = eltype(o.albd)
    @inbounds if mask[p]
        betads_T = oftype(o.albd[p, 1], betads)
        betais_T = oftype(o.albd[p, 1], betais)
        tfrz_T = oftype(o.albd[p, 1], tfrz)
        c = in.column[p]
        cosz = smooth_max(T(0.001), coszen_patch[p])

        for ib in 1:numrad
            # Two-stream parameters omega, betad, betai
            omegal = in.rho_in[p, ib] + in.tau_in[p, ib]
            asu = T(0.5) * omegal * s.gdir[p] / s.temp0[p] * s.temp2[p]
            betadl = (T(1.0) + s.avmu[p] * s.twostext[p]) / (omegal * s.avmu[p] * s.twostext[p]) * asu
            betail = T(0.5) * ((in.rho_in[p, ib] + in.tau_in[p, ib]) + (in.rho_in[p, ib] - in.tau_in[p, ib]) *
                     ((T(1.0) + s.chil[p]) / T(2.0))^2) / omegal

            # Adjust omega, betad, betai for intercepted snow
            if lSFonly || (!snowveg_affects_radiation && in.t_veg[p] > tfrz_T)
                tmp0 = omegal
                tmp1 = betadl
                tmp2 = betail
            else
                if snowveg_affects_radiation
                    tmp0 = (T(1.0) - in.fcansno[p]) * omegal + in.fcansno[p] * in.omegas[ib]
                    tmp1 = ((T(1.0) - in.fcansno[p]) * omegal * betadl + in.fcansno[p] * in.omegas[ib] * betads_T) / tmp0
                    tmp2 = ((T(1.0) - in.fcansno[p]) * omegal * betail + in.fcansno[p] * in.omegas[ib] * betais_T) / tmp0
                else
                    tmp0 = (T(1.0) - in.fwet[p]) * omegal + in.fwet[p] * in.omegas[ib]
                    tmp1 = ((T(1.0) - in.fwet[p]) * omegal * betadl + in.fwet[p] * in.omegas[ib] * betads_T) / tmp0
                    tmp2 = ((T(1.0) - in.fwet[p]) * omegal * betail + in.fwet[p] * in.omegas[ib] * betais_T) / tmp0
                end
            end

            s.omega[p, ib] = tmp0
            betad = tmp1
            betai = tmp2

            # Common terms
            b = T(1.0) - s.omega[p, ib] + s.omega[p, ib] * betai
            c1 = s.omega[p, ib] * betai
            tmp0_val = s.avmu[p] * s.twostext[p]
            d = tmp0_val * s.omega[p, ib] * betad
            f = tmp0_val * s.omega[p, ib] * (T(1.0) - betad)
            tmp1_val = b * b - c1 * c1
            h = sqrt(tmp1_val) / s.avmu[p]
            sigma = tmp0_val * tmp0_val - tmp1_val

            p1 = b + s.avmu[p] * h
            p2 = b - s.avmu[p] * h
            p3 = b + tmp0_val
            p4 = b - tmp0_val

            # Absorbed, reflected, transmitted fluxes for full canopy
            t1 = smooth_min(h * (in.elai[p] + in.esai[p]), T(40.0))
            s1 = exp(-t1)
            t1 = smooth_min(s.twostext[p] * (in.elai[p] + in.esai[p]), T(40.0))
            s2 = exp(-t1)

            # ---- Direct beam ----
            if !lSFonly
                u1 = b - c1 / in.albgrd[c, ib]
                u2 = b - c1 * in.albgrd[c, ib]
                u3 = f + c1 * in.albgrd[c, ib]
            else
                u1 = b - c1 / in.albsod[c, ib]
                u2 = b - c1 * in.albsod[c, ib]
                u3 = f + c1 * in.albsod[c, ib]
            end
            tmp2_val = u1 - s.avmu[p] * h
            tmp3 = u1 + s.avmu[p] * h
            d1 = p1 * tmp2_val / s1 - p2 * tmp3 * s1
            tmp4 = u2 + s.avmu[p] * h
            tmp5 = u2 - s.avmu[p] * h
            d2 = tmp4 / s1 - tmp5 * s1
            h1 = -d * p4 - c1 * f
            tmp6 = d - h1 * p3 / sigma
            tmp7 = (d - c1 - h1 / sigma * (u1 + tmp0_val)) * s2
            h2 = (tmp6 * tmp2_val / s1 - p2 * tmp7) / d1
            h3 = -(tmp6 * tmp3 * s1 - p1 * tmp7) / d1
            h4 = -f * p3 - c1 * d
            tmp8 = h4 / sigma
            tmp9 = (u3 - tmp8 * (u2 - tmp0_val)) * s2
            h5 = -(tmp8 * tmp4 / s1 + tmp9) / d2
            h6 = (tmp8 * tmp5 * s1 + tmp9) / d2

            if !lSFonly
                o.albd[p, ib] = h1 / sigma + h2 + h3
                o.ftid[p, ib] = h4 * s2 / sigma + h5 * s1 + h6 / s1
                o.ftdd[p, ib] = s2
                o.fabd[p, ib] = T(1.0) - o.albd[p, ib] -
                    (T(1.0) - in.albgrd[c, ib]) * o.ftdd[p, ib] -
                    (T(1.0) - in.albgri[c, ib]) * o.ftid[p, ib]
            else
                o.albdSF[p, ib] = h1 / sigma + h2 + h3
            end

            a1 = h1 / sigma * (T(1.0) - s2 * s2) / (T(2.0) * s.twostext[p]) +
                 h2 * (T(1.0) - s2 * s1) / (s.twostext[p] + h) +
                 h3 * (T(1.0) - s2 / s1) / (s.twostext[p] - h)

            a2 = h4 / sigma * (T(1.0) - s2 * s2) / (T(2.0) * s.twostext[p]) +
                 h5 * (T(1.0) - s2 * s1) / (s.twostext[p] + h) +
                 h6 * (T(1.0) - s2 / s1) / (s.twostext[p] - h)

            if !lSFonly
                o.fabd_sun[p, ib] = (T(1.0) - s.omega[p, ib]) *
                    (T(1.0) - s2 + T(1.0) / s.avmu[p] * (a1 + a2))
                o.fabd_sha[p, ib] = o.fabd[p, ib] - o.fabd_sun[p, ib]
            end

            # ---- Diffuse ----
            if !lSFonly
                u1 = b - c1 / in.albgri[c, ib]
                u2 = b - c1 * in.albgri[c, ib]
            else
                u1 = b - c1 / in.albsoi[c, ib]
                u2 = b - c1 * in.albsoi[c, ib]
            end
            tmp2_val = u1 - s.avmu[p] * h
            tmp3 = u1 + s.avmu[p] * h
            d1 = p1 * tmp2_val / s1 - p2 * tmp3 * s1
            tmp4 = u2 + s.avmu[p] * h
            tmp5 = u2 - s.avmu[p] * h
            d2 = tmp4 / s1 - tmp5 * s1
            h7 = (c1 * tmp2_val) / (d1 * s1)
            h8 = (-c1 * tmp3 * s1) / d1
            h9 = tmp4 / (d2 * s1)
            h10 = (-tmp5 * s1) / d2

            if lSFonly
                o.albiSF[p, ib] = h7 + h8
            else
                o.albi[p, ib] = h7 + h8
                o.ftii[p, ib] = h9 * s1 + h10 / s1
                o.fabi[p, ib] = T(1.0) - o.albi[p, ib] -
                    (T(1.0) - in.albgri[c, ib]) * o.ftii[p, ib]

                a1 = h7 * (T(1.0) - s2 * s1) / (s.twostext[p] + h) +
                     h8 * (T(1.0) - s2 / s1) / (s.twostext[p] - h)
                a2 = h9 * (T(1.0) - s2 * s1) / (s.twostext[p] + h) +
                     h10 * (T(1.0) - s2 / s1) / (s.twostext[p] - h)

                o.fabi_sun[p, ib] = (T(1.0) - s.omega[p, ib]) / s.avmu[p] * (a1 + a2)
                o.fabi_sha[p, ib] = o.fabi[p, ib] - o.fabi_sun[p, ib]

                # Per-layer derivatives for PAR (ib == 1)
                if ib == 1
                    if nlevcan == 1
                        # Sun/shade big leaf: single layer
                        o.fsun_z[p, 1] = (T(1.0) - s2) / t1

                        laisum = in.elai[p] + in.esai[p]
                        o.fabd_sun_z[p, 1] = o.fabd_sun[p, ib] /
                            (o.fsun_z[p, 1] * laisum)
                        o.fabi_sun_z[p, 1] = o.fabi_sun[p, ib] /
                            (o.fsun_z[p, 1] * laisum)
                        o.fabd_sha_z[p, 1] = o.fabd_sha[p, ib] /
                            ((T(1.0) - o.fsun_z[p, 1]) * laisum)
                        o.fabi_sha_z[p, 1] = o.fabi_sha[p, ib] /
                            ((T(1.0) - o.fsun_z[p, 1]) * laisum)

                        # Leaf to canopy scaling coefficients
                        extkn = T(0.30)
                        extkb = s.twostext[p]
                        o.vcmaxcintsun[p] = (T(1.0) - exp(-(extkn + extkb) * in.elai[p])) /
                            (extkn + extkb)
                        o.vcmaxcintsha[p] = (T(1.0) - exp(-extkn * in.elai[p])) / extkn -
                            o.vcmaxcintsun[p]
                        if in.elai[p] > T(0.0)
                            o.vcmaxcintsun[p] /= (o.fsun_z[p, 1] * in.elai[p])
                            o.vcmaxcintsha[p] /= ((T(1.0) - o.fsun_z[p, 1]) * in.elai[p])
                        else
                            o.vcmaxcintsun[p] = 0.0
                            o.vcmaxcintsha[p] = 0.0
                        end

                    elseif nlevcan > 1
                        # Multi-layer canopy
                        laisum_l = zero(T)
                        for iv in 1:in.nrad[p]
                            # Cumulative lai+sai at center of layer
                            if iv == 1
                                laisum_l = T(0.5) * (in.tlai_z[p, iv] + in.tsai_z[p, iv])
                            else
                                laisum_l = laisum_l + T(0.5) * ((in.tlai_z[p, iv-1] + in.tsai_z[p, iv-1]) +
                                    (in.tlai_z[p, iv] + in.tsai_z[p, iv]))
                            end

                            t1_l = smooth_min(h * laisum_l, T(40.0))
                            s1_l = exp(-t1_l)
                            t1_l = smooth_min(s.twostext[p] * laisum_l, T(40.0))
                            s2_l = exp(-t1_l)
                            o.fsun_z[p, iv] = s2_l

                            # ---- Direct beam derivatives ----
                            u1_l = b - c1 / in.albgrd[c, ib]
                            u2_l = b - c1 * in.albgrd[c, ib]
                            u3_l = f + c1 * in.albgrd[c, ib]

                            tmp2_l = u1_l - s.avmu[p] * h
                            tmp3_l = u1_l + s.avmu[p] * h
                            d1_l = p1 * tmp2_l / s1_l - p2 * tmp3_l * s1_l
                            tmp4_l = u2_l + s.avmu[p] * h
                            tmp5_l = u2_l - s.avmu[p] * h
                            d2_l = tmp4_l / s1_l - tmp5_l * s1_l

                            h1_l = -d * p4 - c1 * f
                            tmp6_l = d - h1_l * p3 / sigma
                            tmp7_l = (d - c1 - h1_l / sigma * (u1_l + tmp0_val)) * s2_l

                            h2_l = (tmp6_l * tmp2_l / s1_l - p2 * tmp7_l) / d1_l
                            h3_l = -(tmp6_l * tmp3_l * s1_l - p1 * tmp7_l) / d1_l

                            h4_l = -f * p3 - c1 * d
                            tmp8_l = h4_l / sigma
                            tmp9_l = (u3_l - tmp8_l * (u2_l - tmp0_val)) * s2_l
                            h5_l = -(tmp8_l * tmp4_l / s1_l + tmp9_l) / d2_l
                            h6_l = (tmp8_l * tmp5_l * s1_l + tmp9_l) / d2_l

                            # Derivatives
                            v_val = d1_l
                            dv = h * p1 * tmp2_l / s1_l + h * p2 * tmp3_l * s1_l

                            u_val = tmp6_l * tmp2_l / s1_l - p2 * tmp7_l
                            du = h * tmp6_l * tmp2_l / s1_l + s.twostext[p] * p2 * tmp7_l
                            dh2 = (v_val * du - u_val * dv) / (v_val * v_val)

                            u_val = -tmp6_l * tmp3_l * s1_l + p1 * tmp7_l
                            du = h * tmp6_l * tmp3_l * s1_l - s.twostext[p] * p1 * tmp7_l
                            dh3 = (v_val * du - u_val * dv) / (v_val * v_val)

                            v_val = d2_l
                            dv = h * tmp4_l / s1_l + h * tmp5_l * s1_l

                            u_val = -h4_l / sigma * tmp4_l / s1_l - tmp9_l
                            du = -h * h4_l / sigma * tmp4_l / s1_l + s.twostext[p] * tmp9_l
                            dh5 = (v_val * du - u_val * dv) / (v_val * v_val)

                            u_val = h4_l / sigma * tmp5_l * s1_l + tmp9_l
                            du = -h * h4_l / sigma * tmp5_l * s1_l - s.twostext[p] * tmp9_l
                            dh6 = (v_val * du - u_val * dv) / (v_val * v_val)

                            da1 = h1_l / sigma * s2_l * s2_l + h2_l * s2_l * s1_l + h3_l * s2_l / s1_l +
                                  (T(1.0) - s2_l * s1_l) / (s.twostext[p] + h) * dh2 +
                                  (T(1.0) - s2_l / s1_l) / (s.twostext[p] - h) * dh3
                            da2 = h4_l / sigma * s2_l * s2_l + h5_l * s2_l * s1_l + h6_l * s2_l / s1_l +
                                  (T(1.0) - s2_l * s1_l) / (s.twostext[p] + h) * dh5 +
                                  (T(1.0) - s2_l / s1_l) / (s.twostext[p] - h) * dh6

                            d_ftid = -s.twostext[p] * h4_l / sigma * s2_l - h * h5_l * s1_l + h * h6_l / s1_l +
                                     dh5 * s1_l + dh6 / s1_l
                            d_fabd = -(dh2 + dh3) + (T(1.0) - in.albgrd[c, ib]) * s.twostext[p] * s2_l -
                                     (T(1.0) - in.albgri[c, ib]) * d_ftid
                            d_fabd_sun = (T(1.0) - s.omega[p, ib]) *
                                (s.twostext[p] * s2_l + T(1.0) / s.avmu[p] * (da1 + da2))
                            d_fabd_sha = d_fabd - d_fabd_sun

                            o.fabd_sun_z[p, iv] = smooth_max(d_fabd_sun, T(0.0)) /
                                o.fsun_z[p, iv]
                            o.fabd_sha_z[p, iv] = smooth_max(d_fabd_sha, T(0.0)) /
                                (T(1.0) - o.fsun_z[p, iv])

                            # ---- Diffuse derivatives ----
                            u1_l = b - c1 / in.albgri[c, ib]
                            u2_l = b - c1 * in.albgri[c, ib]

                            tmp2_l = u1_l - s.avmu[p] * h
                            tmp3_l = u1_l + s.avmu[p] * h
                            d1_l = p1 * tmp2_l / s1_l - p2 * tmp3_l * s1_l
                            tmp4_l = u2_l + s.avmu[p] * h
                            tmp5_l = u2_l - s.avmu[p] * h
                            d2_l = tmp4_l / s1_l - tmp5_l * s1_l

                            h7_l = (c1 * tmp2_l) / (d1_l * s1_l)
                            h8_l = (-c1 * tmp3_l * s1_l) / d1_l
                            h9_l = tmp4_l / (d2_l * s1_l)
                            h10_l = (-tmp5_l * s1_l) / d2_l

                            a1_l = h7_l * (T(1.0) - s2_l * s1_l) / (s.twostext[p] + h) +
                                   h8_l * (T(1.0) - s2_l / s1_l) / (s.twostext[p] - h)
                            a2_l = h9_l * (T(1.0) - s2_l * s1_l) / (s.twostext[p] + h) +
                                   h10_l * (T(1.0) - s2_l / s1_l) / (s.twostext[p] - h)

                            v_val = d1_l
                            dv = h * p1 * tmp2_l / s1_l + h * p2 * tmp3_l * s1_l

                            u_val = c1 * tmp2_l / s1_l
                            du = h * c1 * tmp2_l / s1_l
                            dh7 = (v_val * du - u_val * dv) / (v_val * v_val)

                            u_val = -c1 * tmp3_l * s1_l
                            du = h * c1 * tmp3_l * s1_l
                            dh8 = (v_val * du - u_val * dv) / (v_val * v_val)

                            v_val = d2_l
                            dv = h * tmp4_l / s1_l + h * tmp5_l * s1_l

                            u_val = tmp4_l / s1_l
                            du = h * tmp4_l / s1_l
                            dh9 = (v_val * du - u_val * dv) / (v_val * v_val)

                            u_val = -tmp5_l * s1_l
                            du = h * tmp5_l * s1_l
                            dh10 = (v_val * du - u_val * dv) / (v_val * v_val)

                            da1 = h7_l * s2_l * s1_l + h8_l * s2_l / s1_l +
                                  (T(1.0) - s2_l * s1_l) / (s.twostext[p] + h) * dh7 +
                                  (T(1.0) - s2_l / s1_l) / (s.twostext[p] - h) * dh8
                            da2 = h9_l * s2_l * s1_l + h10_l * s2_l / s1_l +
                                  (T(1.0) - s2_l * s1_l) / (s.twostext[p] + h) * dh9 +
                                  (T(1.0) - s2_l / s1_l) / (s.twostext[p] - h) * dh10

                            d_ftii = -h * h9_l * s1_l + h * h10_l / s1_l + dh9 * s1_l + dh10 / s1_l
                            d_fabi = -(dh7 + dh8) - (T(1.0) - in.albgri[c, ib]) * d_ftii
                            d_fabi_sun = (T(1.0) - s.omega[p, ib]) / s.avmu[p] * (da1 + da2)
                            d_fabi_sha = d_fabi - d_fabi_sun

                            o.fabi_sun_z[p, iv] = smooth_max(d_fabi_sun, T(0.0)) /
                                o.fsun_z[p, iv]
                            o.fabi_sha_z[p, iv] = smooth_max(d_fabi_sha, T(0.0)) /
                                (T(1.0) - o.fsun_z[p, iv])
                        end  # iv loop
                    end  # nlevcan
                end  # ib == 1
            end  # !lSFonly
        end  # ib loop
    end
end

# --------------------------------------------------------------------------
# two_stream!
# --------------------------------------------------------------------------

"""
    two_stream!(surfalb, patchdata, col, canopystate, temperature,
                waterdiagbulk, coszen_patch, rho, tau,
                mask_vegsol, bounds_patch;
                SFonly=false)

Two-stream fluxes for canopy radiative transfer.

Uses two-stream approximation of Dickinson (1983) and Sellers (1985) to
calculate fluxes absorbed, reflected, and transmitted by vegetation.
Calculates sunlit and shaded fluxes as described by Bonan et al (2011).

Ported from subroutine `TwoStream` in `SurfaceAlbedoMod.F90`.

# Arguments
- `mask_vegsol::BitVector` : mask for vegetated patches with coszen > 0
- `SFonly::Bool`           : if true, only calculate snow-free albedos
"""
function two_stream!(surfalb::SurfaceAlbedoData,
                     patchdata::PatchData,
                     col::ColumnData,
                     canopystate::CanopyStateData,
                     temperature::TemperatureData,
                     waterdiagbulk::WaterDiagnosticBulkData,
                     coszen_patch::AbstractVector{<:Real},
                     rho_in::AbstractMatrix{<:Real},
                     tau_in::AbstractMatrix{<:Real},
                     pftcon_xl::AbstractVector{<:Real},
                     mask_vegsol::AbstractVector{Bool},
                     bounds_patch::UnitRange{Int};
                     SFonly::Bool = false)

    numrad = NUMRAD
    nlevcan = NLEVCAN
    lSFonly = SFonly

    # Aliases
    elai = canopystate.elai_patch
    esai = canopystate.esai_patch
    t_veg = temperature.t_veg_patch
    fwet = waterdiagbulk.fwet_patch
    fcansno = waterdiagbulk.fcansno_patch
    tlai_z = surfalb.tlai_z_patch
    tsai_z = surfalb.tsai_z_patch
    nrad = surfalb.nrad_patch
    albgrd = surfalb.albgrd_col
    albgri = surfalb.albgri_col
    albsod = surfalb.albsod_col
    albsoi = surfalb.albsoi_col

    np = length(bounds_patch)
    np == 0 && return nothing
    T = eltype(coszen_patch)

    # Device-resident per-patch scratch (waveband-independent params + omega),
    # shared between the two kernels. fill!/similar tracks the working backend.
    s = _TsScratch(;
        chil     = fill!(similar(coszen_patch, T, np), zero(T)),
        gdir     = fill!(similar(coszen_patch, T, np), zero(T)),
        twostext = fill!(similar(coszen_patch, T, np), zero(T)),
        avmu     = fill!(similar(coszen_patch, T, np), zero(T)),
        temp0    = fill!(similar(coszen_patch, T, np), zero(T)),
        temp2    = fill!(similar(coszen_patch, T, np), zero(T)),
        omega    = fill!(similar(coszen_patch, T, np, numrad), zero(T)))

    # Module-global single-scattering albedo lookup onto the working backend at
    # working eltype (no-op identity on CPU Arrays).
    omegas_dev = coszen_patch isa Array ? OMEGAS :
        copyto!(similar(coszen_patch, T, length(OMEGAS)), T.(OMEGAS))

    be = _kernel_backend(coszen_patch)

    # Calculate two-stream parameters independent of waveband
    _ts_params_kernel!(be)(s, mask_vegsol, patchdata.itype, coszen_patch,
                           pftcon_xl; ndrange = np)
    KA.synchronize(be)

    # Loop over wavebands (internal `for ib` in the per-patch kernel)
    o = _TsOut(;
        albd = surfalb.albd_patch, ftid = surfalb.ftid_patch,
        ftdd = surfalb.ftdd_patch, fabd = surfalb.fabd_patch,
        albdSF = surfalb.albdSF_patch,
        fabd_sun = surfalb.fabd_sun_patch, fabd_sha = surfalb.fabd_sha_patch,
        albi = surfalb.albi_patch, ftii = surfalb.ftii_patch,
        fabi = surfalb.fabi_patch, albiSF = surfalb.albiSF_patch,
        fabi_sun = surfalb.fabi_sun_patch, fabi_sha = surfalb.fabi_sha_patch,
        fsun_z = surfalb.fsun_z_patch,
        fabd_sun_z = surfalb.fabd_sun_z_patch, fabd_sha_z = surfalb.fabd_sha_z_patch,
        fabi_sun_z = surfalb.fabi_sun_z_patch, fabi_sha_z = surfalb.fabi_sha_z_patch,
        vcmaxcintsun = surfalb.vcmaxcintsun_patch,
        vcmaxcintsha = surfalb.vcmaxcintsha_patch)
    inb = _TsIn(;
        column = patchdata.column, nrad = nrad,
        rho_in = rho_in, tau_in = tau_in,
        fcansno = fcansno, fwet = fwet, t_veg = t_veg,
        albgrd = albgrd, albgri = albgri, albsod = albsod, albsoi = albsoi,
        elai = elai, esai = esai, tlai_z = tlai_z, tsai_z = tsai_z,
        omegas = omegas_dev)

    _ts_waveband_kernel!(be)(o, inb, s, mask_vegsol, coszen_patch, numrad,
                             nlevcan, lSFonly, surfalb_con.snowveg_affects_radiation,
                             T(BETADS), T(BETAIS), T(TFRZ); ndrange = np)
    KA.synchronize(be)

    return nothing
end

# --------------------------------------------------------------------------
# surface_albedo!  (top-level driver)
# --------------------------------------------------------------------------

"""
    surface_albedo!(surfalb, con, grc, col, lun, patchdata,
                    canopystate, temperature, waterstatebulk,
                    waterdiagbulk, lakestate, aerosol,
                    mask_nourbanc, mask_nourbanp,
                    nextsw_cday, declinp1,
                    bounds_grc, bounds_col, bounds_patch,
                    pftcon_rhol, pftcon_rhos, pftcon_taul, pftcon_taus,
                    pftcon_xl;
                    use_SSRE=false, use_snicar_frc=false, use_fates=false,
                    use_subgrid_fluxes=true, do_sno_oc=false,
                    lakepuddling=false,
                    downscale_hillslope_meteorology=false)

Surface albedo and two-stream fluxes (top-level driver).

Ported from subroutine `SurfaceAlbedo` in `SurfaceAlbedoMod.F90`.
Uses SNICAR_RT for snow albedo when optics tables are loaded.
FATES calls are represented as stubs.
"""
function surface_albedo!(surfalb::SurfaceAlbedoData,
                         con::SurfaceAlbedoConstants,
                         grc::GridcellData,
                         col::ColumnData,
                         lun::LandunitData,
                         patchdata::PatchData,
                         canopystate::CanopyStateData,
                         temperature::TemperatureData,
                         waterstatebulk::WaterStateBulkData,
                         waterdiagbulk::WaterDiagnosticBulkData,
                         lakestate::LakeStateData,
                         aerosol::AerosolData,
                         mask_nourbanc::AbstractVector{Bool},
                         mask_nourbanp::AbstractVector{Bool},
                         nextsw_cday::Real,
                         declinp1::Real,
                         bounds_grc::UnitRange{Int},
                         bounds_col::UnitRange{Int},
                         bounds_patch::UnitRange{Int},
                         pftcon_rhol::AbstractMatrix{<:Real},
                         pftcon_rhos::AbstractMatrix{<:Real},
                         pftcon_taul::AbstractMatrix{<:Real},
                         pftcon_taus::AbstractMatrix{<:Real},
                         pftcon_xl::AbstractVector{<:Real},
                         coszen_func::Function;
                         use_SSRE::Bool = false,
                         use_snicar_frc::Bool = false,
                         use_fates::Bool = false,
                         use_subgrid_fluxes::Bool = true,
                         do_sno_oc::Bool = false,
                         lakepuddling::Bool = false,
                         downscale_hillslope_meteorology::Bool = false)

    numrad = NUMRAD
    nlevcan = NLEVCAN
    nlevsno = varpar.nlevsno
    mpe = MPE_ALBEDO

    # Aliases
    coszen_grc = surfalb.coszen_grc
    coszen_col_arr = surfalb.coszen_col
    albgrd = surfalb.albgrd_col
    albgri = surfalb.albgri_col
    albsod = surfalb.albsod_col
    albsoi = surfalb.albsoi_col
    albd = surfalb.albd_patch
    albi = surfalb.albi_patch
    frac_sno = waterdiagbulk.frac_sno_col
    tlai = canopystate.tlai_patch
    tsai = canopystate.tsai_patch
    elai = canopystate.elai_patch
    esai = canopystate.esai_patch
    nrad = surfalb.nrad_patch
    ncan = surfalb.ncan_patch
    tlai_z = surfalb.tlai_z_patch
    tsai_z = surfalb.tsai_z_patch
    fsun_z = surfalb.fsun_z_patch

    # --- Cosine solar zenith angle ---
    # coszen_func is a host closure (not GPU-safe). Evaluate it on the HOST over
    # Array copies of lat/lon (only ng elements), then copy the result into the
    # (possibly device-resident) coszen_grc. Byte-identical to the scalar loop.
    let lat_h = Array(grc.lat), lon_h = Array(grc.lon),
        coszen_grc_h = Array(coszen_grc)
        for g in bounds_grc
            coszen_grc_h[g] = coszen_func(nextsw_cday, lat_h[g], lon_h[g], declinp1)
        end
        copyto!(coszen_grc, coszen_grc_h)
    end

    FT = eltype(coszen_col_arr)
    coszen_patch_arr = fill!(similar(coszen_col_arr, FT, length(bounds_patch)), zero(FT))

    # Hillslope downscaling branch is currently a no-op (both branches identical),
    # so coszen_col[c] = coszen_grc[gridcell[c]] unconditionally.
    surfalb_coszen_col!(coszen_col_arr, col.gridcell, coszen_grc)

    surfalb_coszen_patch!(coszen_patch_arr, mask_nourbanp, patchdata.column, coszen_col_arr)

    # --- Initialize output ---
    let initcol = _SaInitColOut(;
            albsod = albsod, albsoi = albsoi, albgrd = albgrd, albgri = albgri,
            albgrd_pur = surfalb.albgrd_pur_col, albgri_pur = surfalb.albgri_pur_col,
            albgrd_bc = surfalb.albgrd_bc_col, albgri_bc = surfalb.albgri_bc_col,
            albgrd_oc = surfalb.albgrd_oc_col, albgri_oc = surfalb.albgri_oc_col,
            albgrd_dst = surfalb.albgrd_dst_col, albgri_dst = surfalb.albgri_dst_col,
            albgrd_hst = surfalb.albgrd_hst_col, albgri_hst = surfalb.albgri_hst_col,
            albgrd_pur_hst = surfalb.albgrd_pur_hst_col, albgri_pur_hst = surfalb.albgri_pur_hst_col,
            albgrd_bc_hst = surfalb.albgrd_bc_hst_col, albgri_bc_hst = surfalb.albgri_bc_hst_col,
            albgrd_oc_hst = surfalb.albgrd_oc_hst_col, albgri_oc_hst = surfalb.albgri_oc_hst_col,
            albgrd_dst_hst = surfalb.albgrd_dst_hst_col, albgri_dst_hst = surfalb.albgri_dst_hst_col,
            albsnd_hst2 = surfalb.albsnd_hst2_col, albsni_hst2 = surfalb.albsni_hst2_col,
            flx_absdv = surfalb.flx_absdv_col, flx_absdn = surfalb.flx_absdn_col,
            flx_absiv = surfalb.flx_absiv_col, flx_absin = surfalb.flx_absin_col)
        be = _kernel_backend(albsod)
        nflx = size(surfalb.flx_absdv_col, 2)
        if length(bounds_col) * numrad != 0
            _sa_init_col_kernel!(be)(initcol, mask_nourbanc, nflx;
                                     ndrange = (length(bounds_col), numrad))
            KA.synchronize(be)
        end
    end

    let initpch = _SaInitPatchOut(;
            albd = albd, albi = albi,
            albd_hst = surfalb.albd_hst_patch, albi_hst = surfalb.albi_hst_patch,
            albdSF = surfalb.albdSF_patch, albiSF = surfalb.albiSF_patch,
            fabd = surfalb.fabd_patch, fabd_sun = surfalb.fabd_sun_patch, fabd_sha = surfalb.fabd_sha_patch,
            fabi = surfalb.fabi_patch, fabi_sun = surfalb.fabi_sun_patch, fabi_sha = surfalb.fabi_sha_patch,
            ftdd = surfalb.ftdd_patch, ftid = surfalb.ftid_patch, ftii = surfalb.ftii_patch)
        be = _kernel_backend(albd)
        if length(bounds_patch) * numrad != 0
            _sa_init_patch_kernel!(be)(initpch, mask_nourbanp, use_SSRE;
                                       ndrange = (length(bounds_patch), numrad))
            KA.synchronize(be)
        end
    end

    # --- Soil albedo ---
    soil_albedo!(surfalb, con, col, lun,
                 coszen_col_arr, temperature, waterstatebulk, lakestate,
                 mask_nourbanc, bounds_col; lakepuddling=lakepuddling)

    # --- Snow albedos via SNICAR radiative transfer ---
    nc = length(bounds_col)
    nlevsno_val = nlevsno

    # Check if SNICAR optics tables are loaded
    use_snicar = length(snicar_optics.ext_cff_mss_snw_drc) > 0 &&
                 any(x -> x != 0.0, @view snicar_optics.ext_cff_mss_snw_drc[1:min(10, end), 1:min(1, end)])

    # Total snow water content (kernelized; device-resident so downstream kernels
    # can read it).
    h2osno_total_arr = fill!(similar(coszen_col_arr, FT, nc), zero(FT))
    waterstate_calculate_total_h2osno!(waterstatebulk.ws, mask_nourbanc,
                                        bounds_col, col.snl, h2osno_total_arr)

    # Snow albedo / flux-absorption outputs — device-resident (consumed by the
    # ground-albedo, flux-map and history kernels below).
    albsnd = fill!(similar(coszen_col_arr, FT, nc, numrad), zero(FT))
    albsni = fill!(similar(coszen_col_arr, FT, nc, numrad), zero(FT))
    flx_absd_snw = fill!(similar(coszen_col_arr, FT, nc, nlevsno_val + 1, numrad), zero(FT))
    flx_absi_snw = fill!(similar(coszen_col_arr, FT, nc, nlevsno_val + 1, numrad), zero(FT))

    if use_snicar && albsnd isa Array
        # CPU path: snicar_rt! runs in-place on host arrays. (The on-device port
        # snicar_rt_device! handles the GPU branch below; this host path is kept
        # byte-identical for CPU runs and AD.) Build its inputs on the HOST from
        # Array() copies of the (CPU) state, run it on host scratch outputs, then
        # copy the results back into albsnd/albsni/flx_*_snw.
        h2osno_liq_snw = zeros(FT, nc, nlevsno_val)
        h2osno_ice_snw = zeros(FT, nc, nlevsno_val)
        snw_rds_int = fill(round(Int, snicar_params.snw_rds_min), nc, nlevsno_val)
        mss_cnc_aer = zeros(FT, nc, nlevsno_val, SNO_NBR_AER)
        albsfc_snw = zeros(FT, nc, numrad)

        mask_h = Array(mask_nourbanc)
        h2osoi_liq_h = Array(waterstatebulk.ws.h2osoi_liq_col)
        h2osoi_ice_h = Array(waterstatebulk.ws.h2osoi_ice_col)
        snw_rds_h = Array(waterdiagbulk.snw_rds_col)
        albsoi_h = Array(albsoi)
        aer_avail = size(aerosol.mss_cnc_bcphi_col, 1)
        bcphi_h = Array(aerosol.mss_cnc_bcphi_col); bcpho_h = Array(aerosol.mss_cnc_bcpho_col)
        ocphi_h = Array(aerosol.mss_cnc_ocphi_col); ocpho_h = Array(aerosol.mss_cnc_ocpho_col)
        dst1_h = Array(aerosol.mss_cnc_dst1_col); dst2_h = Array(aerosol.mss_cnc_dst2_col)
        dst3_h = Array(aerosol.mss_cnc_dst3_col); dst4_h = Array(aerosol.mss_cnc_dst4_col)

        for c in bounds_col
            mask_h[c] || continue
            # Snow layer water content (first nlevsno slots are snow layers)
            for j in 1:nlevsno_val
                h2osno_liq_snw[c, j] = smooth_max(0.0, h2osoi_liq_h[c, j])
                h2osno_ice_snw[c, j] = smooth_max(0.0, h2osoi_ice_h[c, j])
            end
            # Snow grain radius (Float64 → Int, microns)
            for j in 1:nlevsno_val
                rds = snw_rds_h[c, j]
                if !isnan(rds) && rds > 0.0
                    snw_rds_int[c, j] = round(Int, clamp(rds, SNW_RDS_MIN_TBL, SNW_RDS_MAX_TBL))
                end
            end
            # Aerosol mass concentrations (8 species: bcphi, bcpho, ocphi, ocpho, dst1-4)
            if aer_avail >= c
                for j in 1:nlevsno_val
                    mss_cnc_aer[c, j, 1] = bcphi_h[c, j]
                    mss_cnc_aer[c, j, 2] = bcpho_h[c, j]
                    mss_cnc_aer[c, j, 3] = ocphi_h[c, j]
                    mss_cnc_aer[c, j, 4] = ocpho_h[c, j]
                    mss_cnc_aer[c, j, 5] = dst1_h[c, j]
                    mss_cnc_aer[c, j, 6] = dst2_h[c, j]
                    mss_cnc_aer[c, j, 7] = dst3_h[c, j]
                    mss_cnc_aer[c, j, 8] = dst4_h[c, j]
                end
            end
            # Underlying surface albedo (diffuse soil albedo)
            for ib in 1:numrad
                albsfc_snw[c, ib] = albsoi_h[c, ib]
            end
        end

        coszen_h = Array(coszen_col_arr)
        h2osno_total_h = Array(h2osno_total_arr)
        snl_h = Array(col.snl)
        frac_sno_h = Array(frac_sno)
        mask_bv = mask_h isa BitVector ? mask_h : BitVector(mask_h)

        albsnd_h = zeros(FT, nc, numrad)
        albsni_h = zeros(FT, nc, numrad)
        flx_absd_h = zeros(FT, nc, nlevsno_val + 1, numrad)
        flx_absi_h = zeros(FT, nc, nlevsno_val + 1, numrad)

        # Call SNICAR_RT for direct beam (flg_slr=1)
        snicar_rt!(coszen_h, 1,
                   h2osno_liq_snw, h2osno_ice_snw, h2osno_total_h,
                   snw_rds_int, mss_cnc_aer, albsfc_snw,
                   snl_h, frac_sno_h,
                   albsnd_h, flx_absd_h, nlevsno_val;
                   mask_nourbanc=mask_bv)

        # Call SNICAR_RT for diffuse (flg_slr=2)
        snicar_rt!(coszen_h, 2,
                   h2osno_liq_snw, h2osno_ice_snw, h2osno_total_h,
                   snw_rds_int, mss_cnc_aer, albsfc_snw,
                   snl_h, frac_sno_h,
                   albsni_h, flx_absi_h, nlevsno_val;
                   mask_nourbanc=mask_bv)

        copyto!(albsnd, albsnd_h)
        copyto!(albsni, albsni_h)
        copyto!(flx_absd_snw, flx_absd_h)
        copyto!(flx_absi_snw, flx_absi_h)
    elseif use_snicar
        # GPU path: build the SNICAR snow-layer inputs on-device (prep kernel),
        # then run the on-device Adding-Doubling solver snicar_rt_device! straight
        # into the device-resident albsnd/albsni/flx_*_snw. No host roundtrip.
        h2osno_liq_snw = fill!(similar(albsnd, FT, nc, nlevsno_val), zero(FT))
        h2osno_ice_snw = fill!(similar(albsnd, FT, nc, nlevsno_val), zero(FT))
        snw_rds_dev    = fill!(similar(col.snl, eltype(col.snl), nc, nlevsno_val), zero(eltype(col.snl)))
        mss_cnc_dev    = fill!(similar(albsnd, FT, nc, nlevsno_val, SNO_NBR_AER), zero(FT))
        aer_avail      = size(aerosol.mss_cnc_bcphi_col, 1)
        rds_min_int    = round(Int, snicar_params.snw_rds_min)

        _launch!(_snicar_prep_kernel!, h2osno_liq_snw, h2osno_ice_snw, snw_rds_dev, mss_cnc_dev,
                 mask_nourbanc, waterstatebulk.ws.h2osoi_liq_col, waterstatebulk.ws.h2osoi_ice_col,
                 waterdiagbulk.snw_rds_col,
                 aerosol.mss_cnc_bcphi_col, aerosol.mss_cnc_bcpho_col,
                 aerosol.mss_cnc_ocphi_col, aerosol.mss_cnc_ocpho_col,
                 aerosol.mss_cnc_dst1_col, aerosol.mss_cnc_dst2_col,
                 aerosol.mss_cnc_dst3_col, aerosol.mss_cnc_dst4_col,
                 nlevsno_val, aer_avail, rds_min_int; ndrange = nc)

        # Underlying surface albedo for SNICAR is the diffuse soil albedo (albsoi).
        snicar_rt_device!(coszen_col_arr, 1, h2osno_liq_snw, h2osno_ice_snw, h2osno_total_arr,
                          snw_rds_dev, mss_cnc_dev, albsoi, col.snl, frac_sno,
                          albsnd, flx_absd_snw, nlevsno_val; mask = mask_nourbanc)
        snicar_rt_device!(coszen_col_arr, 2, h2osno_liq_snw, h2osno_ice_snw, h2osno_total_arr,
                          snw_rds_dev, mss_cnc_dev, albsoi, col.snl, frac_sno,
                          albsni, flx_absi_snw, nlevsno_val; mask = mask_nourbanc)
    else
        # Fallback: age-based snow albedo (used when SNICAR optics not loaded)
        snow_persist = waterstatebulk.snow_persistence_col
        _launch!(_sa_snow_fallback_kernel!, albsnd, albsni, mask_nourbanc,
                 coszen_col_arr, frac_sno, snow_persist, length(snow_persist);
                 ndrange = nc)
    end

    # --- Ground albedos (snow-fraction weighting) ---
    let go = _SaGrndOut(;
            albgrd = albgrd, albgri = albgri,
            albgrd_bc = surfalb.albgrd_bc_col, albgri_bc = surfalb.albgri_bc_col,
            albgrd_dst = surfalb.albgrd_dst_col, albgri_dst = surfalb.albgri_dst_col,
            albgrd_pur = surfalb.albgrd_pur_col, albgri_pur = surfalb.albgri_pur_col,
            albgrd_oc = surfalb.albgrd_oc_col, albgri_oc = surfalb.albgri_oc_col),
        gin = _SaGrndIn(; albsod = albsod, albsoi = albsoi,
                        albsnd = albsnd, albsni = albsni, frac_sno = frac_sno)
        be = _kernel_backend(albgrd)
        if length(bounds_col) * numrad != 0
            _sa_grnd_kernel!(be)(go, gin, mask_nourbanc, coszen_col_arr,
                                 use_snicar_frc, do_sno_oc;
                                 ndrange = (length(bounds_col), numrad))
            KA.synchronize(be)
        end
    end

    # --- Map SNICAR flux absorption to surfalb fields ---
    if use_snicar
        fm = _SaFluxMap(;
            flx_absdv = surfalb.flx_absdv_col, flx_absdn = surfalb.flx_absdn_col,
            flx_absiv = surfalb.flx_absiv_col, flx_absin = surfalb.flx_absin_col,
            flx_absd_snw = flx_absd_snw, flx_absi_snw = flx_absi_snw,
            albsnd = albsnd, albsni = albsni, albsod = albsod, albsoi = albsoi,
            frac_sno = frac_sno, landunit = col.landunit, lakpoi = lun.lakpoi)
        be = _kernel_backend(surfalb.flx_absdv_col)
        if nc != 0
            _sa_fluxmap_kernel!(be)(fm, mask_nourbanc, nlevsno_val + 1,
                                    use_subgrid_fluxes, IVIS, INIR; ndrange = nc)
            KA.synchronize(be)
        end
    end

    # --- Snow albedo history diagnostics ---
    _launch!(_sa_snowhist_kernel!, surfalb.albsnd_hst_col,
             surfalb.albsni_hst_col, mask_nourbanc, coszen_col_arr,
             h2osno_total_arr, albsnd, albsni;
             ndrange = (length(bounds_col), numrad))

    # --- Canopy layering ---
    dincmax = 0.25
    let cl = _SaCanLayer(; nrad = nrad, ncan = ncan, tlai_z = tlai_z, tsai_z = tsai_z,
                         elai = elai, esai = esai, tlai = tlai, tsai = tsai)
        be = _kernel_backend(nrad)
        if length(bounds_patch) != 0
            _sa_canlayer_kernel!(be)(cl, mask_nourbanp, nlevcan,
                                     eltype(tlai_z)(dincmax), eltype(tlai_z)(mpe);
                                     ndrange = length(bounds_patch))
            KA.synchronize(be)
        end
    end

    # --- Zero fluxes for active canopy layers ---
    _launch!(_sa_zerolayer_kernel!, surfalb.fabd_sun_z_patch,
             surfalb.fabd_sha_z_patch, surfalb.fabi_sun_z_patch,
             surfalb.fabi_sha_z_patch, fsun_z, mask_nourbanp, nrad;
             ndrange = length(bounds_patch))

    # --- Default vcmax scaling (coszen <= 0 case) ---
    extkn = 0.30
    _launch!(_sa_vcmaxdef_kernel!, surfalb.vcmaxcintsun_patch,
             surfalb.vcmaxcintsha_patch, mask_nourbanp, elai, nlevcan,
             eltype(surfalb.vcmaxcintsun_patch)(extkn);
             ndrange = length(bounds_patch))

    # --- Create solar-vegetated filter ---
    # mask_vegsol/mask_novegsol are BitVector scratch consumed by two_stream! and
    # the non-veg patch kernel. BitVector has no device backend (packed bits), so
    # this classification stays on the HOST over Array() copies of the device state
    # it reads. Byte-identical (Array() of a CPU Array is a no-op).
    mask_vegsol = falses(length(bounds_patch))
    mask_novegsol = falses(length(bounds_patch))

    let maskp_h = Array(mask_nourbanp), coszp_h = Array(coszen_patch_arr),
        landunit_h = Array(patchdata.landunit), itype_h = Array(lun.itype),
        elai_h = Array(elai), esai_h = Array(esai)
        for p in bounds_patch
            maskp_h[p] || continue
            if coszp_h[p] > 0.0
                l_p = landunit_h[p]
                if (itype_h[l_p] == ISTSOIL || itype_h[l_p] == ISTCROP) &&
                   (elai_h[p] + esai_h[p]) > 0.0
                    mask_vegsol[p] = true
                else
                    mask_novegsol[p] = true
                end
            end
        end
    end

    # --- Weight reflectance/transmittance by LAI and SAI ---
    np = length(bounds_patch)
    wl = fill!(similar(coszen_col_arr, FT, np), zero(FT))
    ws_arr = fill!(similar(coszen_col_arr, FT, np), zero(FT))
    rho = fill!(similar(coszen_col_arr, FT, np, numrad), zero(FT))
    tau = fill!(similar(coszen_col_arr, FT, np, numrad), zero(FT))

    mpe_FT = FT(mpe)
    # Device copy of the (host) vegsol mask for the kernels; the host BitVector is
    # kept for the two_stream! host solver below.
    mask_vegsol_dev = coszen_col_arr isa Array ? mask_vegsol :
                      copyto!(similar(coszen_col_arr, Bool, length(mask_vegsol)),
                              collect(Bool, mask_vegsol))  # collect: copyto!(device, BitVector) scalar-indexes packed bits
    surfalb_lai_weight!(wl, ws_arr, mask_vegsol_dev, elai, esai, mpe_FT)

    surfalb_rho_tau!(rho, tau, mask_vegsol_dev, patchdata.itype, wl, ws_arr,
                     pftcon_rhol, pftcon_rhos, pftcon_taul, pftcon_taus, mpe_FT, numrad)

    # --- Two-stream calculation ---
    if !use_fates
        two_stream!(surfalb, patchdata, col, canopystate, temperature,
                    waterdiagbulk, coszen_patch_arr, rho, tau,
                    pftcon_xl, mask_vegsol_dev, bounds_patch)

        if use_SSRE
            if nlevcan > 1
                error("use_ssre option was NOT developed with allowance for multi-layer canopy")
            end
            two_stream!(surfalb, patchdata, col, canopystate, temperature,
                        waterdiagbulk, coszen_patch_arr, rho, tau,
                        pftcon_xl, mask_vegsol_dev, bounds_patch; SFonly=true)
        end
    end
    # (FATES canopy radiation stub — not yet ported)

    # --- Non-vegetated patches where coszen > 0 ---
    # mask_novegsol is a BitVector (no device backend); the kernel writes device
    # output arrays, so copy the mask onto the output backend as a Bool vector
    # (no-op-equivalent on CPU). Byte-identical.
    let nv = _SaNovegOut(;
            fabd = surfalb.fabd_patch, fabd_sun = surfalb.fabd_sun_patch, fabd_sha = surfalb.fabd_sha_patch,
            fabi = surfalb.fabi_patch, fabi_sun = surfalb.fabi_sun_patch, fabi_sha = surfalb.fabi_sha_patch,
            ftdd = surfalb.ftdd_patch, ftid = surfalb.ftid_patch, ftii = surfalb.ftii_patch,
            albd = albd, albi = albi, albdSF = surfalb.albdSF_patch, albiSF = surfalb.albiSF_patch)
        be = _kernel_backend(surfalb.fabd_patch)
        if length(bounds_patch) * numrad != 0
            mask_nv = similar(surfalb.fabd_patch, Bool, length(mask_novegsol))
            copyto!(mask_nv, collect(mask_novegsol))
            _sa_noveg_kernel!(be)(nv, mask_nv, patchdata.column,
                                  albgrd, albgri, albsod, albsoi, use_SSRE;
                                  ndrange = (length(bounds_patch), numrad))
            KA.synchronize(be)
        end
    end

    # --- History output variables ---
    let ho = _SaHistColOut(;
            albgrd_hst = surfalb.albgrd_hst_col, albgri_hst = surfalb.albgri_hst_col,
            albgrd_pur_hst = surfalb.albgrd_pur_hst_col, albgri_pur_hst = surfalb.albgri_pur_hst_col,
            albgrd_bc_hst = surfalb.albgrd_bc_hst_col, albgri_bc_hst = surfalb.albgri_bc_hst_col,
            albgrd_oc_hst = surfalb.albgrd_oc_hst_col, albgri_oc_hst = surfalb.albgri_oc_hst_col,
            albgrd_dst_hst = surfalb.albgrd_dst_hst_col, albgri_dst_hst = surfalb.albgri_dst_hst_col,
            albsnd_hst2 = surfalb.albsnd_hst2_col, albsni_hst2 = surfalb.albsni_hst2_col),
        hin = _SaHistColIn(;
            albgrd = albgrd, albgri = albgri,
            albgrd_pur = surfalb.albgrd_pur_col, albgri_pur = surfalb.albgri_pur_col,
            albgrd_bc = surfalb.albgrd_bc_col, albgri_bc = surfalb.albgri_bc_col,
            albgrd_oc = surfalb.albgrd_oc_col, albgri_oc = surfalb.albgri_oc_col,
            albgrd_dst = surfalb.albgrd_dst_col, albgri_dst = surfalb.albgri_dst_col,
            albsnd_hst = surfalb.albsnd_hst_col, albsni_hst = surfalb.albsni_hst_col)
        be = _kernel_backend(albgrd)
        if length(bounds_col) * numrad != 0
            _sa_histcol_kernel!(be)(ho, hin, mask_nourbanc, coszen_col_arr, h2osno_total_arr;
                                    ndrange = (length(bounds_col), numrad))
            KA.synchronize(be)
        end
    end

    _launch!(_sa_histpatch_kernel!, surfalb.albd_hst_patch, surfalb.albi_hst_patch,
             mask_nourbanp, coszen_patch_arr, albd, albi;
             ndrange = (length(bounds_patch), numrad))

    return nothing
end
