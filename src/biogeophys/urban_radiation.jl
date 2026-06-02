# ==========================================================================
# Ported from: src/biogeophys/UrbanRadiationMod.F90
# Calculate solar and longwave radiation, and turbulent fluxes for urban landunit
#
# Public functions:
#   urban_radiation!  — Solar fluxes absorbed and reflected by roof and canyon
#   net_longwave!     — Net longwave radiation for road and walls in urban canyon
# ==========================================================================

# --- Module-level constants ---
const MPE_URBAN_RAD = 1.0e-06   # prevents overflow for division by zero
const SNOEM         = 0.97      # snow emissivity

# ==========================================================================
# GPU kernelization — every per-landunit / per-patch loop in this module runs
# WHOLE-FUNCTION on the device. Each urban landunit is fully independent in the
# canyon longwave solve, so the ~80 Fortran scratch vectors of `net_longwave`
# collapse to thread-local SCALARS in a single per-landunit kernel (no scratch
# arrays cross the launch). State arrays are bundled into immutable device-view
# structs (Adapt.@adapt_structure'd) so each launch stays under Metal's ~31-arg
# cap; the kernel bodies mirror the original Julia paths so CPU stays byte-
# identical. Scalar constants (SB) are passed at the output's element type so no
# Float64 reaches the Float32-only Metal backend.
# ==========================================================================

# --- net_longwave! device-view bundles ---------------------------------------

# Inputs read by the canyon longwave solve (landunit-level).
Base.@kwdef struct _ULWInDV{V}
    canyon_hwr::V; wtroad_perv::V; lwdown::V
    em_roof::V; em_improad::V; em_perroad::V; em_wall::V
    t_roof::V; t_improad::V; t_perroad::V; t_sunwall::V; t_shadewall::V
    vf_sr::V; vf_wr::V; vf_sw::V; vf_rw::V; vf_ww::V
end
Adapt.@adapt_structure _ULWInDV

# Outputs written by the canyon longwave solve (landunit-level).
Base.@kwdef struct _ULWOutDV{V}
    lwnet_roof::V; lwnet_improad::V; lwnet_perroad::V
    lwnet_sunwall::V; lwnet_shadewall::V; lwnet_canyon::V
    lwup_roof::V; lwup_improad::V; lwup_perroad::V
    lwup_sunwall::V; lwup_shadewall::V; lwup_canyon::V
end
Adapt.@adapt_structure _ULWOutDV

# Net longwave radiation for road, both walls, roof, and canyon totals — one
# fully independent thread per urban landunit. All Fortran scratch is local.
@kernel function _urad_netlw_kernel!(@Const(_out), @Const(mask), inp, out,
                                     imin::Int, imax::Int, niter::Int, SBc)
    l = @index(Global)
    @inbounds if imin <= l <= imax && mask[l]
        T = eltype(out.lwnet_roof)
        z = zero(T); o = one(T)

        canyon_hwr  = inp.canyon_hwr[l]
        wtroad_perv = inp.wtroad_perv[l]
        wtroad_imperv = o - wtroad_perv
        lwdown      = inp.lwdown[l]
        em_improad  = inp.em_improad[l]
        em_perroad  = inp.em_perroad[l]
        em_wall     = inp.em_wall[l]
        t_improad   = inp.t_improad[l]
        t_perroad   = inp.t_perroad[l]
        t_sunwall   = inp.t_sunwall[l]
        t_shadewall = inp.t_shadewall[l]
        vf_sr = inp.vf_sr[l]; vf_wr = inp.vf_wr[l]
        vf_sw = inp.vf_sw[l]; vf_rw = inp.vf_rw[l]; vf_ww = inp.vf_ww[l]

        # Atmospheric longwave radiation incident on walls and road in canyon
        lwdown_road      = lwdown * vf_sr
        lwdown_sunwall   = lwdown * vf_sw
        lwdown_shadewall = lwdown * vf_sw

        # Initial absorption, reflection, and emission
        road_a = z; road_r = z; road_e = z

        improad_a           =       em_improad  * lwdown_road
        improad_r           = (o - em_improad) * lwdown_road
        improad_e           = em_improad * SBc * (t_improad^4)
        road_a += improad_a * wtroad_imperv
        road_r += improad_r * wtroad_imperv
        road_e += improad_e * wtroad_imperv

        perroad_a           =       em_perroad  * lwdown_road
        perroad_r           = (o - em_perroad) * lwdown_road
        perroad_e           = em_perroad * SBc * (t_perroad^4)
        road_a += perroad_a * wtroad_perv
        road_r += perroad_r * wtroad_perv
        road_e += perroad_e * wtroad_perv

        road_r_sunwall      = road_r * vf_wr
        road_r_shadewall    = road_r * vf_wr
        road_e_sunwall      = road_e * vf_wr
        road_e_shadewall    = road_e * vf_wr

        sunwall_a           = em_wall * lwdown_sunwall
        sunwall_r           = (o - em_wall) * lwdown_sunwall
        sunwall_r_road      = sunwall_r * vf_rw
        sunwall_r_shadewall = sunwall_r * vf_ww
        sunwall_e           = em_wall * SBc * (t_sunwall^4)
        sunwall_e_road      = sunwall_e * vf_rw
        sunwall_e_shadewall = sunwall_e * vf_ww

        shadewall_a           = em_wall * lwdown_shadewall
        shadewall_r           = (o - em_wall) * lwdown_shadewall
        shadewall_r_road      = shadewall_r * vf_rw
        shadewall_r_sunwall   = shadewall_r * vf_ww
        shadewall_e           = em_wall * SBc * (t_shadewall^4)
        shadewall_e_road      = shadewall_e * vf_rw
        shadewall_e_sunwall   = shadewall_e * vf_ww

        # Initialize sum of net and upward longwave radiation
        lwnet_improad   = improad_e   - improad_a
        lwnet_perroad   = perroad_e   - perroad_a
        lwnet_sunwall   = sunwall_e   - sunwall_a
        lwnet_shadewall = shadewall_e - shadewall_a

        lwup_improad   = improad_r   * vf_sr + improad_e   * vf_sr
        lwup_perroad   = perroad_r   * vf_sr + perroad_e   * vf_sr
        lwup_sunwall   = sunwall_r   * vf_sw + sunwall_e   * vf_sw
        lwup_shadewall = shadewall_r * vf_sw + shadewall_e * vf_sw

        # Multiple reflection iterations within canyon
        for iter in 1:niter
            # step (1): absorption and reflection
            lwtot = (sunwall_r_road + sunwall_e_road +
                     shadewall_r_road + shadewall_e_road) * canyon_hwr
            road_a = z; road_r = z
            improad_r = (o - em_improad) * lwtot
            improad_a =       em_improad  * lwtot
            road_a += improad_a * wtroad_imperv
            road_r += improad_r * wtroad_imperv
            perroad_r = (o - em_perroad) * lwtot
            perroad_a =       em_perroad  * lwtot
            road_a += perroad_a * wtroad_perv
            road_r += perroad_r * wtroad_perv

            lwtot = (road_r_sunwall + road_e_sunwall) / canyon_hwr +
                    (shadewall_r_sunwall + shadewall_e_sunwall)
            sunwall_a =       em_wall  * lwtot
            sunwall_r = (o - em_wall) * lwtot

            lwtot = (road_r_shadewall + road_e_shadewall) / canyon_hwr +
                    (sunwall_r_shadewall + sunwall_e_shadewall)
            shadewall_a =       em_wall  * lwtot
            shadewall_r = (o - em_wall) * lwtot

            # Zero emission terms after first iteration
            sunwall_e_road      = z
            shadewall_e_road    = z
            road_e_sunwall      = z
            shadewall_e_sunwall = z
            road_e_shadewall    = z
            sunwall_e_shadewall = z

            # step (2): add net longwave for ith reflection to total
            lwnet_improad   -= improad_a
            lwnet_perroad   -= perroad_a
            lwnet_sunwall   -= sunwall_a
            lwnet_shadewall -= shadewall_a

            # step (3): distribute reflected radiation
            improad_r_sky       = improad_r * vf_sr
            road_r_sunwall      = road_r * vf_wr
            road_r_shadewall    = road_r * vf_wr
            perroad_r_sky       = perroad_r * vf_sr

            sunwall_r_sky       = sunwall_r * vf_sw
            sunwall_r_road      = sunwall_r * vf_rw
            sunwall_r_shadewall = sunwall_r * vf_ww

            shadewall_r_sky     = shadewall_r * vf_sw
            shadewall_r_road    = shadewall_r * vf_rw
            shadewall_r_sunwall = shadewall_r * vf_ww

            # step (4): add upward longwave radiation to sky
            lwup_improad   += improad_r_sky
            lwup_perroad   += perroad_r_sky
            lwup_sunwall   += sunwall_r_sky
            lwup_shadewall += shadewall_r_sky

            # step (5): convergence check
            crit = max(road_a, sunwall_a, shadewall_a)
            crit < oftype(crit, 0.001) && break
        end

        # Total net longwave radiation for canyon
        lwnet_canyon = lwnet_improad * wtroad_imperv +
                       lwnet_perroad * wtroad_perv +
                       (lwnet_sunwall + lwnet_shadewall) * canyon_hwr

        # Total emitted longwave for canyon
        lwup_canyon = lwup_improad * wtroad_imperv +
                      lwup_perroad * wtroad_perv +
                      (lwup_sunwall + lwup_shadewall) * canyon_hwr

        # Net longwave radiation for roof
        em_roof = inp.em_roof[l]; t_roof = inp.t_roof[l]
        lwup_roof  = em_roof * SBc * (t_roof^4) + (o - em_roof) * lwdown
        lwnet_roof = lwup_roof - lwdown

        out.lwnet_improad[l]   = lwnet_improad
        out.lwnet_perroad[l]   = lwnet_perroad
        out.lwnet_sunwall[l]   = lwnet_sunwall
        out.lwnet_shadewall[l] = lwnet_shadewall
        out.lwnet_canyon[l]    = lwnet_canyon
        out.lwnet_roof[l]      = lwnet_roof
        out.lwup_improad[l]    = lwup_improad
        out.lwup_perroad[l]    = lwup_perroad
        out.lwup_sunwall[l]    = lwup_sunwall
        out.lwup_shadewall[l]  = lwup_shadewall
        out.lwup_canyon[l]     = lwup_canyon
        out.lwup_roof[l]       = lwup_roof
    end
end

"""
    net_longwave!(canyon_hwr, wtroad_perv, lwdown, em_roof, em_improad,
                  em_perroad, em_wall, t_roof, t_improad, t_perroad,
                  t_sunwall, t_shadewall,
                  lwnet_roof, lwnet_improad, lwnet_perroad,
                  lwnet_sunwall, lwnet_shadewall, lwnet_canyon,
                  lwup_roof, lwup_improad, lwup_perroad,
                  lwup_sunwall, lwup_shadewall, lwup_canyon,
                  urbanparams, mask_urbanl, bounds)

Net longwave radiation for road and both walls in urban canyon allowing for
multiple reflection. Also net longwave radiation for urban roof.

Ported from `net_longwave` in `UrbanRadiationMod.F90`.

# Arguments
- `canyon_hwr`    : ratio of building height to street width (landunit)
- `wtroad_perv`   : weight of pervious road wrt total road (landunit)
- `lwdown`        : atmospheric longwave radiation (W/m²) (landunit)
- `em_roof`       : roof emissivity with snow effects (landunit)
- `em_improad`    : impervious road emissivity with snow effects (landunit)
- `em_perroad`    : pervious road emissivity with snow effects (landunit)
- `em_wall`       : wall emissivity (landunit)
- `t_roof`        : roof temperature (K) (landunit)
- `t_improad`     : impervious road temperature (K) (landunit)
- `t_perroad`     : pervious road temperature (K) (landunit)
- `t_sunwall`     : sunlit wall temperature (K) (landunit)
- `t_shadewall`   : shaded wall temperature (K) (landunit)
- `lwnet_*`       : [output] net (outgoing-incoming) longwave radiation (W/m²)
- `lwup_*`        : [output] upward longwave radiation (W/m²)
- `urbanparams`   : UrbanParamsData (view factors)
- `mask_urbanl`   : BitVector mask for urban landunits
- `bounds`        : landunit index range
"""
function net_longwave!(canyon_hwr::AbstractVector{<:Real},
                        wtroad_perv::AbstractVector{<:Real},
                        lwdown::AbstractVector{<:Real},
                        em_roof::AbstractVector{<:Real},
                        em_improad::AbstractVector{<:Real},
                        em_perroad::AbstractVector{<:Real},
                        em_wall::AbstractVector{<:Real},
                        t_roof::AbstractVector{<:Real},
                        t_improad::AbstractVector{<:Real},
                        t_perroad::AbstractVector{<:Real},
                        t_sunwall::AbstractVector{<:Real},
                        t_shadewall::AbstractVector{<:Real},
                        lwnet_roof::AbstractVector{<:Real},
                        lwnet_improad::AbstractVector{<:Real},
                        lwnet_perroad::AbstractVector{<:Real},
                        lwnet_sunwall::AbstractVector{<:Real},
                        lwnet_shadewall::AbstractVector{<:Real},
                        lwnet_canyon::AbstractVector{<:Real},
                        lwup_roof::AbstractVector{<:Real},
                        lwup_improad::AbstractVector{<:Real},
                        lwup_perroad::AbstractVector{<:Real},
                        lwup_sunwall::AbstractVector{<:Real},
                        lwup_shadewall::AbstractVector{<:Real},
                        lwup_canyon::AbstractVector{<:Real},
                        urbanparams::UrbanParamsData,
                        mask_urbanl::AbstractVector{Bool},
                        bounds::UnitRange{Int})

    n = 50  # number of iterations

    inp = _ULWInDV(; canyon_hwr, wtroad_perv, lwdown,
                     em_roof, em_improad, em_perroad, em_wall,
                     t_roof, t_improad, t_perroad, t_sunwall, t_shadewall,
                     vf_sr = urbanparams.vf_sr, vf_wr = urbanparams.vf_wr,
                     vf_sw = urbanparams.vf_sw, vf_rw = urbanparams.vf_rw,
                     vf_ww = urbanparams.vf_ww)
    out = _ULWOutDV(; lwnet_roof, lwnet_improad, lwnet_perroad,
                      lwnet_sunwall, lwnet_shadewall, lwnet_canyon,
                      lwup_roof, lwup_improad, lwup_perroad,
                      lwup_sunwall, lwup_shadewall, lwup_canyon)

    # SB at the output element type → no Float64 reaches a Float32-only backend.
    SBc = convert(eltype(lwnet_roof), SB)
    _launch!(_urad_netlw_kernel!, lwnet_roof, mask_urbanl, inp, out,
             first(bounds), last(bounds), n, SBc; ndrange = length(lwnet_roof))
    return nothing
end

# ==========================================================================
# urban_radiation! kernels — three fully-independent per-landunit / per-patch
# passes (SPVAL fill, snow-effect input forcing, history/atm output).
# ==========================================================================

# Pass A: non-urban restart fields set to SPVAL (per-landunit, masked).
Base.@kwdef struct _USpvalDV{M}
    sabs_roof_dir_lun::M; sabs_roof_dif_lun::M
    sabs_sunwall_dir_lun::M; sabs_sunwall_dif_lun::M
    sabs_shadewall_dir_lun::M; sabs_shadewall_dif_lun::M
    sabs_improad_dir_lun::M; sabs_improad_dif_lun::M
    sabs_perroad_dir_lun::M; sabs_perroad_dif_lun::M
end
Adapt.@adapt_structure _USpvalDV

@kernel function _urad_spval_kernel!(@Const(_out), @Const(mask), sa,
                                     imin::Int, imax::Int, spval)
    l = @index(Global)
    @inbounds if imin <= l <= imax && mask[l]
        for ib in axes(sa.sabs_roof_dir_lun, 2)
            sa.sabs_roof_dir_lun[l, ib]      = spval
            sa.sabs_roof_dif_lun[l, ib]      = spval
            sa.sabs_sunwall_dir_lun[l, ib]   = spval
            sa.sabs_sunwall_dif_lun[l, ib]   = spval
            sa.sabs_shadewall_dir_lun[l, ib] = spval
            sa.sabs_shadewall_dif_lun[l, ib] = spval
            sa.sabs_improad_dir_lun[l, ib]   = spval
            sa.sabs_improad_dif_lun[l, ib]   = spval
            sa.sabs_perroad_dir_lun[l, ib]   = spval
            sa.sabs_perroad_dif_lun[l, ib]   = spval
        end
    end
end

# Pass B: snow-effect input forcing (per-landunit; internal sequential column
# loop over coli:colf). Writes per-landunit temps / emissivities / lwdown.
Base.@kwdef struct _UInOutDV{V}
    t_roof_l::V; t_improad_l::V; t_perroad_l::V; t_sunwall_l::V; t_shadewall_l::V
    em_roof_s::V; em_improad_s::V; em_perroad_s::V; lwdown_l::V
end
Adapt.@adapt_structure _UInOutDV

Base.@kwdef struct _UInSrcDV{V,VI}
    em_roof::V; em_improad::V; em_perroad::V
    t_grnd_col::V; frac_sno_col::V; forc_lwrad::V
    itype::VI; coli::VI; colf::VI; gridcell::VI
end
Adapt.@adapt_structure _UInSrcDV

@kernel function _urad_input_kernel!(@Const(_out), @Const(mask), o, s,
        imin::Int, imax::Int, tdef, snoem,
        icol_roof::Int, icol_improad::Int, icol_perv::Int,
        icol_sun::Int, icol_shade::Int)
    l = @index(Global)
    @inbounds if imin <= l <= imax && mask[l]
        T = eltype(o.t_roof_l)
        on = one(T)
        o.t_roof_l[l]      = tdef
        o.t_sunwall_l[l]   = tdef
        o.t_shadewall_l[l] = tdef
        o.t_improad_l[l]   = tdef
        o.t_perroad_l[l]   = tdef
        o.em_roof_s[l]    = s.em_roof[l]
        o.em_improad_s[l] = s.em_improad[l]
        o.em_perroad_s[l] = s.em_perroad[l]
        for c in s.coli[l]:s.colf[l]
            it = s.itype[c]
            fsno = s.frac_sno_col[c]
            if it == icol_roof
                o.t_roof_l[l]  = s.t_grnd_col[c]
                o.em_roof_s[l] = s.em_roof[l] * (on - fsno) + snoem * fsno
            elseif it == icol_improad
                o.t_improad_l[l]  = s.t_grnd_col[c]
                o.em_improad_s[l] = s.em_improad[l] * (on - fsno) + snoem * fsno
            elseif it == icol_perv
                o.t_perroad_l[l]  = s.t_grnd_col[c]
                o.em_perroad_s[l] = s.em_perroad[l] * (on - fsno) + snoem * fsno
            elseif it == icol_sun
                o.t_sunwall_l[l] = s.t_grnd_col[c]
            elseif it == icol_shade
                o.t_shadewall_l[l] = s.t_grnd_col[c]
            end
        end
        o.lwdown_l[l] = s.forc_lwrad[s.gridcell[l]]
    end
end

# Pass C: history / atm-communication output (per-patch, masked, branchy on
# column type). Bundles the SolarAbsorbedData sabs matrices + outputs, the
# EnergyFluxData outputs, and the per-landunit lwup/lwnet results.
Base.@kwdef struct _UHistSolDV{V,M}
    sabg_patch::V; sabv_patch::V; fsa_patch::V; fsa_u_patch::V
    sabs_roof_dir_lun::M; sabs_roof_dif_lun::M
    sabs_sunwall_dir_lun::M; sabs_sunwall_dif_lun::M
    sabs_shadewall_dir_lun::M; sabs_shadewall_dif_lun::M
    sabs_improad_dir_lun::M; sabs_improad_dif_lun::M
    sabs_perroad_dir_lun::M; sabs_perroad_dif_lun::M
end
Adapt.@adapt_structure _UHistSolDV

Base.@kwdef struct _UHistLwDV{V}
    lwup_roof::V; lwup_sunwall::V; lwup_shadewall::V; lwup_perroad::V; lwup_improad::V
    lwnet_roof::V; lwnet_sunwall::V; lwnet_shadewall::V; lwnet_perroad::V; lwnet_improad::V
    eflx_lwrad_out_patch::V; eflx_lwrad_net_patch::V; eflx_lwrad_net_u_patch::V
end
Adapt.@adapt_structure _UHistLwDV

@kernel function _urad_history_kernel!(@Const(_out), @Const(mask), sol, lw,
        @Const(column), @Const(landunit), @Const(gridcell), @Const(itype),
        @Const(forc_solad), @Const(forc_solai),
        icol_roof::Int, icol_sun::Int, icol_shade::Int,
        icol_perv::Int, icol_improad::Int)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(sol.sabg_patch)
        c = column[p]; l = landunit[p]; g = gridcell[p]
        it = itype[c]
        if it == icol_roof
            lw.eflx_lwrad_out_patch[p]   = lw.lwup_roof[l]
            lw.eflx_lwrad_net_patch[p]   = lw.lwnet_roof[l]
            lw.eflx_lwrad_net_u_patch[p] = lw.lwnet_roof[l]
            sol.sabg_patch[p] = sol.sabs_roof_dir_lun[l, 1] * forc_solad[g, 1] +
                sol.sabs_roof_dif_lun[l, 1] * forc_solai[g, 1] +
                sol.sabs_roof_dir_lun[l, 2] * forc_solad[g, 2] +
                sol.sabs_roof_dif_lun[l, 2] * forc_solai[g, 2]
        elseif it == icol_sun
            lw.eflx_lwrad_out_patch[p]   = lw.lwup_sunwall[l]
            lw.eflx_lwrad_net_patch[p]   = lw.lwnet_sunwall[l]
            lw.eflx_lwrad_net_u_patch[p] = lw.lwnet_sunwall[l]
            sol.sabg_patch[p] = sol.sabs_sunwall_dir_lun[l, 1] * forc_solad[g, 1] +
                sol.sabs_sunwall_dif_lun[l, 1] * forc_solai[g, 1] +
                sol.sabs_sunwall_dir_lun[l, 2] * forc_solad[g, 2] +
                sol.sabs_sunwall_dif_lun[l, 2] * forc_solai[g, 2]
        elseif it == icol_shade
            lw.eflx_lwrad_out_patch[p]   = lw.lwup_shadewall[l]
            lw.eflx_lwrad_net_patch[p]   = lw.lwnet_shadewall[l]
            lw.eflx_lwrad_net_u_patch[p] = lw.lwnet_shadewall[l]
            sol.sabg_patch[p] = sol.sabs_shadewall_dir_lun[l, 1] * forc_solad[g, 1] +
                sol.sabs_shadewall_dif_lun[l, 1] * forc_solai[g, 1] +
                sol.sabs_shadewall_dir_lun[l, 2] * forc_solad[g, 2] +
                sol.sabs_shadewall_dif_lun[l, 2] * forc_solai[g, 2]
        elseif it == icol_perv
            lw.eflx_lwrad_out_patch[p]   = lw.lwup_perroad[l]
            lw.eflx_lwrad_net_patch[p]   = lw.lwnet_perroad[l]
            lw.eflx_lwrad_net_u_patch[p] = lw.lwnet_perroad[l]
            sol.sabg_patch[p] = sol.sabs_perroad_dir_lun[l, 1] * forc_solad[g, 1] +
                sol.sabs_perroad_dif_lun[l, 1] * forc_solai[g, 1] +
                sol.sabs_perroad_dir_lun[l, 2] * forc_solad[g, 2] +
                sol.sabs_perroad_dif_lun[l, 2] * forc_solai[g, 2]
        elseif it == icol_improad
            lw.eflx_lwrad_out_patch[p]   = lw.lwup_improad[l]
            lw.eflx_lwrad_net_patch[p]   = lw.lwnet_improad[l]
            lw.eflx_lwrad_net_u_patch[p] = lw.lwnet_improad[l]
            sol.sabg_patch[p] = sol.sabs_improad_dir_lun[l, 1] * forc_solad[g, 1] +
                sol.sabs_improad_dif_lun[l, 1] * forc_solai[g, 1] +
                sol.sabs_improad_dir_lun[l, 2] * forc_solad[g, 2] +
                sol.sabs_improad_dif_lun[l, 2] * forc_solai[g, 2]
        end
        sol.sabv_patch[p]  = zero(T)
        sol.fsa_patch[p]   = sol.sabv_patch[p] + sol.sabg_patch[p]
        sol.fsa_u_patch[p] = sol.fsa_patch[p]
    end
end

"""
    urban_radiation!(solarabs, energyflux, col, lun, pch, urbanparams, temperature,
                     waterdiag, forc_lwrad, forc_solad, forc_solai,
                     mask_nourbanl, mask_urbanl, mask_urbanc, mask_urbanp,
                     bounds_lun, bounds_patch)

Solar fluxes absorbed and reflected by roof and canyon (walls, road).
Also net and upward longwave fluxes.

Ported from `UrbanRadiation` in `UrbanRadiationMod.F90`.

# Arguments
- `solarabs`    : SolarAbsorbedData (input/output)
- `energyflux`  : EnergyFluxData (output: eflx_lwrad_out, eflx_lwrad_net, eflx_lwrad_net_u)
- `col`         : ColumnData (input: itype)
- `lun`         : LandunitData (input: coli, colf, canyon_hwr, wtroad_perv, gridcell)
- `pch`         : PatchData (input: column, landunit, gridcell)
- `urbanparams` : UrbanParamsData (input: em_roof, em_improad, em_perroad, em_wall, view factors)
- `temperature` : TemperatureData (input: t_grnd_col)
- `waterdiag`   : WaterDiagnosticBulkData (input: frac_sno_col)
- `forc_lwrad`  : downward longwave radiation per gridcell (Vector)
- `forc_solad`  : direct beam radiation per gridcell (Matrix: ngrc × numrad)
- `forc_solai`  : diffuse radiation per gridcell (Matrix: ngrc × numrad)
- `mask_nourbanl` : BitVector mask for non-urban landunits
- `mask_urbanl`   : BitVector mask for urban landunits
- `mask_urbanc`   : BitVector mask for urban columns (unused, kept for interface compat)
- `mask_urbanp`   : BitVector mask for urban patches
- `bounds_lun`    : landunit index range
- `bounds_patch`  : patch index range
"""
function urban_radiation!(solarabs::SolarAbsorbedData,
                            energyflux::EnergyFluxData,
                            col::ColumnData,
                            lun::LandunitData,
                            pch::PatchData,
                            urbanparams::UrbanParamsData,
                            temperature::TemperatureData,
                            waterdiag::WaterDiagnosticBulkData,
                            forc_lwrad::AbstractVector{<:Real},
                            forc_solad::AbstractMatrix{<:Real},
                            forc_solai::AbstractMatrix{<:Real},
                            mask_nourbanl::AbstractVector{Bool},
                            mask_urbanl::AbstractVector{Bool},
                            mask_urbanc::AbstractVector{Bool},
                            mask_urbanp::AbstractVector{Bool},
                            bounds_lun::UnitRange{Int},
                            bounds_patch::UnitRange{Int})

    nl = last(bounds_lun)

    # Landunit-level scratch — device-resident (mirror the output backend) so
    # the kernels below run wholly on-device. Use `similar`+`fill!`, not zeros().
    proto = temperature.t_grnd_col
    FT = eltype(proto)
    mkl() = (s = similar(proto, nl); fill!(s, zero(FT)); s)
    lwnet_roof      = mkl(); lwnet_improad   = mkl(); lwnet_perroad   = mkl()
    lwnet_sunwall   = mkl(); lwnet_shadewall = mkl(); lwnet_canyon    = mkl()
    lwup_roof       = mkl(); lwup_improad    = mkl(); lwup_perroad    = mkl()
    lwup_sunwall    = mkl(); lwup_shadewall  = mkl(); lwup_canyon     = mkl()
    t_roof_l        = mkl(); t_improad_l     = mkl(); t_perroad_l     = mkl()
    t_sunwall_l     = mkl(); t_shadewall_l   = mkl(); lwdown_l        = mkl()
    em_roof_s       = mkl(); em_improad_s    = mkl(); em_perroad_s    = mkl()

    imin = first(bounds_lun); imax = last(bounds_lun)

    # Pass A: non-urban restart fields set to SPVAL (per-landunit, masked).
    spval_dv = _USpvalDV(;
        sabs_roof_dir_lun = solarabs.sabs_roof_dir_lun,
        sabs_roof_dif_lun = solarabs.sabs_roof_dif_lun,
        sabs_sunwall_dir_lun = solarabs.sabs_sunwall_dir_lun,
        sabs_sunwall_dif_lun = solarabs.sabs_sunwall_dif_lun,
        sabs_shadewall_dir_lun = solarabs.sabs_shadewall_dir_lun,
        sabs_shadewall_dif_lun = solarabs.sabs_shadewall_dif_lun,
        sabs_improad_dir_lun = solarabs.sabs_improad_dir_lun,
        sabs_improad_dif_lun = solarabs.sabs_improad_dif_lun,
        sabs_perroad_dir_lun = solarabs.sabs_perroad_dir_lun,
        sabs_perroad_dif_lun = solarabs.sabs_perroad_dif_lun)
    _launch!(_urad_spval_kernel!, solarabs.sabg_patch, mask_nourbanl, spval_dv,
             imin, imax, FT(SPVAL); ndrange = nl)

    # Pass B: snow-effect input forcing (per-landunit; internal column loop).
    in_o = _UInOutDV(; t_roof_l, t_improad_l, t_perroad_l, t_sunwall_l, t_shadewall_l,
                       em_roof_s, em_improad_s, em_perroad_s, lwdown_l)
    in_s = _UInSrcDV(; em_roof = urbanparams.em_roof, em_improad = urbanparams.em_improad,
                       em_perroad = urbanparams.em_perroad,
                       t_grnd_col = temperature.t_grnd_col,
                       frac_sno_col = waterdiag.frac_sno_col, forc_lwrad = forc_lwrad,
                       itype = col.itype, coli = lun.coli, colf = lun.colf,
                       gridcell = lun.gridcell)
    _launch!(_urad_input_kernel!, lwdown_l, mask_urbanl, in_o, in_s,
             imin, imax, FT(19.0 + TFRZ), FT(SNOEM),
             ICOL_ROOF, ICOL_ROAD_IMPERV, ICOL_ROAD_PERV, ICOL_SUNWALL, ICOL_SHADEWALL;
             ndrange = nl)

    # Net longwave radiation for road and both walls in urban canyon.
    net_longwave!(lun.canyon_hwr, lun.wtroad_perv,
                   lwdown_l, em_roof_s, em_improad_s, em_perroad_s,
                   urbanparams.em_wall,
                   t_roof_l, t_improad_l, t_perroad_l,
                   t_sunwall_l, t_shadewall_l,
                   lwnet_roof, lwnet_improad, lwnet_perroad,
                   lwnet_sunwall, lwnet_shadewall, lwnet_canyon,
                   lwup_roof, lwup_improad, lwup_perroad,
                   lwup_sunwall, lwup_shadewall, lwup_canyon,
                   urbanparams, mask_urbanl, bounds_lun)

    # Pass C: history output and communication with atm (per-patch, masked).
    hist_sol = _UHistSolDV(;
        sabg_patch = solarabs.sabg_patch, sabv_patch = solarabs.sabv_patch,
        fsa_patch = solarabs.fsa_patch, fsa_u_patch = solarabs.fsa_u_patch,
        sabs_roof_dir_lun = solarabs.sabs_roof_dir_lun,
        sabs_roof_dif_lun = solarabs.sabs_roof_dif_lun,
        sabs_sunwall_dir_lun = solarabs.sabs_sunwall_dir_lun,
        sabs_sunwall_dif_lun = solarabs.sabs_sunwall_dif_lun,
        sabs_shadewall_dir_lun = solarabs.sabs_shadewall_dir_lun,
        sabs_shadewall_dif_lun = solarabs.sabs_shadewall_dif_lun,
        sabs_improad_dir_lun = solarabs.sabs_improad_dir_lun,
        sabs_improad_dif_lun = solarabs.sabs_improad_dif_lun,
        sabs_perroad_dir_lun = solarabs.sabs_perroad_dir_lun,
        sabs_perroad_dif_lun = solarabs.sabs_perroad_dif_lun)
    hist_lw = _UHistLwDV(;
        lwup_roof, lwup_sunwall, lwup_shadewall, lwup_perroad, lwup_improad,
        lwnet_roof, lwnet_sunwall, lwnet_shadewall, lwnet_perroad, lwnet_improad,
        eflx_lwrad_out_patch = energyflux.eflx_lwrad_out_patch,
        eflx_lwrad_net_patch = energyflux.eflx_lwrad_net_patch,
        eflx_lwrad_net_u_patch = energyflux.eflx_lwrad_net_u_patch)
    _launch!(_urad_history_kernel!, solarabs.sabg_patch, mask_urbanp, hist_sol, hist_lw,
             pch.column, pch.landunit, pch.gridcell, col.itype,
             forc_solad, forc_solai,
             ICOL_ROOF, ICOL_SUNWALL, ICOL_SHADEWALL, ICOL_ROAD_PERV, ICOL_ROAD_IMPERV;
             ndrange = length(solarabs.sabg_patch))

    return nothing
end
