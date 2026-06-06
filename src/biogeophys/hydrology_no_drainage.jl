# ==========================================================================
# Ported from: src/biogeophys/HydrologyNoDrainageMod.F90 (~727 lines)
# Calculate snow and soil hydrology without drainage.
#
# Public functions:
#   update_snow_persistence!               — Update snow persistence counters
#   accumulate_snow_ice_liq!               — Sum snowice/snowliq over snow layers
#   update_snowdp!                         — Column-average snow depth
#   compute_snow_internal_temperature!     — Snow internal temperature (t_sno_mul_mss)
#   update_ground_and_soil_temperatures!   — Ground temp, tsl, t_soi_10cm, tsoi17
#   update_h2osoi_vol!                     — Volumetric soil water content
#   update_soilpsi!                        — Soil water potential (MPa)
#   update_smp_l!                          — Matric potential for history/ch4
#   compute_wf!                            — Soil water fraction of WHC (top 0.05m)
#   compute_wf2!                           — Soil water fraction of WHC (top 0.17m)
#   update_snow_top_layer_diagnostics!     — Top-layer snow diagnostics
#   hydrology_no_drainage!                 — Main orchestrator (calls sub-functions)
#   calc_and_withdraw_irrigation_fluxes!   — Irrigation withdrawal (stub)
#   handle_new_snow!                       — Handle new snow falling on ground
# ==========================================================================

# =========================================================================
# KernelAbstractions kernels for simple independent per-element loops.
# (`@kernel`, `@index`, `@Const`, `_launch!` come from infrastructure/kernels.jl.)
# Each kernel mirrors the scalar loop it replaces exactly (same masks, offsets,
# branches) to preserve Fortran parity and AD behavior.
# =========================================================================

# --- update_snow_persistence! : per-column, fully independent ---
@kernel function _hydrond_snow_persistence_kernel!(snow_persistence, @Const(mask_snow),
                                                   @Const(mask_nosnow), dtime, cmin::Int)
    i = @index(Global)
    @inbounds begin
        c = cmin + i - 1
        if mask_snow[c]
            snow_persistence[c] = snow_persistence[c] + dtime
        elseif mask_nosnow[c]
            snow_persistence[c] = 0.0
        end
    end
end
hydrond_snow_persistence!(snow_persistence, mask_snow, mask_nosnow, dtime, bounds) =
    _launch!(_hydrond_snow_persistence_kernel!, snow_persistence, mask_snow, mask_nosnow,
             eltype(snow_persistence)(dtime), first(bounds); ndrange = length(bounds))

# --- update_snowdp! : per-column, fully independent ---
@kernel function _hydrond_snowdp_kernel!(snowdp, @Const(snow_depth), @Const(frac_sno_eff),
                                         cmin::Int)
    i = @index(Global)
    @inbounds begin
        c = cmin + i - 1
        snowdp[c] = snow_depth[c] * frac_sno_eff[c]
    end
end
hydrond_snowdp!(snowdp, snow_depth, frac_sno_eff, bounds) =
    _launch!(_hydrond_snowdp_kernel!, snowdp, snow_depth, frac_sno_eff,
             first(bounds); ndrange = length(bounds))

# --- update_h2osoi_vol! : 2D (column, layer), each [c,j] independent ---
# Non-urban (non sunwall/shadewall/roof) branch.
@kernel function _hydrond_h2osoi_vol_nonurban_kernel!(h2osoi_vol, @Const(h2osoi_liq),
        @Const(h2osoi_ice), @Const(dz), @Const(col_itype), @Const(mask_nolake),
        joff::Int, cmin::Int, denh2o, denice)
    i, j = @index(Global, NTuple)
    @inbounds begin
        c = cmin + i - 1
        if mask_nolake[c] &&
           col_itype[c] != ICOL_SUNWALL && col_itype[c] != ICOL_SHADEWALL &&
           col_itype[c] != ICOL_ROOF
            h2osoi_vol[c, j] = h2osoi_liq[c, j + joff] / (dz[c, j + joff] * denh2o) +
                               h2osoi_ice[c, j + joff] / (dz[c, j + joff] * denice)
        end
    end
end
hydrond_h2osoi_vol_nonurban!(h2osoi_vol, h2osoi_liq, h2osoi_ice, dz, col_itype,
        mask_nolake, joff, bounds, nlevgrnd) =
    _launch!(_hydrond_h2osoi_vol_nonurban_kernel!, h2osoi_vol, h2osoi_liq, h2osoi_ice,
             dz, col_itype, mask_nolake, joff, first(bounds),
             eltype(h2osoi_vol)(DENH2O), eltype(h2osoi_vol)(DENICE);
             ndrange = (length(bounds), nlevgrnd))

# Urban (sunwall/shadewall/roof) branch.
@kernel function _hydrond_h2osoi_vol_urban_kernel!(h2osoi_vol, @Const(h2osoi_liq),
        @Const(h2osoi_ice), @Const(dz), @Const(col_itype), @Const(mask_urban),
        joff::Int, cmin::Int, denh2o, denice)
    i, j = @index(Global, NTuple)
    @inbounds begin
        c = cmin + i - 1
        if mask_urban[c] &&
           (col_itype[c] == ICOL_SUNWALL || col_itype[c] == ICOL_SHADEWALL ||
            col_itype[c] == ICOL_ROOF)
            h2osoi_vol[c, j] = h2osoi_liq[c, j + joff] / (dz[c, j + joff] * denh2o) +
                               h2osoi_ice[c, j + joff] / (dz[c, j + joff] * denice)
        end
    end
end
hydrond_h2osoi_vol_urban!(h2osoi_vol, h2osoi_liq, h2osoi_ice, dz, col_itype,
        mask_urban, joff, bounds, nlevurb) =
    _launch!(_hydrond_h2osoi_vol_urban_kernel!, h2osoi_vol, h2osoi_liq, h2osoi_ice,
             dz, col_itype, mask_urban, joff, first(bounds),
             eltype(h2osoi_vol)(DENH2O), eltype(h2osoi_vol)(DENICE);
             ndrange = (length(bounds), nlevurb))

# --- update_soilpsi! : 2D (column, layer), each [c,j] independent ---
@kernel function _hydrond_soilpsi_kernel!(soilpsi, @Const(h2osoi_liq), @Const(dz),
        @Const(watsat), @Const(sucsat), @Const(bsw), @Const(mask_hydrology),
        joff::Int, cmin::Int, denh2o)
    i, j = @index(Global, NTuple)
    @inbounds begin
        T = eltype(soilpsi)
        c = cmin + i - 1
        if mask_hydrology[c]
            if h2osoi_liq[c, j + joff] > zero(T)
                vwc = h2osoi_liq[c, j + joff] / (dz[c, j + joff] * denh2o)
                fsattmp = max(vwc / watsat[c, j], T(0.001))
                psi = sucsat[c, j] * T(-9.8e-6) * (fsattmp)^(-bsw[c, j])
                soilpsi[c, j] = min(max(psi, T(-15.0)), zero(T))
            else
                soilpsi[c, j] = -15.0
            end
        end
    end
end
hydrond_soilpsi!(soilpsi, h2osoi_liq, dz, watsat, sucsat, bsw, mask_hydrology,
        joff, bounds, nlevgrnd) =
    _launch!(_hydrond_soilpsi_kernel!, soilpsi, h2osoi_liq, dz, watsat, sucsat, bsw,
             mask_hydrology, joff, first(bounds), eltype(soilpsi)(DENH2O);
             ndrange = (length(bounds), nlevgrnd))

# --- update_smp_l! : 2D (column, layer), each [c,j] independent ---
@kernel function _hydrond_smp_l_kernel!(smp_l, @Const(h2osoi_vol), @Const(watsat),
        @Const(sucsat), @Const(bsw), @Const(smpmin), @Const(mask_hydrology), cmin::Int)
    i, j = @index(Global, NTuple)
    @inbounds begin
        T = eltype(smp_l)
        c = cmin + i - 1
        if mask_hydrology[c]
            s_node = max(h2osoi_vol[c, j] / watsat[c, j], T(0.01))
            s_node = min(one(T), s_node)
            smp_l[c, j] = -sucsat[c, j] * s_node^(-bsw[c, j])
            smp_l[c, j] = max(smpmin[c], smp_l[c, j])
        end
    end
end
hydrond_smp_l!(smp_l, h2osoi_vol, watsat, sucsat, bsw, smpmin, mask_hydrology,
        bounds, nlevgrnd) =
    _launch!(_hydrond_smp_l_kernel!, smp_l, h2osoi_vol, watsat, sucsat, bsw, smpmin,
             mask_hydrology, first(bounds); ndrange = (length(bounds), nlevgrnd))

# =========================================================================
# update_snow_persistence!
# =========================================================================

"""
    update_snow_persistence!(snow_persistence, dtime,
        mask_snow, mask_nosnow, bounds)

For columns with snow, accumulate time-covered-by-snow counter.
For columns without snow, reset counter to zero.

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function update_snow_persistence!(
    snow_persistence::AbstractVector{<:Real},
    dtime::Real,
    mask_snow::AbstractVector{Bool},
    mask_nosnow::AbstractVector{Bool},
    bounds::UnitRange{Int}
)
    hydrond_snow_persistence!(snow_persistence, mask_snow, mask_nosnow, dtime, bounds)
    return nothing
end

# =========================================================================
# accumulate_snow_ice_liq!
# =========================================================================

"""
    accumulate_snow_ice_liq!(snowice, snowliq, h2osoi_ice, h2osoi_liq,
        snl, mask_nolake, mask_snow, bounds, nlevsno)

Vertically sum h2osoi_ice and h2osoi_liq over all snow layers for history output.

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
@kernel function _hydrond_accum_snow_iceliq_kernel!(snowice, snowliq,
        @Const(h2osoi_ice), @Const(h2osoi_liq), @Const(snl),
        @Const(mask_nolake), @Const(mask_snow), nlevsno::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    T = eltype(snowice)
    @inbounds if cmin <= c <= cmax && mask_nolake[c]
        si = zero(T); sl = zero(T)
        if mask_snow[c]
            for j_fortran in (-nlevsno + 1):0
                j = j_fortran + nlevsno
                if j_fortran >= snl[c] + 1
                    si += h2osoi_ice[c, j]
                    sl += h2osoi_liq[c, j]
                end
            end
        end
        snowice[c] = si
        snowliq[c] = sl
    end
end

function accumulate_snow_ice_liq!(
    snowice::AbstractVector{<:Real},
    snowliq::AbstractVector{<:Real},
    h2osoi_ice::AbstractMatrix{<:Real},
    h2osoi_liq::AbstractMatrix{<:Real},
    snl::AbstractVector{<:Integer},
    mask_nolake::AbstractVector{Bool},
    mask_snow::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsno::Int
)
    _launch!(_hydrond_accum_snow_iceliq_kernel!, snowice, snowliq, h2osoi_ice,
             h2osoi_liq, snl, mask_nolake, mask_snow, nlevsno,
             first(bounds), last(bounds); ndrange = last(bounds))
    return nothing
end

# =========================================================================
# update_snowdp!
# =========================================================================

"""
    update_snowdp!(snowdp, snow_depth, frac_sno_eff, bounds)

Calculate column average snow depth = snow_depth * frac_sno_eff.

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function update_snowdp!(
    snowdp::AbstractVector{<:Real},
    snow_depth::AbstractVector{<:Real},
    frac_sno_eff::AbstractVector{<:Real},
    bounds::UnitRange{Int}
)
    hydrond_snowdp!(snowdp, snow_depth, frac_sno_eff, bounds)
    return nothing
end

# =========================================================================
# compute_snow_internal_temperature!
# =========================================================================

"""
    compute_snow_internal_temperature!(t_sno_mul_mss, h2osoi_ice, h2osoi_liq,
        t_soisno, snl, mask_nolake, mask_snow, bounds, nlevsno, tfrz)

Compute the snow internal temperature diagnostic: sum of layer mass × temperature.

The snow internal temperature (SIT) is the weighted average temperature of the
snowpack, weighted by layer mass. This function computes the numerator:
  Sum_i [ (ice_mass(i) * T(i)) + (liq_mass(i) * tfrz) ]

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
# Per-column kernel with an internal sequential snow-layer accumulation. The non-lake
# zero-init and the snow-layer sum are fused: non-snow non-lake columns get 0, snow
# columns get the mass-weighted temperature sum.
@kernel function _hydrond_snow_internal_temp_kernel!(t_sno_mul_mss,
        @Const(h2osoi_ice), @Const(h2osoi_liq), @Const(t_soisno), @Const(snl),
        @Const(mask_nolake), @Const(mask_snow), nlevsno::Int, tfrz, cmin::Int, cmax::Int)
    c = @index(Global)
    T = eltype(t_sno_mul_mss)
    @inbounds if cmin <= c <= cmax && mask_nolake[c]
        v = zero(T)
        if mask_snow[c]
            for j_fortran in (-nlevsno + 1):0
                j = j_fortran + nlevsno
                if j_fortran >= snl[c] + 1
                    v += h2osoi_ice[c, j] * t_soisno[c, j]
                    v += h2osoi_liq[c, j] * tfrz
                end
            end
        end
        t_sno_mul_mss[c] = v
    end
end

function compute_snow_internal_temperature!(
    t_sno_mul_mss::AbstractVector{<:Real},
    h2osoi_ice::AbstractMatrix{<:Real},
    h2osoi_liq::AbstractMatrix{<:Real},
    t_soisno::AbstractMatrix{<:Real},
    snl::AbstractVector{<:Integer},
    mask_nolake::AbstractVector{Bool},
    mask_snow::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsno::Int,
    tfrz::Real
)
    _launch!(_hydrond_snow_internal_temp_kernel!, t_sno_mul_mss, h2osoi_ice,
             h2osoi_liq, t_soisno, snl, mask_nolake, mask_snow, nlevsno,
             eltype(t_sno_mul_mss)(tfrz), first(bounds), last(bounds);
             ndrange = last(bounds))
    return nothing
end

# =========================================================================
# update_ground_and_soil_temperatures!
# =========================================================================

"""
    update_ground_and_soil_temperatures!(t_grnd, t_grnd_u, t_grnd_r,
        tsl, t_soi_10cm, tsoi17,
        t_soisno, t_h2osfc, frac_sno_eff, frac_h2osfc,
        snl, dz, zi, col_landunit, col_itype,
        lun_urbpoi, lun_itype,
        mask_nolake, bounds, nlevsoi, nlevsno)

Compute ground temperature, near-surface soil temperature, and
depth-weighted soil temperatures for top 10cm and 17cm.

- t_grnd: weighted average of exposed soil, snow surface, and surface water
- t_grnd_u: urban ground temperature (= top snow/soil layer)
- t_grnd_r: rural ground temperature (for soil/crop landunits)
- tsl: temperature of near-surface soil layer (layer 1)
- t_soi_10cm: depth-weighted soil temperature in top 10cm
- tsoi17: depth-weighted soil temperature in top 17cm

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
# Per-column kernel fusing the three host loops (init / depth-weighted accumulate /
# finalize+ground-temp). The 10cm and 17cm sums are accumulated into thread-locals
# and written once; ground temperature is computed inline. Non-urban columns get the
# depth-weighted average, urban columns get the top-layer ground temp.
@kernel function _hydrond_ground_soil_temps_kernel!(t_grnd, t_grnd_u, t_grnd_r, tsl,
        t_soi_10cm, tsoi17, @Const(t_soisno), @Const(t_h2osfc), @Const(frac_sno_eff),
        @Const(frac_h2osfc), @Const(snl), @Const(dz), @Const(zi),
        @Const(col_landunit), @Const(col_itype), @Const(lun_urbpoi), @Const(lun_itype),
        @Const(mask_nolake), nlevsoi::Int, joff::Int, joff_zi::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    T = eltype(t_grnd)
    @inbounds if cmin <= c <= cmax && mask_nolake[c]
        l = col_landunit[c]
        if !lun_urbpoi[l]
            acc10 = zero(T)
            acc17 = zero(T)
            for j in 1:nlevsoi
                # Near-surface soil layer temperature
                if j == 1
                    tsl[c] = t_soisno[c, j + joff]
                end

                # Soil T at top 17 cm (F. Li and S. Levis)
                if zi[c, j + joff_zi] <= T(0.17)
                    fracl = one(T)
                    acc17 += t_soisno[c, j + joff] * dz[c, j + joff] * fracl
                else
                    if zi[c, j + joff_zi] > T(0.17) && zi[c, j + joff] < T(0.17)
                        fracl = (T(0.17) - zi[c, j + joff]) / dz[c, j + joff]
                        acc17 += t_soisno[c, j + joff] * dz[c, j + joff] * fracl
                    end
                end

                # Soil T at top 10 cm
                if zi[c, j + joff_zi] <= T(0.1)
                    fracl = one(T)
                    acc10 += t_soisno[c, j + joff] * dz[c, j + joff] * fracl
                else
                    if zi[c, j + joff_zi] > T(0.1) && zi[c, j + joff] < T(0.1)
                        fracl = (T(0.1) - zi[c, j + joff]) / dz[c, j + joff]
                        acc10 += t_soisno[c, j + joff] * dz[c, j + joff] * fracl
                    end
                end
            end
            t_soi_10cm[c] = acc10
            tsoi17[c] = acc17
        end

        # t_grnd is weighted average of exposed soil and snow
        if snl[c] < 0
            t_grnd[c] = frac_sno_eff[c] * t_soisno[c, snl[c] + 1 + joff] +
                         (one(T) - frac_sno_eff[c] - frac_h2osfc[c]) * t_soisno[c, 1 + joff] +
                         frac_h2osfc[c] * t_h2osfc[c]
        else
            t_grnd[c] = (one(T) - frac_h2osfc[c]) * t_soisno[c, 1 + joff] +
                         frac_h2osfc[c] * t_h2osfc[c]
        end

        if lun_urbpoi[l]
            t_grnd_u[c] = t_soisno[c, snl[c] + 1 + joff]
        else
            t_soi_10cm[c] = t_soi_10cm[c] / T(0.1)
            tsoi17[c] = tsoi17[c] / T(0.17)
        end

        if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
            t_grnd_r[c] = t_soisno[c, snl[c] + 1 + joff]
        end
    end
end

function update_ground_and_soil_temperatures!(
    t_grnd::AbstractVector{<:Real},
    t_grnd_u::AbstractVector{<:Real},
    t_grnd_r::AbstractVector{<:Real},
    tsl::AbstractVector{<:Real},
    t_soi_10cm::AbstractVector{<:Real},
    tsoi17::AbstractVector{<:Real},
    t_soisno::AbstractMatrix{<:Real},
    t_h2osfc::AbstractVector{<:Real},
    frac_sno_eff::AbstractVector{<:Real},
    frac_h2osfc::AbstractVector{<:Real},
    snl::AbstractVector{<:Integer},
    dz::AbstractMatrix{<:Real},
    zi::AbstractMatrix{<:Real},
    col_landunit::AbstractVector{<:Integer},
    col_itype::AbstractVector{<:Integer},
    lun_urbpoi::AbstractVector{Bool},
    lun_itype::AbstractVector{<:Integer},
    mask_nolake::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsoi::Int,
    nlevsno::Int
)
    # Julia offset for z/dz/t_soisno: Fortran j → Julia j + nlevsno
    # (arrays dimensioned -nlevsno+1:nlevmaxurbgrnd in Fortran)
    joff = nlevsno
    # Julia offset for zi: Fortran j → Julia j + nlevsno + 1
    # (zi dimensioned -nlevsno:nlevmaxurbgrnd in Fortran, one extra entry)
    joff_zi = nlevsno + 1

    _launch!(_hydrond_ground_soil_temps_kernel!, t_grnd, t_grnd_u, t_grnd_r, tsl,
             t_soi_10cm, tsoi17, t_soisno, t_h2osfc, frac_sno_eff, frac_h2osfc, snl,
             dz, zi, col_landunit, col_itype, lun_urbpoi, lun_itype, mask_nolake,
             nlevsoi, joff, joff_zi, first(bounds), last(bounds);
             ndrange = last(bounds))

    return nothing
end

# =========================================================================
# update_h2osoi_vol!
# =========================================================================

"""
    update_h2osoi_vol!(h2osoi_vol, h2osoi_liq, h2osoi_ice,
        dz, col_itype, mask_nolake, mask_urban, bounds,
        nlevgrnd, nlevurb, nlevsno)

Update volumetric soil water content from liquid and ice masses.

For non-urban non-wall/roof columns: over nlevgrnd layers.
For urban wall/roof columns: over nlevurb layers.

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function update_h2osoi_vol!(
    h2osoi_vol::AbstractMatrix{<:Real},
    h2osoi_liq::AbstractMatrix{<:Real},
    h2osoi_ice::AbstractMatrix{<:Real},
    dz::AbstractMatrix{<:Real},
    col_itype::AbstractVector{<:Integer},
    mask_nolake::AbstractVector{Bool},
    mask_urban::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevgrnd::Int,
    nlevurb::Int,
    nlevsno::Int
)
    joff = nlevsno

    # Non-urban columns (excluding sunwall, shadewall, roof)
    # h2osoi_vol is soil-only (1:nlevmaxurbgrnd), no offset needed
    # h2osoi_liq/ice and dz are snow+soil arrays, need joff offset for soil layer j
    hydrond_h2osoi_vol_nonurban!(h2osoi_vol, h2osoi_liq, h2osoi_ice, dz, col_itype,
                                 mask_nolake, joff, bounds, nlevgrnd)

    # Urban columns (sunwall, shadewall, roof) — only nlevurb layers
    hydrond_h2osoi_vol_urban!(h2osoi_vol, h2osoi_liq, h2osoi_ice, dz, col_itype,
                              mask_urban, joff, bounds, nlevurb)

    return nothing
end

# =========================================================================
# update_soilpsi!
# =========================================================================

"""
    update_soilpsi!(soilpsi, h2osoi_liq, dz, watsat, sucsat, bsw,
        mask_hydrology, bounds, nlevgrnd, nlevsno)

Update soil water potential (soilpsi) in each soil layer [MPa].

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function update_soilpsi!(
    soilpsi::AbstractMatrix{<:Real},
    h2osoi_liq::AbstractMatrix{<:Real},
    dz::AbstractMatrix{<:Real},
    watsat::AbstractMatrix{<:Real},
    sucsat::AbstractMatrix{<:Real},
    bsw::AbstractMatrix{<:Real},
    mask_hydrology::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevgrnd::Int,
    nlevsno::Int
)
    joff = nlevsno

    # h2osoi_liq and dz are snow+soil arrays, need joff offset for soil layer j
    # soilpsi, watsat, sucsat, bsw are soil-only arrays (1-based)
    hydrond_soilpsi!(soilpsi, h2osoi_liq, dz, watsat, sucsat, bsw, mask_hydrology,
                     joff, bounds, nlevgrnd)

    return nothing
end

# =========================================================================
# update_smp_l!
# =========================================================================

"""
    update_smp_l!(smp_l, h2osoi_vol, watsat, sucsat, bsw, smpmin,
        mask_hydrology, bounds, nlevgrnd)

Update soil matric potential (smp_l) for history output and ch4Mod [mm].

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function update_smp_l!(
    smp_l::AbstractMatrix{<:Real},
    h2osoi_vol::AbstractMatrix{<:Real},
    watsat::AbstractMatrix{<:Real},
    sucsat::AbstractMatrix{<:Real},
    bsw::AbstractMatrix{<:Real},
    smpmin::AbstractVector{<:Real},
    mask_hydrology::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevgrnd::Int
)
    hydrond_smp_l!(smp_l, h2osoi_vol, watsat, sucsat, bsw, smpmin, mask_hydrology,
                   bounds, nlevgrnd)

    return nothing
end

# =========================================================================
# compute_wf!
# =========================================================================

"""
    compute_wf!(wf, h2osoi_vol, watsat, sucsat, bsw,
        z, dz, mask_hydrology, bounds, nlevgrnd, nlevsno,
        depth_limit)

Compute soil water as fraction of water holding capacity (WHC)
integrated to a given depth_limit [m].

Available soil water = h2osoi_vol - watdry (air-dry water content)
Potentially available soil water (WHC) = watsat - watdry

Returns wf = available / potentially_available for each column.

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
# Per-column kernel fusing the three host loops; rwat/swat/rz become thread-locals
# (no scratch allocations). The depth-weighted sums are accumulated then the WHC
# fraction finalized inline, replicating the rz==0 fallback exactly.
@kernel function _hydrond_compute_wf_kernel!(wf_out, @Const(h2osoi_vol),
        @Const(watsat), @Const(sucsat), @Const(bsw), @Const(z), @Const(dz),
        @Const(mask_hydrology), nlevgrnd::Int, joff::Int, depth_limit, cmin::Int, cmax::Int)
    c = @index(Global)
    T = eltype(wf_out)
    @inbounds if cmin <= c <= cmax && mask_hydrology[c]
        rwat = zero(T); swat = zero(T); rz = zero(T)
        for j in 1:nlevgrnd
            if z[c, j + joff] + T(0.5) * dz[c, j + joff] <= depth_limit
                watdry = watsat[c, j] * (T(316230.0) / sucsat[c, j])^(T(-1.0) / bsw[c, j])
                rwat += (h2osoi_vol[c, j] - watdry) * dz[c, j + joff]
                swat += (watsat[c, j]      - watdry) * dz[c, j + joff]
                rz   += dz[c, j + joff]
            end
        end
        if rz != zero(T)
            tsw  = rwat / rz
            stsw = swat / rz
        else
            watdry = watsat[c, 1] * (T(316230.0) / sucsat[c, 1])^(T(-1.0) / bsw[c, 1])
            tsw  = h2osoi_vol[c, 1] - watdry
            stsw = watsat[c, 1] - watdry
        end
        wf_out[c] = tsw / stsw
    end
end

function compute_wf!(
    wf_out::AbstractVector{<:Real},
    h2osoi_vol::AbstractMatrix{<:Real},
    watsat::AbstractMatrix{<:Real},
    sucsat::AbstractMatrix{<:Real},
    bsw::AbstractMatrix{<:Real},
    z::AbstractMatrix{<:Real},
    dz::AbstractMatrix{<:Real},
    mask_hydrology::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevgrnd::Int,
    nlevsno::Int,
    depth_limit::Real
)
    joff = nlevsno

    _launch!(_hydrond_compute_wf_kernel!, wf_out, h2osoi_vol, watsat, sucsat, bsw,
             z, dz, mask_hydrology, nlevgrnd, joff, eltype(wf_out)(depth_limit),
             first(bounds), last(bounds); ndrange = last(bounds))

    return nothing
end

# =========================================================================
# update_snow_top_layer_diagnostics!
# =========================================================================

"""
    update_snow_top_layer_diagnostics!(h2osno_top, snw_rds, snot_top,
        dTdz_top, snw_rds_top, sno_liq_top,
        h2osoi_ice, h2osoi_liq, snl,
        mask_snow, mask_nosnow, bounds, nlevsno, spval)

Update top-layer snow diagnostics:
- h2osno_top: mass of snow in top layer
- Zero snow variables and set spval diagnostics for columns without snow.

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
# Per-column kernel fusing the snow / nosnow loops. mask_snow and mask_nosnow are
# disjoint, so each column takes at most one branch. The spval stores are bare
# assignments (no arithmetic), so spval may flow through as a scalar arg.
@kernel function _hydrond_snow_top_diag_kernel!(h2osno_top, snw_rds, snot_top,
        dTdz_top, snw_rds_top, sno_liq_top, @Const(h2osoi_ice), @Const(h2osoi_liq),
        @Const(snl), @Const(mask_snow), @Const(mask_nosnow), nlevsno::Int, spval,
        cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax
        if mask_snow[c]
            j_top = snl[c] + 1 + nlevsno  # Julia index of top snow layer
            h2osno_top[c] = h2osoi_ice[c, j_top] + h2osoi_liq[c, j_top]
        elseif mask_nosnow[c]
            h2osno_top[c] = 0.0
            for j in 1:nlevsno
                snw_rds[c, j] = 0.0
            end
            snot_top[c]    = spval
            dTdz_top[c]    = spval
            snw_rds_top[c] = spval
            sno_liq_top[c] = spval
        end
    end
end

function update_snow_top_layer_diagnostics!(
    h2osno_top::AbstractVector{<:Real},
    snw_rds::AbstractMatrix{<:Real},
    snot_top::AbstractVector{<:Real},
    dTdz_top::AbstractVector{<:Real},
    snw_rds_top::AbstractVector{<:Real},
    sno_liq_top::AbstractVector{<:Real},
    h2osoi_ice::AbstractMatrix{<:Real},
    h2osoi_liq::AbstractMatrix{<:Real},
    snl::AbstractVector{<:Integer},
    mask_snow::AbstractVector{Bool},
    mask_nosnow::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsno::Int,
    spval::Real
)
    _launch!(_hydrond_snow_top_diag_kernel!, h2osno_top, snw_rds, snot_top, dTdz_top,
             snw_rds_top, sno_liq_top, h2osoi_ice, h2osoi_liq, snl, mask_snow,
             mask_nosnow, nlevsno, eltype(snot_top)(spval),
             first(bounds), last(bounds); ndrange = last(bounds))

    return nothing
end

# =========================================================================
# calc_and_withdraw_irrigation_fluxes!
# =========================================================================

"""
    calc_and_withdraw_irrigation_fluxes!(soilhydrology, soilstate,
        waterfluxbulk, waterstatebulk, water,
        mask_soil, bounds, nlevsno, nlevgrnd;
        use_groundwater_irrigation=false)

Calculates irrigation withdrawal fluxes and withdraws from groundwater.

In Fortran this calls `irrigation_inst%CalcIrrigationFluxes` and
`WithdrawGroundwaterIrrigation`. Since IrrigationMod is not yet ported,
this is a stub that only performs the groundwater withdrawal step when
`use_groundwater_irrigation` is true.

Ported from `CalcAndWithdrawIrrigationFluxes` in `HydrologyNoDrainageMod.F90`.
"""
function calc_and_withdraw_irrigation_fluxes!(
    waterflux::WaterFluxData,
    waterstate::WaterStateData,
    mask_soil::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsoi::Int,
    dtime::Real;
    use_groundwater_irrigation::Bool = false
)
    # Stub: CalcIrrigationFluxes requires IrrigationMod (not yet ported)
    # When IrrigationMod is ported, add the irrigation flux calculation here.

    # Groundwater irrigation withdrawal
    if use_groundwater_irrigation
        withdraw_groundwater_irrigation!(
            waterflux, waterstate,
            mask_soil, bounds, nlevsoi, dtime)
    end

    return nothing
end

# =========================================================================
# handle_new_snow!
# =========================================================================

"""
    handle_new_snow!(temperature, waterstatebulk, waterdiagbulk,
        col, lun, mask_nolake, bounds, dtime, nlevsno;
        forc_t, qflx_snow_grnd)

Handle new snow falling on the ground.

Coordinates:
1. Update quantities for new snow (bulk density, diagnostics, state)
2. Remove snow from thawed wetlands
3. Initialize explicit snow pack where warranted

Ported from `HandleNewSnow` in `HydrologyNoDrainageMod.F90`.
"""
# Fallback temperature-only new-snow bulk density (used only when no wind forcing
# is supplied). Mirrors the scalar branches; literals carried at the working
# precision so no Float64 reaches a Float32-only backend (Metal). Byte-identical on
# Float64.
@kernel function _hydrond_bifall_fallback_kernel!(bifall, @Const(mask_nolake),
        @Const(forc_t), tfrz, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_nolake[c]
        T = eltype(bifall)
        ft = forc_t[c]
        if ft > tfrz + T(2.0)
            bifall[c] = T(50.0) + T(1.7) * T(17.0)^T(1.5)
        elseif ft > tfrz - T(15.0)
            bifall[c] = T(50.0) + T(1.7) * (ft - tfrz + T(15.0))^T(1.5)
        else
            bifall[c] = T(50.0)
        end
    end
end

# Total H2O in the snow pack per column: the no-layer reservoir plus the sum of
# ice+liquid over the resolved snow layers [snl+1 .. 0]. This is a per-column kernel
# with an INTERNAL sequential snow-layer loop (loop-carried accumulation), the
# established pattern for variable-snl reductions.
@kernel function _hydrond_h2osno_total_kernel!(h2osno_total, @Const(mask_nolake),
        @Const(snl), @Const(h2osno_no_layers), @Const(h2osoi_ice), @Const(h2osoi_liq),
        nlevsno::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_nolake[c]
        tot = h2osno_no_layers[c]
        for j in (snl[c] + 1):0
            jj = j + nlevsno
            tot += h2osoi_ice[c, jj] + h2osoi_liq[c, jj]
        end
        h2osno_total[c] = tot
    end
end

# Gather per-column landunit type and urban flag (landunit-indexed reads) so the
# downstream diagnostics see flat column-indexed inputs.
@kernel function _hydrond_gather_lun_kernel!(lun_itype_col, @Const(mask_nolake),
        @Const(col_landunit), @Const(lun_itype), @Const(lun_urbpoi), urbpoi_col,
        cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_nolake[c]
        l = col_landunit[c]
        lun_itype_col[c] = lun_itype[l]
        urbpoi_col[c] = lun_urbpoi[l]
    end
end

function handle_new_snow!(
    temperature::TemperatureData,
    waterstatebulk::WaterStateBulkData,
    waterdiagbulk::WaterDiagnosticBulkData,
    col::ColumnData,
    lun::LandunitData,
    mask_nolake::AbstractVector{Bool},
    bounds::UnitRange{Int},
    dtime::Real,
    nlevsno::Int;
    forc_t::AbstractVector{<:Real},
    forc_wind::AbstractVector{<:Real} = Float64[],
    qflx_snow_grnd::AbstractVector{<:Real},
    qflx_snow_drain::AbstractVector{<:Real} = zeros(length(mask_nolake)),
    int_snow::AbstractVector{<:Real} = zeros(length(mask_nolake)),
    scf_method::SnowCoverFractionBase = SnowCoverFractionSwensonLawrence2012()
)
    nc = length(mask_nolake)
    FT = eltype(forc_t)

    # --- 1. Update quantities for new snow ---
    # Step 1a: Compute bulk density of new snow
    bifall = fill!(similar(forc_t, FT, nc), FT(50.0))
    if !isempty(forc_wind)
        new_snow_bulk_density!(bifall, forc_t, forc_wind, col.gridcell,
                               mask_nolake, bounds)
    else
        # Fallback: use temperature-only bulk density
        _launch!(_hydrond_bifall_fallback_kernel!, bifall, mask_nolake, forc_t,
                 FT(TFRZ), first(bounds), last(bounds))
    end

    # Step 1b: Calculate total H2O in snow (h2osno_total) — per-column kernel with
    # an internal sequential snow-layer loop (loop-carried accumulation).
    FT_hs = eltype(forc_t)
    h2osno_total = fill!(similar(forc_t, FT_hs, nc), zero(FT_hs))
    _launch!(_hydrond_h2osno_total_kernel!, h2osno_total, mask_nolake, col.snl,
             waterstatebulk.ws.h2osno_no_layers_col,
             waterstatebulk.ws.h2osoi_ice_col, waterstatebulk.ws.h2osoi_liq_col,
             nlevsno, first(bounds), last(bounds))

    # Step 1c: Update snow diagnostics (snow_depth, frac_sno, int_snow, swe_old)
    lun_itype_col = similar(col.landunit, Int, nc)
    urbpoi_col = similar(col.landunit, Bool, nc)
    _launch!(_hydrond_gather_lun_kernel!, lun_itype_col, mask_nolake, col.landunit,
             lun.itype, lun.urbpoi, urbpoi_col, first(bounds), last(bounds))

    bulkdiag_new_snow_diagnostics!(
        col.dz, int_snow, waterdiagbulk.swe_old_col,
        waterdiagbulk.frac_sno_col, waterdiagbulk.frac_sno_eff_col,
        waterdiagbulk.snow_depth_col, waterdiagbulk.snomelt_accum_col,
        dtime, lun_itype_col, urbpoi_col, col.snl, bifall,
        h2osno_total,
        waterstatebulk.ws.h2osoi_ice_col, waterstatebulk.ws.h2osoi_liq_col,
        qflx_snow_grnd, qflx_snow_drain,
        mask_nolake, bounds, nlevsno;
        scf_method=scf_method)

    # Step 1d: Add new snow mass to appropriate state variable
    update_state_add_new_snow!(
        waterstatebulk.ws.h2osno_no_layers_col,
        waterstatebulk.ws.h2osoi_ice_col,
        dtime, col.snl, qflx_snow_grnd,
        mask_nolake, bounds, nlevsno)

    # --- 2. Remove snow from thawed wetlands ---
    # Reuse lun_itype_col from Step 1c (already populated for all active columns)

    mask_thawed_wetland = fill!(similar(mask_nolake, Bool, nc), false)
    build_filter_thawed_wetland_thin_snowpack!(
        mask_thawed_wetland,
        temperature.t_grnd_col, lun_itype_col, col.snl,
        mask_nolake, bounds)

    update_state_remove_snow_thawed_wetlands!(
        waterstatebulk.ws.h2osno_no_layers_col,
        mask_thawed_wetland, bounds)

    bulk_remove_snow_thawed_wetlands!(
        waterdiagbulk.snow_depth_col,
        mask_thawed_wetland, bounds)

    # --- 3. Initialize explicit snow pack ---
    mask_init_snowpack = fill!(similar(mask_nolake, Bool, nc), false)
    build_filter_snowpack_initialized!(
        mask_init_snowpack,
        col.snl, lun_itype_col,
        waterdiagbulk.frac_sno_eff_col,
        waterdiagbulk.snow_depth_col,
        qflx_snow_grnd,
        mask_nolake, bounds)

    update_state_initialize_snow_pack!(
        waterstatebulk.ws.h2osno_no_layers_col,
        waterstatebulk.ws.h2osoi_ice_col,
        waterstatebulk.ws.h2osoi_liq_col,
        mask_init_snowpack, bounds, nlevsno)

    bulk_initialize_snow_pack!(
        col.snl, col.zi, col.dz, col.z,
        temperature.t_soisno_col,
        waterdiagbulk.frac_iceold_col,
        waterdiagbulk.snomelt_accum_col,
        forc_t, waterdiagbulk.snow_depth_col,
        mask_init_snowpack, bounds, nlevsno)

    return nothing
end

# =========================================================================
# hydrology_no_drainage! — Main orchestrator
# =========================================================================

"""
    hydrology_no_drainage!(temperature, soilstate, waterstatebulk, waterfluxbulk,
        waterdiagbulk, soilhydrology, col, lun, water,
        mask_nolake, mask_hydrology, mask_urban, mask_snow, mask_nosnow,
        bounds, dtime, nlevsno, nlevsoi, nlevgrnd, nlevurb)

Main subroutine to execute the calculation of soil/snow hydrology
without drainage.

This is the top-level orchestrator that calls sub-functions for:
  - Snow persistence tracking
  - Snow ice/liquid accumulation
  - Ground temperature calculation
  - Volumetric soil water update
  - Soil water potential update
  - Soil water fraction diagnostics
  - Top-layer snow diagnostics

Note: Many sub-calls from the Fortran original (SnowWater, SoilWater,
WaterTable, etc.) are called separately in the driver. This function
focuses on the inline diagnostic/state updates that follow.

Ported from `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function hydrology_no_drainage!(
    temperature::TemperatureData,
    soilstate::SoilStateData,
    waterstatebulk::WaterStateBulkData,
    waterdiagbulk::WaterDiagnosticBulkData,
    col::ColumnData,
    lun::LandunitData,
    mask_nolake::AbstractVector{Bool},
    mask_hydrology::AbstractVector{Bool},
    mask_urban::AbstractVector{Bool},
    mask_snow::AbstractVector{Bool},
    mask_nosnow::AbstractVector{Bool},
    bounds::UnitRange{Int},
    dtime::Real,
    nlevsno::Int,
    nlevsoi::Int,
    nlevgrnd::Int,
    nlevurb::Int;
    soilhydrology::Union{SoilHydrologyData, Nothing} = nothing
)
    # Note: set_soil_water_fractions! is now called in the driver before physics

    # --- Snow persistence ---
    update_snow_persistence!(
        waterstatebulk.snow_persistence_col,
        dtime, mask_snow, mask_nosnow, bounds)

    # --- Accumulate snowice and snowliq over snow layers ---
    accumulate_snow_ice_liq!(
        waterdiagbulk.snowice_col,
        waterdiagbulk.snowliq_col,
        waterstatebulk.ws.h2osoi_ice_col,
        waterstatebulk.ws.h2osoi_liq_col,
        col.snl, mask_nolake, mask_snow, bounds, nlevsno)

    # --- Column-average snow depth ---
    update_snowdp!(
        waterdiagbulk.snowdp_col,
        waterdiagbulk.snow_depth_col,
        waterdiagbulk.frac_sno_eff_col,
        bounds)

    # --- Snow internal temperature ---
    compute_snow_internal_temperature!(
        temperature.t_sno_mul_mss_col,
        waterstatebulk.ws.h2osoi_ice_col,
        waterstatebulk.ws.h2osoi_liq_col,
        temperature.t_soisno_col,
        col.snl, mask_nolake, mask_snow, bounds, nlevsno, TFRZ)

    # --- Ground and soil temperatures ---
    update_ground_and_soil_temperatures!(
        temperature.t_grnd_col,
        temperature.t_grnd_u_col,
        temperature.t_grnd_r_col,
        temperature.tsl_col,
        temperature.t_soi10cm_col,
        temperature.t_soi17cm_col,
        temperature.t_soisno_col,
        temperature.t_h2osfc_col,
        waterdiagbulk.frac_sno_eff_col,
        waterdiagbulk.frac_h2osfc_col,
        col.snl, col.dz, col.zi,
        col.landunit, col.itype,
        lun.urbpoi, lun.itype,
        mask_nolake, bounds, nlevsoi, nlevsno)

    # --- Volumetric soil water ---
    update_h2osoi_vol!(
        waterstatebulk.ws.h2osoi_vol_col,
        waterstatebulk.ws.h2osoi_liq_col,
        waterstatebulk.ws.h2osoi_ice_col,
        col.dz, col.itype,
        mask_nolake, mask_urban, bounds,
        nlevgrnd, nlevurb, nlevsno)

    # --- Soil water potential ---
    update_soilpsi!(
        soilstate.soilpsi_col,
        waterstatebulk.ws.h2osoi_liq_col,
        col.dz,
        soilstate.watsat_col,
        soilstate.sucsat_col,
        soilstate.bsw_col,
        mask_hydrology, bounds, nlevgrnd, nlevsno)

    # --- Matric potential for history/ch4 ---
    update_smp_l!(
        soilstate.smp_l_col,
        waterstatebulk.ws.h2osoi_vol_col,
        soilstate.watsat_col,
        soilstate.sucsat_col,
        soilstate.bsw_col,
        soilstate.smpmin_col,
        mask_hydrology, bounds, nlevgrnd)

    # --- Soil water as fraction of WHC (top 0.05m) ---
    compute_wf!(
        waterdiagbulk.wf_col,
        waterstatebulk.ws.h2osoi_vol_col,
        soilstate.watsat_col,
        soilstate.sucsat_col,
        soilstate.bsw_col,
        col.z, col.dz,
        mask_hydrology, bounds, nlevgrnd, nlevsno, 0.05)

    # --- Soil water as fraction of WHC (top 0.17m) ---
    compute_wf!(
        waterdiagbulk.wf2_col,
        waterstatebulk.ws.h2osoi_vol_col,
        soilstate.watsat_col,
        soilstate.sucsat_col,
        soilstate.bsw_col,
        col.z, col.dz,
        mask_hydrology, bounds, nlevgrnd, nlevsno, 0.17)

    # --- Top-layer snow diagnostics ---
    update_snow_top_layer_diagnostics!(
        waterdiagbulk.h2osno_top_col,
        waterdiagbulk.snw_rds_col,
        temperature.snot_top_col,
        temperature.dTdz_top_col,
        waterdiagbulk.snw_rds_top_col,
        waterdiagbulk.sno_liq_top_col,
        waterstatebulk.ws.h2osoi_ice_col,
        waterstatebulk.ws.h2osoi_liq_col,
        col.snl,
        mask_snow, mask_nosnow, bounds, nlevsno, SPVAL)

    return nothing
end
